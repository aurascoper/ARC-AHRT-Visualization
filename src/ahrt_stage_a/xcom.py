"""NIST XCOM data import and log-log interpolation.

The application deliberately contains no physical XCOM coefficient constants.
A normalized CSV export and provenance manifest are required at runtime.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REQUIRED_COLUMNS = (
    "energy_MeV",
    "mu_coherent_cm2_g",
    "mu_incoherent_cm2_g",
    "mu_photoelectric_cm2_g",
    "mu_pair_nuclear_cm2_g",
    "mu_pair_electron_cm2_g",
    "mu_total_with_coherent_cm2_g",
    "mu_total_without_coherent_cm2_g",
)


class XCOMDataError(ValueError):
    """Raised when XCOM data or its provenance is incomplete or inconsistent."""


@dataclass(frozen=True)
class XCOMCoefficients:
    energy_MeV: float
    coherent: float
    incoherent: float
    photoelectric: float
    pair_nuclear: float
    pair_electron: float
    total_with_coherent: float
    total_without_coherent: float

    @property
    def pair_total(self) -> float:
        return self.pair_nuclear + self.pair_electron

    @property
    def partial_sum_without_coherent(self) -> float:
        return self.incoherent + self.photoelectric + self.pair_total


class XCOMTable:
    def __init__(self, rows: Iterable[XCOMCoefficients], source_path: Path):
        self.rows = sorted(rows, key=lambda r: r.energy_MeV)
        self.source_path = source_path
        if len(self.rows) < 2:
            raise XCOMDataError("XCOM table needs at least two energy rows.")
        energies = [r.energy_MeV for r in self.rows]
        if any(not math.isfinite(e) or e <= 0 for e in energies):
            raise XCOMDataError("All XCOM energies must be finite and positive.")
        if any(b <= a for a, b in zip(energies, energies[1:])):
            raise XCOMDataError("XCOM energies must be unique and strictly increasing.")

    @classmethod
    def from_csv(cls, path: str | Path) -> "XCOMTable":
        file_path = Path(path)
        if not file_path.is_file():
            raise FileNotFoundError(file_path)
        rows: list[XCOMCoefficients] = []
        with file_path.open("r", newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise XCOMDataError("XCOM CSV has no header row.")
            missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
            if missing:
                raise XCOMDataError(f"XCOM CSV is missing columns: {', '.join(missing)}")
            for line_no, raw in enumerate(reader, start=2):
                if not raw or all(v is None or not str(v).strip() for v in raw.values()):
                    continue
                try:
                    vals = {k: float(raw[k]) for k in REQUIRED_COLUMNS}
                except (TypeError, ValueError) as exc:
                    raise XCOMDataError(f"Invalid numeric value at CSV line {line_no}.") from exc
                for key, value in vals.items():
                    if not math.isfinite(value) or (key != "energy_MeV" and value < 0):
                        raise XCOMDataError(f"Invalid value for {key} at CSV line {line_no}.")
                rows.append(
                    XCOMCoefficients(
                        energy_MeV=vals["energy_MeV"],
                        coherent=vals["mu_coherent_cm2_g"],
                        incoherent=vals["mu_incoherent_cm2_g"],
                        photoelectric=vals["mu_photoelectric_cm2_g"],
                        pair_nuclear=vals["mu_pair_nuclear_cm2_g"],
                        pair_electron=vals["mu_pair_electron_cm2_g"],
                        total_with_coherent=vals["mu_total_with_coherent_cm2_g"],
                        total_without_coherent=vals["mu_total_without_coherent_cm2_g"],
                    )
                )
        return cls(rows, file_path)

    def at(self, energy_MeV: float) -> XCOMCoefficients:
        """Interpolate coefficients in log-log space between adjacent rows.

        XCOM uses coefficient-specific interpolation procedures. This compact
        implementation uses log-log linear interpolation as a transparent
        first-stage approximation; it should be benchmarked against a direct
        XCOM query at the energy of interest. It does not interpolate across
        identified absorption edges.
        """
        if not math.isfinite(energy_MeV) or energy_MeV <= 0:
            raise ValueError("Requested energy must be finite and positive.")
        if energy_MeV < self.rows[0].energy_MeV or energy_MeV > self.rows[-1].energy_MeV:
            raise XCOMDataError(
                f"Energy {energy_MeV:g} MeV outside XCOM table range "
                f"[{self.rows[0].energy_MeV:g}, {self.rows[-1].energy_MeV:g}] MeV."
            )

        for row in self.rows:
            if math.isclose(row.energy_MeV, energy_MeV, rel_tol=1e-12, abs_tol=0.0):
                return row

        lo = self.rows[0]
        hi = self.rows[-1]
        for a, b in zip(self.rows, self.rows[1:]):
            if a.energy_MeV < energy_MeV < b.energy_MeV:
                lo, hi = a, b
                break

        fraction = (math.log(energy_MeV) - math.log(lo.energy_MeV)) / (
            math.log(hi.energy_MeV) - math.log(lo.energy_MeV)
        )

        def interp(v0: float, v1: float) -> float:
            if v0 == 0.0 and v1 == 0.0:
                return 0.0
            if v0 <= 0.0 or v1 <= 0.0:
                # Do not invent a power law through a process threshold.
                return max(0.0, v0 + fraction * (v1 - v0))
            return math.exp(math.log(v0) + fraction * (math.log(v1) - math.log(v0)))

        names = (
            "coherent", "incoherent", "photoelectric", "pair_nuclear",
            "pair_electron",
        )
        values = {name: interp(getattr(lo, name), getattr(hi, name)) for name in names}
        # Reconstruct totals from interpolated partials so process selection
        # and free-path sampling share exactly the same normalized coefficient.
        partial_without = values["incoherent"] + values["photoelectric"] + values["pair_nuclear"] + values["pair_electron"]
        values["total_without_coherent"] = partial_without
        values["total_with_coherent"] = partial_without + values["coherent"]
        return XCOMCoefficients(energy_MeV=energy_MeV, **values)

    def max_row_closure_errors(self) -> dict:
        """Check raw imported rows, before any interpolation can mask errors."""
        worst_without = (0.0, None)
        worst_with = (0.0, None)
        for row in self.rows:
            err_without, err_with = coefficient_closure_error(row)
            if err_without > worst_without[0]:
                worst_without = (err_without, row.energy_MeV)
            if err_with > worst_with[0]:
                worst_with = (err_with, row.energy_MeV)
        return {
            "max_relative_error_without_coherent": worst_without[0],
            "worst_energy_without_coherent_MeV": worst_without[1],
            "max_relative_error_with_coherent": worst_with[0],
            "worst_energy_with_coherent_MeV": worst_with[1],
        }


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_and_verify_manifest(table_path: str | Path, manifest_path: str | Path) -> dict:
    table_path = Path(table_path)
    manifest_path = Path(manifest_path)
    if not manifest_path.is_file():
        raise XCOMDataError(f"Missing XCOM provenance manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    required = ("source", "material", "units", "retrieved_utc", "sha256")
    missing = [k for k in required if not manifest.get(k)]
    if missing:
        raise XCOMDataError(f"Manifest missing fields: {', '.join(missing)}")
    if "nist.gov" not in manifest["source"].lower():
        raise XCOMDataError("XCOM source must identify an official NIST URL.")
    if manifest["units"] != "cm2/g; energy_MeV":
        raise XCOMDataError("Unexpected units; expected 'cm2/g; energy_MeV'.")
    actual = sha256_file(table_path)
    if actual.lower() != str(manifest["sha256"]).lower():
        raise XCOMDataError("XCOM CSV SHA-256 does not match its manifest.")
    return manifest


def coefficient_closure_error(coeff: XCOMCoefficients) -> tuple[float, float]:
    """Return relative closure errors for totals with/without coherent scatter."""
    total_no = coeff.total_without_coherent
    total_yes = coeff.total_with_coherent
    partial_no = coeff.partial_sum_without_coherent
    partial_yes = partial_no + coeff.coherent
    err_no = abs(total_no - partial_no) / max(abs(total_no), 1e-300)
    err_yes = abs(total_yes - partial_yes) / max(abs(total_yes), 1e-300)
    return err_no, err_yes
