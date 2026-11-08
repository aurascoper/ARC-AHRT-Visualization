"""Stage A Monte Carlo photon transport in a homogeneous water cylinder.

No electron transport, track-structure physics, Manning/PB calculation, or
biological response is implemented here. Compton recoil electrons are emitted
as weighted phase-space source records for a later transport stage.
"""

from __future__ import annotations

import csv
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path

from .geometry import Cylinder
from .klein_nishina import (
    electron_direction_from_compton,
    rotate_direction,
    sample_compton,
)
from .xcom import XCOMTable


@dataclass(frozen=True)
class TransportConfig:
    photon_energy_MeV: float = 0.661657
    density_g_cm3: float = 1.0
    macro_radius_cm: float = 0.3
    macro_length_cm: float = 0.6
    cell_radius_cm: float = 5.0e-4
    cell_length_cm: float = 2.0e-3
    photon_cutoff_MeV: float = 0.001
    max_photon_interactions: int = 100
    seed: int = 20261103
    histories: int = 100_000

    def validate(self) -> None:
        numeric_positive = {
            "photon_energy_MeV": self.photon_energy_MeV,
            "density_g_cm3": self.density_g_cm3,
            "macro_radius_cm": self.macro_radius_cm,
            "macro_length_cm": self.macro_length_cm,
            "cell_radius_cm": self.cell_radius_cm,
            "cell_length_cm": self.cell_length_cm,
            "photon_cutoff_MeV": self.photon_cutoff_MeV,
        }
        for name, value in numeric_positive.items():
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive.")
        if self.cell_radius_cm >= self.macro_radius_cm or self.cell_length_cm >= self.macro_length_cm:
            raise ValueError("Cell cylinder must be smaller than macro cylinder.")
        if self.histories <= 0 or self.max_photon_interactions <= 0:
            raise ValueError("histories and max_photon_interactions must be positive.")
        if self.photon_energy_MeV >= 1.022 and self.photon_energy_MeV > 0:
            # Pair production is outside this Stage A event generator.
            raise ValueError(
                "Stage A v0.1 is scoped to the 661.657 keV Cs-137 line; "
                "pair-production-capable source spectra require a later release."
            )


@dataclass
class RunTally:
    histories: int = 0
    photon_interactions: int = 0
    compton_interactions: int = 0
    photoelectric_interactions: int = 0
    photon_escapes: int = 0
    cutoff_terminations: int = 0
    max_interaction_cap_terminations: int = 0
    direct_cell_photon_interactions: int = 0
    secondary_electron_sources: int = 0
    secondary_electron_energy_MeV: float = 0.0
    photoelectric_local_energy_MeV: float = 0.0
    escaped_photon_energy_MeV: float = 0.0
    below_cutoff_unresolved_energy_MeV: float = 0.0
    capped_unresolved_photon_energy_MeV: float = 0.0


class PhotonTransport:
    def __init__(self, table: XCOMTable, config: TransportConfig):
        config.validate()
        raw_closure = table.max_row_closure_errors()
        if max(
            raw_closure["max_relative_error_without_coherent"],
            raw_closure["max_relative_error_with_coherent"],
        ) > 0.01:
            raise ValueError("XCOM table fails the 1% raw-row coefficient closure gate; run validate-xcom first.")
        self.table = table
        self.config = config
        self.domain = Cylinder(config.macro_radius_cm, config.macro_length_cm)
        self.cell = Cylinder(config.cell_radius_cm, config.cell_length_cm)

    def run(self, output_dir: str | Path) -> dict:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        rng = random.Random(self.config.seed)
        tally = RunTally(histories=self.config.histories)
        phase_path = out / "compton_electron_sources.csv"
        interaction_path = out / "photon_interactions.csv"

        with phase_path.open("w", newline="", encoding="utf-8") as phase_stream, interaction_path.open(
            "w", newline="", encoding="utf-8"
        ) as interaction_stream:
            phase_writer = csv.DictWriter(
                phase_stream,
                fieldnames=(
                    "history_id", "interaction_id", "x_cm", "y_cm", "z_cm",
                    "electron_kinetic_energy_MeV", "electron_dx", "electron_dy", "electron_dz",
                    "parent_photon_energy_MeV", "statistical_weight",
                ),
            )
            interaction_writer = csv.DictWriter(
                interaction_stream,
                fieldnames=(
                    "history_id", "interaction_id", "process", "x_cm", "y_cm", "z_cm",
                    "photon_energy_in_MeV", "photon_energy_out_MeV", "inside_cell_roi",
                ),
            )
            phase_writer.writeheader()
            interaction_writer.writeheader()

            for history_id in range(self.config.histories):
                p = self.domain.sample_parallel_beam_entry(rng)
                d = (1.0, 0.0, 0.0)
                energy = self.config.photon_energy_MeV
                for interaction_id in range(self.config.max_photon_interactions):
                    coeff = self.table.at(energy)
                    # NIST XCOM total attenuation without coherent scattering is
                    # the transport coefficient used in this first-order model.
                    mu_mass = coeff.total_without_coherent
                    if mu_mass <= 0:
                        raise RuntimeError(f"Non-positive attenuation at {energy:g} MeV.")
                    mu_linear = mu_mass * self.config.density_g_cm3  # cm^-1
                    free_path = -math.log(max(rng.random(), 1e-300)) / mu_linear
                    boundary_distance = self.domain.distance_to_exit(p, d)
                    if free_path >= boundary_distance:
                        tally.photon_escapes += 1
                        tally.escaped_photon_energy_MeV += energy
                        break

                    p = tuple(p[i] + free_path * d[i] for i in range(3))
                    tally.photon_interactions += 1
                    inside_cell = self.cell.contains(p)
                    if inside_cell:
                        tally.direct_cell_photon_interactions += 1

                    # Coherent (Rayleigh) scattering is intentionally excluded
                    # by using XCOM's total attenuation without coherent scatter.
                    comps = (
                        ("compton", coeff.incoherent),
                        ("photoelectric", coeff.photoelectric),
                        ("pair", coeff.pair_total),
                    )
                    total = sum(v for _, v in comps)
                    draw = rng.random() * total
                    cumulative = 0.0
                    process = "compton"
                    for label, value in comps:
                        cumulative += value
                        if draw <= cumulative:
                            process = label
                            break

                    if process == "pair":
                        raise RuntimeError("Pair production encountered; unsupported in AHRT Stage A v0.1.")
                    if process == "photoelectric":
                        tally.photoelectric_interactions += 1
                        tally.photoelectric_local_energy_MeV += energy
                        interaction_writer.writerow({
                            "history_id": history_id, "interaction_id": interaction_id,
                            "process": process, "x_cm": p[0], "y_cm": p[1], "z_cm": p[2],
                            "photon_energy_in_MeV": energy, "photon_energy_out_MeV": 0.0,
                            "inside_cell_roi": inside_cell,
                        })
                        break

                    # Free-electron KN angular law. XCOM's integrated incoherent
                    # coefficient includes bound-electron scattering corrections;
                    # that mismatch is a declared v0.1 approximation.
                    tally.compton_interactions += 1
                    sample = sample_compton(energy, rng)
                    d_out = rotate_direction(d, sample.cos_theta, sample.phi_rad)
                    e_dir = electron_direction_from_compton(
                        d, d_out, energy, sample.photon_energy_out_MeV
                    )
                    e_kin = sample.electron_kinetic_energy_MeV
                    tally.secondary_electron_sources += 1
                    tally.secondary_electron_energy_MeV += e_kin
                    phase_writer.writerow({
                        "history_id": history_id,
                        "interaction_id": interaction_id,
                        "x_cm": p[0], "y_cm": p[1], "z_cm": p[2],
                        "electron_kinetic_energy_MeV": e_kin,
                        "electron_dx": e_dir[0], "electron_dy": e_dir[1], "electron_dz": e_dir[2],
                        "parent_photon_energy_MeV": energy,
                        "statistical_weight": 1.0,
                    })
                    interaction_writer.writerow({
                        "history_id": history_id, "interaction_id": interaction_id,
                        "process": process, "x_cm": p[0], "y_cm": p[1], "z_cm": p[2],
                        "photon_energy_in_MeV": energy,
                        "photon_energy_out_MeV": sample.photon_energy_out_MeV,
                        "inside_cell_roi": inside_cell,
                    })
                    energy = sample.photon_energy_out_MeV
                    d = d_out
                    if energy < self.config.photon_cutoff_MeV:
                        tally.cutoff_terminations += 1
                        tally.below_cutoff_unresolved_energy_MeV += energy
                        break
                else:
                    tally.max_interaction_cap_terminations += 1
                    tally.capped_unresolved_photon_energy_MeV += energy

        summary = {
            "model": "AHRT Stage A v0.1",
            "scope": "primary photon transport and Compton recoil-electron source phase space only",
            "physics_omissions": [
                "electron transport and track structure",
                "coherent Rayleigh scattering",
                "photoelectron and atomic-relaxation cascades",
                "pair production",
                "water radiolysis and chemical transport",
                "Manning/PB and ion binding",
                "Arc-specific molecular damage or biological response",
            ],
            "config": asdict(self.config),
            "tallies": asdict(tally),
            "energy_accounting": {
                "incident_energy_MeV": self.config.histories * self.config.photon_energy_MeV,
                "escaped_photon_energy_MeV": tally.escaped_photon_energy_MeV,
                "compton_electron_source_energy_MeV": tally.secondary_electron_energy_MeV,
                "photoelectric_local_energy_MeV": tally.photoelectric_local_energy_MeV,
                "below_cutoff_unresolved_energy_MeV": tally.below_cutoff_unresolved_energy_MeV,
                "interaction_cap_unresolved_energy_MeV": tally.capped_unresolved_photon_energy_MeV,
                "note": "Compton recoil-electron energy is emitted as phase-space source energy, not transported or deposited dose.",
            },
        }
        (out / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return summary
