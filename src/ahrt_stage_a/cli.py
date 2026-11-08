"""Command line interface for AHRT Stage A."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .transport import PhotonTransport, TransportConfig
from .validation import validate_kn, validate_xcom
from .xcom import XCOMDataError, XCOMTable, load_and_verify_manifest


def _write_json(path: str | Path, obj: dict) -> None:
    Path(path).write_text(json.dumps(obj, indent=2), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ahrt", description="AHRT Stage A reference prototype")
    sub = parser.add_subparsers(dest="command", required=True)

    xcom = sub.add_parser("validate-xcom", help="Validate NIST XCOM import and coefficient closure")
    xcom.add_argument("--table", required=True)
    xcom.add_argument("--manifest", required=True)
    xcom.add_argument("--energy-MeV", type=float, default=0.661657)
    xcom.add_argument("--density-g-cm3", type=float, default=1.0)
    xcom.add_argument("--slab-thickness-cm", type=float, default=1.0)
    xcom.add_argument("--histories", type=int, default=3_000_000)
    xcom.add_argument("--seed", type=int, default=445566)
    xcom.add_argument("--out", default="results/g1_xcom.json")

    kn = sub.add_parser("validate-kn", help="Validate sampled Compton energy-angle law")
    kn.add_argument("--energy-MeV", type=float, default=0.661657)
    kn.add_argument("--samples", type=int, default=150_000)
    kn.add_argument("--seed", type=int, default=81723)
    kn.add_argument("--out", default="results/g1_kn.json")

    run = sub.add_parser("run", help="Run Stage A primary photon transport")
    run.add_argument("--table", required=True)
    run.add_argument("--manifest", required=True)
    run.add_argument("--histories", type=int, default=100_000)
    run.add_argument("--seed", type=int, default=20261103)
    run.add_argument("--energy-MeV", type=float, default=0.661657)
    run.add_argument("--density-g-cm3", type=float, default=1.0)
    run.add_argument("--outdir", default="results/cs137_reference")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "validate-xcom":
            result = validate_xcom(
                args.table, args.manifest, args.energy_MeV, args.density_g_cm3,
                args.slab_thickness_cm, args.histories, args.seed
            )
            Path(args.out).parent.mkdir(parents=True, exist_ok=True)
            _write_json(args.out, result)
            print(json.dumps(result, indent=2))
            return 0 if result["status"] == "PASS" else 2
        if args.command == "validate-kn":
            result = validate_kn(args.energy_MeV, args.samples, args.seed)
            Path(args.out).parent.mkdir(parents=True, exist_ok=True)
            _write_json(args.out, result)
            print(json.dumps(result, indent=2))
            return 0 if result["status"] == "PASS" else 2
        if args.command == "run":
            manifest = load_and_verify_manifest(args.table, args.manifest)
            table = XCOMTable.from_csv(args.table)
            config = TransportConfig(
                photon_energy_MeV=args.energy_MeV,
                density_g_cm3=args.density_g_cm3,
                histories=args.histories,
                seed=args.seed,
            )
            summary = PhotonTransport(table, config).run(args.outdir)
            summary["xcom_provenance"] = {
                "source": manifest["source"],
                "material": manifest["material"],
                "sha256": manifest["sha256"],
                "retrieved_utc": manifest["retrieved_utc"],
            }
            _write_json(Path(args.outdir) / "run_summary.json", summary)
            print(json.dumps(summary, indent=2))
            return 0
    except (ValueError, OSError, XCOMDataError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
