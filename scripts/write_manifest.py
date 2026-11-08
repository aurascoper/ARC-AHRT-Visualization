#!/usr/bin/env python3
"""Create a provenance manifest for a normalized XCOM CSV without editing it."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_path", type=Path)
    ap.add_argument("--raw-source", default="")
    ap.add_argument("--source", default="https://physics.nist.gov/PhysRefData/Xcom/")
    ap.add_argument("--material", default="H2O")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    if not args.csv_path.is_file():
        ap.error(f"CSV not found: {args.csv_path}")
    manifest = {
        "source": args.source,
        "material": args.material,
        "units": "cm2/g; energy_MeV",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "sha256": digest(args.csv_path),
        "raw_source_file": args.raw_source,
        "notes": "Normalized NIST XCOM output; verify coefficient columns and retain raw source output.",
    }
    target = args.out or args.csv_path.with_name("xcom_water_manifest.json")
    target.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(target)
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
