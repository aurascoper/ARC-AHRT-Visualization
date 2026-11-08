"""AHRT-OpenMC-Crosscheck-01 — independent photon-physics comparison vs Stage A.

Geometry/source matched exactly to ~/Developer/AHRT/config/cs137_reference.toml.
Source phase space is written with openmc.write_source_file from explicit
numpy sampling — x = -sqrt(R^2 - y^2), y~U(-R,R), z~U(-L/2,L/2), dir=+x —
so the beam map matches Stage A's sampler exactly (no stock-distribution
approximation).

Modes: --mode full (macro beam, unconditional estimator) | aimed (cell-footprint beam, conditional estimator)
Photons only. OpenMC electron treatment: led (local energy deposition) — declared scope.
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import openmc

E_SRC_MEV = 0.661657
RHO = 1.0
MACRO_R, MACRO_L = 0.3, 0.6
CELL_R, CELL_L = 5e-4, 2e-3
BASE = "/Users/aurascoper/Developer/AHRT"
OUT = os.path.join(BASE, "results", "openmc_crosscheck")
XS = os.path.join(BASE, "data", "openmc_xs", "official", "cross_sections_min.xml")


def write_beam_source(path, n, radius, half_length, seed):
    rng = np.random.default_rng(seed)
    y = rng.uniform(-radius, radius, n)
    z = rng.uniform(-half_length, half_length, n)
    x = -np.sqrt(np.maximum(0.0, radius**2 - y**2))
    E = np.full(n, E_SRC_MEV * 1e6)
    parts = [
        openmc.SourceParticle(
            r=(float(x[i]), float(y[i]), float(z[i])),
            u=(1.0, 0.0, 0.0),
            E=float(E[i]),
            wgt=1.0,
            particle=openmc.ParticleType.PHOTON,
        )
        for i in range(n)
    ]
    openmc.write_source_file(parts, path)


def build_model(mode: str, n_histories: int, seed: int):
    os.environ["OPENMC_CROSS_SECTIONS"] = XS
    h2o = openmc.Material(name="water")
    h2o.add_nuclide("H1", 2.0); h2o.add_nuclide("O16", 1.0)
    h2o.set_density("g/cm3", RHO)
    mats = openmc.Materials([h2o])

    cyl = openmc.ZCylinder(r=MACRO_R, boundary_type="vacuum")
    ztop = openmc.ZPlane(z0=MACRO_L / 2, boundary_type="vacuum")
    zbot = openmc.ZPlane(z0=-MACRO_L / 2, boundary_type="vacuum")

    cell_cyl = openmc.ZCylinder(r=CELL_R)
    cell_top = openmc.ZPlane(z0=CELL_L / 2)
    cell_bot = openmc.ZPlane(z0=-CELL_L / 2)
    cell_roi = openmc.Cell(fill=h2o, region=-cell_cyl & -cell_top & +cell_bot, name="cell_roi")
    macro_cell = openmc.Cell(fill=h2o, region=(-cyl & -ztop & +zbot) & ~cell_roi.region, name="macro")
    geom = openmc.Geometry([cell_roi, macro_cell])

    ns = min(n_histories, 200_000)
    srcfile = os.path.join(OUT, f"source_phase_{mode}.h5")
    if mode == "full":
        write_beam_source(srcfile, ns, MACRO_R, MACRO_L / 2, seed)
    else:
        write_beam_source(srcfile, ns, CELL_R if mode == "aimed" else MACRO_R, CELL_L / 2 if mode == "aimed" else MACRO_L / 2, seed)
    src = openmc.FileSource(srcfile); src.particle = "photon"

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.batches = 10
    settings.particles = n_histories  # total across batches
    settings.seed = seed
    settings.photon_transport = True
    settings.electron_treatment = "led"
    settings.source = src
    if mode == "tracks":
        settings.track = [(1, 1, i) for i in range(1, 26)]  # first 25 primaries

    scores = ["scatter", "coherent-scatter", "incoherent-scatter", "photoelectric", "pair-production"]
    t_cell = openmc.Tally(name="cell"); t_cell.filters = [openmc.CellFilter(cell_roi)]; t_cell.scores = scores + ["heating"]
    t_macro = openmc.Tally(name="macro"); t_macro.filters = [openmc.CellFilter(macro_cell)]; t_macro.scores = scores + ["heating"]
    return openmc.Model(geometry=geom, materials=mats, settings=settings, tallies=openmc.Tallies([t_cell, t_macro]))


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["full", "aimed", "tracks"], default="full")
    ap.add_argument("--histories", type=int, default=100_000)
    ap.add_argument("--seed", type=int, default=20261104)  # independent of Stage A's 20261103
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    model = build_model(args.mode, args.histories, args.seed)
    spath = model.run(cwd=OUT)
    print("statepoint:", spath)

    sp = openmc.StatePoint(spath)
    result = {"mode": args.mode, "histories": args.histories, "seed": args.seed, "tallies": {}}
    for tname in ("cell", "macro"):
        t = sp.get_tally(name=tname)
        df = t.get_pandas_dataframe()
        cols = list(df.columns)
        stdcol = "std. dev." if "std. dev." in cols else None
        result["tallies"][tname] = {
            r["score"]: {"mean": float(r["mean"]), "std": float(r[stdcol]) if stdcol else None}
            for _, r in df.iterrows()}
    outname = os.path.join(OUT, f"tallies_{args.mode}.json")
    json.dump(result, open(outname, "w"), indent=2)
    print("wrote", outname)


if __name__ == "__main__":
    main()
