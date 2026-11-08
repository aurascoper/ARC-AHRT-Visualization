# Portable ChimeraX Python script to generate the three AHRT scenes:
# 1. arc_structural_reference.cxs + render
# 2. arc_heterarchy.cxs + render
# 3. arc_radiation_ready.cxs + render
#
# Run with:
#   /Applications/ChimeraX-1.12.app/Contents/MacOS/ChimeraX --nogui --offscreen workbench/scene_ahrt_01_portable.py

import json
from pathlib import Path
from chimerax.core.commands import run

def resolve_repo_root() -> Path:
    candidates = []
    if "__file__" in globals():
        p = Path(__file__).resolve()
        candidates.extend([p.parent.parent, p.parent])
    cwd = Path.cwd().resolve()
    candidates.extend([cwd, cwd.parent])
    for cand in candidates:
        if (cand / "structures" / "6TAP.cif").exists():
            return cand
    raise RuntimeError("Could not locate structures/6TAP.cif from script path or working directory.")


ROOT = resolve_repo_root()
STRUCT = ROOT / "structures"
WB = ROOT / "workbench"
RENDERS = WB / "renders"
RENDERS.mkdir(parents=True, exist_ok=True)


def setup_scene():
    run(session, "close session")
    run(session, "set bgColor #080d14")
    run(session, "lighting soft")
    run(session, "graphics silhouettes true")


# ==============================================================================
# SCENE 1: arc_structural_reference
# ==============================================================================
setup_scene()
run(session, f'open "{STRUCT / "6GSE.cif"}"')
run(session, f'open "{STRUCT / "4X3H.cif"}"')
run(session, "hide atoms")
run(session, "show cartoons")
run(session, "color #1 #38bdf8")
run(session, "color #2 #f97316")
run(session, "tile #1-2 columns 2 spacingFactor 1.6")
run(session, "view orient")
run(
    session,
    '2dlabels create title text "Scene 1: Mammalian Arc Structural Reference (Domain Only)" '
    "xpos 0.02 ypos 0.95 size 16 color #f8fafc",
)
run(
    session,
    '2dlabels create sub1 text "Left: 6GSE (Rat Arc CA domain monomer, solution NMR, C1) | '
    'Right: 4X3H (Rat Arc N-lobe + Stargazin peptide, X-ray 2.40 A)" '
    "xpos 0.02 ypos 0.90 size 13 color #94a3b8",
)
run(
    session,
    '2dlabels create caveat text "EPISTEMIC NOTE: Neither entry is an icosahedral capsid lattice. '
    'No high-resolution mammalian Arc full capsid structure is deposited in the PDB." '
    "xpos 0.02 ypos 0.03 size 12 color #fbbf24",
)
run(session, f'save "{WB / "arc_structural_reference.cxs"}"')
run(
    session,
    f'save "{RENDERS / "arc_structural_reference.png"}" width 1920 height 1080 supersample 3',
)

# ==============================================================================
# SCENE 2: arc_heterarchy
# ==============================================================================
setup_scene()
run(session, f'open "{STRUCT / "6GSE.cif"}"')
run(session, f'open "{STRUCT / "6TAP.cif"}"')
run(session, f'open "{STRUCT / "6R24.cif"}"')
run(session, f'open "{STRUCT / "1OHF.cif"}"')
run(session, "hide atoms")
run(session, "show cartoons")
run(session, "color #1 #a855f7")
run(session, "color #2 #38bdf8")
run(session, "color #3 #f97316")
run(session, "color #4 #22c55e")
run(session, "tile #1-4 columns 4 spacingFactor 1.5")
run(session, "view orient")
run(
    session,
    '2dlabels create title text "Scene 2: Arc Heterarchy and Capsid Comparison (ASU View)" '
    "xpos 0.02 ypos 0.95 size 16 color #f8fafc",
)
run(
    session,
    '2dlabels create sub1 text "6GSE (Rat monomer) | 6TAP (dArc1 ASU: 4 CA + 1 NC + 1 Zn) | '
    '6R24 (Ty3/gypsy retrotransposon ASU) | 1OHF (NwV T=4 reference ASU)" '
    "xpos 0.02 ypos 0.90 size 12 color #94a3b8",
)
run(
    session,
    '2dlabels create caveat text "EPISTEMIC NOTE: Geometric comparison only. '
    '6TAP has 240 CA shell subunits (T=4) plus 60 interior NC zinc-knuckle fragments (chain E) and 60 Zn ions." '
    "xpos 0.02 ypos 0.03 size 12 color #fbbf24",
)
run(session, f'save "{WB / "arc_heterarchy.cxs"}"')
run(
    session,
    f'save "{RENDERS / "arc_heterarchy.png"}" width 1920 height 1080 supersample 3',
)

# ==============================================================================
# SCENE 3: arc_radiation_ready (6TAP Full Capsid + Interior NC/Zn Target)
# ==============================================================================
setup_scene()
run(session, f'open "{STRUCT / "6TAP.cif"}"')

with open(STRUCT / "6TAP_operators.json") as f:
    ops = json.load(f)

run(session, "hide #1 atoms")
run(session, "show #1 cartoons")
run(session, "color #1/A-D #38bdf8")
run(session, "color #1/E #ec4899")
run(session, "show #1/F atoms")
run(session, "style #1/F sphere")
run(session, "color #1/F #facc15")
run(session, "size #1/F atomRadius 3.5")

for i, op in enumerate(ops[1:], start=2):
    m = op["matrix"]
    v = op["vector"]
    mat_str = (
        f"{m[0]:.6f},{m[1]:.6f},{m[2]:.6f},{v[0]:.4f},"
        f"{m[3]:.6f},{m[4]:.6f},{m[5]:.6f},{v[1]:.4f},"
        f"{m[6]:.6f},{m[7]:.6f},{m[8]:.6f},{v[2]:.4f}"
    )
    run(session, "combine #1 modelId " + str(i) + ' name "op_' + str(i) + '"')
    run(session, f"view matrix models #{i},{mat_str}")
    if m[2] * 0 + m[5] * 0 + m[8] * 0 + v[2] > 350.0:
        run(session, f"transparency #{i}/A-D 75 cartoons")
    else:
        run(session, f"transparency #{i}/A-D 20 cartoons")

run(session, "view orient")
run(session, "zoom 1.05")
run(
    session,
    '2dlabels create title text "Scene 3: AHRT v0.1 Radiation Target (6TAP dArc1 Capsid + Interior NC/Zn)" '
    "xpos 0.02 ypos 0.95 size 16 color #f8fafc",
)
run(
    session,
    '2dlabels create sub1 text "Outer Shell: 240 CA subunits (Chains A-D, blue, R_in=15.5 nm, R_out=19.9 nm) | '
    'Interior: 60 NC Zn-knuckle fragments (Chain E, pink) + 60 Zn ions (Chain F, yellow spheres)" '
    "xpos 0.02 ypos 0.90 size 12 color #94a3b8",
)
run(
    session,
    '2dlabels create caveat text "AHRT v0.1 FROZEN BASELINE: Stage A (XCOM/KN photon) -> Stage B (TOPAS-nBio track/IRT) -> '
    'Stage C (damage manifold). RNA core is NOT resolved in 6TAP (parameterized separately)." '
    "xpos 0.02 ypos 0.03 size 11 color #fbbf24",
)
run(session, f'save "{WB / "arc_radiation_ready.cxs"}"')
run(
    session,
    f'save "{RENDERS / "arc_radiation_ready.png"}" width 1920 height 1080 supersample 3',
)
