# AHRT-ChimeraX-Scene-01 -- chain-class-aware 6TAP dArc1 assembly
# Runs headless: ChimeraX --nogui --offscreen --exit --script scene_ahrt_01.py
# Produces three scene states + verification JSON printed to stdout.
# Epistemic status: deposited coordinates + deposit's own assembly operators only.
# NO fabricated trajectories, NO rescaled distances, NO mammalian-Arc claims.

import json, sys
from chimerax.core.commands import run

STRUCT = "/Users/aurascoper/Developer/ArcAdvisor/structures/6TAP.cif"
OUT = "/Users/aurascoper/Developer/ArcAdvisor/workbench"
LOG = {}

run(session, "close")
run(session, f'open "{STRUCT}"')
run(session, "sym #1 assembly 1 copies true")

# ---- inventory: count chain-class instances across COPIES only (exclude deposited #1) ----
from chimerax.atomic import AtomicStructure
all_models = session.models.list(type=AtomicStructure)
# deposited model is #1; sym copies are #2.x children
copies = [m for m in all_models if m.id[0] != 1]
run(session, "hide #1 models")  # original counted separately; assembly copies carry the scene
ca_models = 0
nc_models = 0
zn_atoms = 0
asym_chain_ids = set()
for m in copies:
    for c in m.chains:
        asym_chain_ids.add(c.chain_id)
        if c.chain_id in ("A", "B", "C", "D"):
            ca_models += 1
        elif c.chain_id == "E":
            nc_models += 1
    # zinc: element Zn, not chain label
    try:
        zn_atoms += int((m.atoms.elements.number == 30).sum())
    except Exception:
        pass
LOG["num_assembly_copy_models"] = len(copies)
LOG["asym_chain_ids_seen"] = sorted(asym_chain_ids)
LOG["ca_chain_instances"] = ca_models
LOG["nc_fragment_instances"] = nc_models
LOG["zn_sites"] = zn_atoms
LOG["polymer_chain_instances_total"] = ca_models + nc_models
LOG["expected"] = {"ca": 240, "nc": 60, "zn": 60}
LOG["counts_match_expected"] = (ca_models == 240 and nc_models == 60)

# ---- zinc sites: deposit has ONE Zn (chain F, coords E:Cys227 SG at 2.333 A, cif metalc1).
# ChimeraX sym does not replicate HETATM into copies -> place 60 Zn markers manually
# by applying the deposit's own 60 assembly operators (pdbx_struct_oper_list) to the
# deposited Zn coordinate. Method declared, not hidden.
import re as _re, json as _json
_cif = open(STRUCT).read()
_zn_m = _re.search(r'^HETATM\s+\d+\s+ZN\s+ZN.*?(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', _cif, _re.M)
_zn = (float(_zn_m.group(1)), float(_zn_m.group(2)), float(_zn_m.group(3)))
_matrices = [(o["matrix"], o["vector"]) for o in _json.load(open(STRUCT.replace("6TAP.cif","6TAP_operators.json")))]
LOG["n_assembly_operators_parsed"] = len(_matrices)

from chimerax.markers import MarkerSet
_zset = MarkerSet(session, "zinc_sites_manual_from_operators")
_yellow = (255, 255, 80, 255)
for (_m, _t) in _matrices:
    _p = (
        _m[0]*_zn[0] + _m[1]*_zn[1] + _m[2]*_zn[2] + _t[0],
        _m[3]*_zn[0] + _m[4]*_zn[1] + _m[5]*_zn[2] + _t[1],
        _m[6]*_zn[0] + _m[7]*_zn[1] + _m[8]*_zn[2] + _t[2],
    )
    _zset.create_marker(_p, _yellow, 2.2)
session.models.add([_zset])
LOG["zn_deposit_sites"] = 1
LOG["zn_coordination_reference"] = "metalc1: chain E Cys227 SG, 2.333 A (cif)"
LOG["zn_assembly_markers_placed"] = len(_matrices)
LOG["zn_method_note"] = "single deposited Zn x60 deposit operators; ChimeraX does not replicate HETATM in sym copies"

# residue range sanity for one NC fragment (chain E: 224-250)
e_res = []
for m in copies:
    for c in m.chains:
        if c.chain_id == "E":
            e_res = sorted(set(r.number for r in c.residues if r is not None))
            break
    if e_res:
        break
LOG["chain_E_residue_range"] = [e_res[0], e_res[-1]] if e_res else None
LOG["chain_E_range_ok"] = bool(e_res and e_res[0] >= 220 and e_res[-1] <= 255)



def try_png(session, cmd, log_dict, key):
    from chimerax.core.errors import LimitationError
    try:
        run(session, cmd)
        log_dict[key] = "saved"
    except LimitationError as e:
        log_dict[key] = f"SKIPPED: {e}"

# ---- styling: chain-class colors ----
run(session, "hide atoms")
run(session, "cartoon")
run(session, "color /A-D cornflower blue")
run(session, "color /E magenta")
run(session, "color /F yellow")
run(session, "style /F sphere")
run(session, "size atomRadius 1.6")
run(session, "set bgColor white")
run(session, "lighting soft")
run(session, "windowsize 1920 1080")

# ---- named selections (chain classes + sites) ----
run(session, 'name ca_shell /A-D')
run(session, 'name nc_fragments /E')
run(session, 'name zinc_sites /F')

# ---- scene 1: structural reference ----
run(session, "clip off")
run(session, "view")
try_png(session, f'save "{OUT}/renders/arc_structural_reference.png" width 1920 height 1080 supersample 2', LOG, "png_arc_structural_reference")
run(session, f'save "{OUT}/arc_structural_reference.cxs" includeMaps false')

# ---- scene 2: heterarchy (contacts network, coordinate-derived only) ----
run(session, "contacts nc_fragments restrict cross interModel false overlapCutoff -0.4 hbondAllowance 0.0 reveal true color medium purple")
LOG["contacts_scene_command"] = "contacts nc_fragments restrict cross"
try_png(session, f'save "{OUT}/renders/arc_heterarchy.png" width 1920 height 1080 supersample 2', LOG, "png_arc_heterarchy")
run(session, f'save "{OUT}/arc_heterarchy.cxs" includeMaps false')

# ---- scene 3: radiation-ready (clean labels, empty ROI annotation layer stub) ----
run(session, "hide pseudobonds")
run(session, "hide cartoons")
run(session, "surface /A-D resolution 8")
run(session, "color /A-D #87ceeb70")  # translucent CA surface
run(session, "show /E ribbons")
run(session, "show /F atoms")
# empty marker placeholder layer for future tracks/exposure (no fabricated geometry)
from chimerax.markers import MarkerSet
roi = MarkerSet(session, "ahrt_roi_empty_annotations")
session.models.add([roi])
run(session, "view")
try_png(session, f'save "{OUT}/renders/arc_radiation_ready.png" width 1920 height 1080 supersample 2', LOG, "png_arc_radiation_ready")
run(session, f'save "{OUT}/arc_radiation_ready.cxs" includeMaps false')

LOG["scenes_saved"] = [
    "arc_structural_reference.cxs",
    "arc_heterarchy.cxs",
    "arc_radiation_ready.cxs",
]
LOG["renders"] = [
    "renders/arc_structural_reference.png",
    "renders/arc_heterarchy.png",
    "renders/arc_radiation_ready.png",
]
print("AHRT_SCENE01_VERIFY " + json.dumps(LOG, sort_keys=True))
sys.stdout.flush()
