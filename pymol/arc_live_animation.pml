# ==============================================================================
# pymol/arc_live_animation.pml
# Title: Live Interactive PyMOL 360-Degree Turntable for 6TAP Radiation Target
# Epistemic Status: [MEASURED GEOMETRY FROM DEPOSITED COORDINATES - NO MECHANISM CLAIMED]
# ==============================================================================

python
import json
from pathlib import Path
from pymol import cmd, util

cwd = Path.cwd().resolve()
root = cwd if (cwd / "structures" / "6TAP.cif").exists() else cwd.parent
struct_dir = root / "structures"

cmd.reinitialize()
cmd.bg_color("black")
cmd.set("orthoscopic", 1)
cmd.set("antialias", 2)

cmd.load(str(struct_dir / "6TAP.cif"), "asu_base")
with open(struct_dir / "6TAP_operators.json") as f:
    ops = json.load(f)

cmd.select("base_ca_shell", "asu_base and chain A+B+C+D and name CA")
cmd.select("base_nc_core", "asu_base and chain E and name CA")
cmd.select("base_zn_ions", "asu_base and (chain F or resn ZN)")

cmd.create("op_1_ca", "base_ca_shell")
cmd.create("op_1_nc", "base_nc_core")
cmd.create("op_1_zn", "base_zn_ions")
cmd.delete("asu_base")

front_ops = []
back_ops = []

for idx, op in enumerate(ops, start=1):
    m = op["matrix"]
    v = op["vector"]
    if idx > 1:
        mat16 = [
            m[0], m[1], m[2], v[0],
            m[3], m[4], m[5], v[1],
            m[6], m[7], m[8], v[2],
            0.0,  0.0,  0.0,  1.0,
        ]
        ca_name = f"op_{idx}_ca"
        nc_name = f"op_{idx}_nc"
        zn_name = f"op_{idx}_zn"
        cmd.copy(ca_name, "op_1_ca")
        cmd.copy(nc_name, "op_1_nc")
        cmd.copy(zn_name, "op_1_zn")
        cmd.transform_object(ca_name, mat16)
        cmd.transform_object(nc_name, mat16)
        cmd.transform_object(zn_name, mat16)
    else:
        ca_name = "op_1_ca"

    if v[2] > 350.0:
        front_ops.append(ca_name)
    else:
        back_ops.append(ca_name)

cmd.group("capsid_ca_back", " ".join(back_ops))
cmd.group("capsid_ca_front", " ".join(front_ops))
cmd.group("capsid_nc_all", "op_*_nc")
cmd.group("capsid_zn_all", "op_*_zn")

cmd.hide("everything", "all")
cmd.show("spheres", "capsid_ca_back")
cmd.set("sphere_scale", 1.8, "capsid_ca_back")
cmd.color("deepteal", "capsid_ca_back")

cmd.show("spheres", "capsid_ca_front")
cmd.set("sphere_scale", 1.2, "capsid_ca_front")
cmd.color("skyblue", "capsid_ca_front")
cmd.set("sphere_transparency", 0.82, "capsid_ca_front")

cmd.show("spheres", "capsid_nc_all")
cmd.set("sphere_scale", 2.3, "capsid_nc_all")
cmd.color("hotpink", "capsid_nc_all")

cmd.show("spheres", "capsid_zn_all")
cmd.set("sphere_scale", 3.6, "capsid_zn_all")
cmd.color("tv_yellow", "capsid_zn_all")

cmd.orient("all")
cmd.zoom("all", buffer=8.0)

# Continuous 360-degree turntable movie loop in PyMOL GUI
cmd.mset("1 x360")
util.mroll(1, 360, 1)
cmd.set("movie_loop", 1)
cmd.mplay()
python end
