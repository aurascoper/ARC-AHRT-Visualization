# ==============================================================================
# pymol/arc_ahrt_capsid_scene.py
# Title: PyMOL 3D Structural & Radiation-Target Visualization for ARC & AHRT
# Epistemic Status: [MEASURED GEOMETRY FROM DEPOSITED COORDINATES - NO MECHANISM CLAIMED]
#
# Usage (Headless 1920x1080 ray-traced export):
#   /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq pymol/arc_ahrt_capsid_scene.py
# Usage (Interactive GUI):
#   /Applications/PyMOL.app/Contents/MacOS/PyMOL pymol/arc_ahrt_capsid_scene.pml
# ==============================================================================

import hashlib
import json
from pathlib import Path
from pymol import cmd

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
OUT_DIR = ROOT / "pymol" / "renders"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def render_domain_and_asu_grid() -> Path:
    cmd.reinitialize()
    cmd.bg_color("black")
    cmd.set("ray_opaque_background", 1)
    cmd.set("orthoscopic", 1)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("grid_mode", 1)

    cmd.load(str(STRUCT / "6GSE.cif"), "rat_6GSE_monomer")
    cmd.load(str(STRUCT / "4X3H.cif"), "rat_4X3H_nlobe")
    cmd.load(str(STRUCT / "6TAP.cif"), "darc1_6TAP_asu")
    cmd.load(str(STRUCT / "6R24.cif"), "ty3_6R24_asu")

    cmd.hide("everything", "all")
    cmd.show("cartoon", "all")

    cmd.color("violetpurple", "rat_6GSE_monomer")
    cmd.color("orange", "rat_4X3H_nlobe")
    cmd.color("marine", "darc1_6TAP_asu and chain A+B+C+D")
    cmd.color("hotpink", "darc1_6TAP_asu and chain E")
    cmd.show("spheres", "darc1_6TAP_asu and (chain F or resn ZN)")
    cmd.color("tv_yellow", "darc1_6TAP_asu and (chain F or resn ZN)")
    cmd.set("sphere_scale", 2.2, "darc1_6TAP_asu and (chain F or resn ZN)")
    cmd.color("forest", "ty3_6R24_asu")

    cmd.orient("all")
    cmd.zoom("all", buffer=4.0)

    out_png = OUT_DIR / "pymol_arc_domain_and_asu_comparison.png"
    cmd.png(str(out_png), width=1920, height=1080, dpi=300, ray=1)
    return out_png


def render_6tap_icosahedral_radiation_target() -> Path:
    cmd.reinitialize()
    cmd.bg_color("black")
    cmd.set("ray_opaque_background", 1)
    cmd.set("orthoscopic", 1)
    cmd.set("antialias", 2)
    cmd.set("sphere_mode", 0)

    cmd.load(str(STRUCT / "6TAP.cif"), "asu_base")
    with open(STRUCT / "6TAP_operators.json") as f:
        ops = json.load(f)

    # Extract only C-alpha atoms for chains A-D, all backbone/C-alpha for chain E, and ZN ions
    # so the 60-copy full icosahedral capsid renders fast and cleanly at 1920x1080.
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

    # Front hemisphere CA shell made highly transparent so interior NC/Zn target is visible
    cmd.show("spheres", "capsid_ca_front")
    cmd.set("sphere_scale", 1.2, "capsid_ca_front")
    cmd.color("skyblue", "capsid_ca_front")
    cmd.set("sphere_transparency", 0.82, "capsid_ca_front")

    # Interior 60 NC zinc-knuckle chains (Chain E) and 60 Zn ions (Chain F)
    cmd.show("spheres", "capsid_nc_all")
    cmd.set("sphere_scale", 2.3, "capsid_nc_all")
    cmd.color("hotpink", "capsid_nc_all")

    cmd.show("spheres", "capsid_zn_all")
    cmd.set("sphere_scale", 3.6, "capsid_zn_all")
    cmd.color("tv_yellow", "capsid_zn_all")

    cmd.orient("all")
    cmd.zoom("all", buffer=8.0)

    out_png = OUT_DIR / "pymol_6tap_icosahedral_radiation_target.png"
    cmd.png(str(out_png), width=1920, height=1080, dpi=300, ray=1)
    return out_png


png1 = render_domain_and_asu_grid()
png2 = render_6tap_icosahedral_radiation_target()

receipt = {
    "schema": "arc-ahrt-pymol-render-receipt-v1",
    "operator_count": 60,
    "ca_subunits": 240,
    "nc_subunits": 60,
    "zn_ions": 60,
    "resolution_px": [1920, 1080],
    "renders": {
        png1.name: {
            "path": f"pymol/renders/{png1.name}",
            "bytes": png1.stat().st_size,
            "sha256": sha256_file(png1),
        },
        png2.name: {
            "path": f"pymol/renders/{png2.name}",
            "bytes": png2.stat().st_size,
            "sha256": sha256_file(png2),
        },
    },
}

with open(OUT_DIR / "pymol_render_receipt.json", "w") as f:
    json.dump(receipt, f, indent=2)

print("PyMOL renders complete:", json.dumps(receipt, indent=2))
