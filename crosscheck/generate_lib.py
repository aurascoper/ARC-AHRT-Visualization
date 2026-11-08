"""Generate minimal H/O photon library from ENDF/B-VII.1 (provenance-pinned)."""
import os, hashlib, json, datetime
import openmc.data

S = os.path.expanduser("~/.hermes/cache/scratch")
OUT = "data/openmc_xs"
manifest = {"generated_utc": datetime.datetime.utcnow().isoformat()+"Z",
            "eval_source": "ENDF/B-VII.1 photoat + atomic_relax (NNDC endf-b7.1/zips, downloaded 2026-11-04)", "files": {}}
for nuc, pf, af in [("H1", "photoat/photoat-001_H_000.endf", "atomic_relax/atom-001_H_000.endf"),
                    ("O16", "photoat/photoat-008_O_000.endf", "atomic_relax/atom-008_O_000.endf")]:
    pp = os.path.join(S, "endf_photoat", pf); ap = os.path.join(S, "endf_atomrelax", af)
    photon = openmc.data.IncidentPhoton.from_endf(pp)
    relax = openmc.data.AtomicRelaxation.from_endf(ap)
    photon.atomic_relaxation = relax
    path = os.path.join(OUT, f"{nuc}.h5")
    photon.export_to_hdf5(path, "w")
    manifest["files"][nuc] = {"path": path, "sha256": hashlib.sha256(open(path,'rb').read()).hexdigest()[:16]}
    print(nuc, "->", path)
lib = openmc.data.DataLibrary()
for nuc in ["H1", "O16"]:
    lib.register_file(os.path.abspath(os.path.join(OUT, f"{nuc}.h5")))
lib.export_to_xml(os.path.join(OUT, "cross_sections.xml"))
json.dump(manifest, open(os.path.join(OUT, "xs_provenance.json"), "w"), indent=2)
print("library xml written")
