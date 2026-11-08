#!/usr/bin/env python3
"""M7ArcClusterBridge Level 1 v2 - corrected per PR #6 review.
Fixes: real inventory validation, residue-level interface contacts, sign-invariant orientation,
k-sweep with full diagnostics, symmetry-aware agreement, standard silhouette, vendored/CLI inputs.
evidence_class=CALCULATED_EXPLORATORY.
"""
import argparse, hashlib, json
import numpy as np
from scipy.spatial import cKDTree

p = argparse.ArgumentParser()
p.add_argument("--cif", default="structures/6TAP.cif")
p.add_argument("--ops", default="structures/6TAP_operators.json")
p.add_argument("--out", default="artifacts")
p.add_argument("--seed", type=int, default=20261105)
p.add_argument("--seeds", type=int, default=25)
p.add_argument("--ks", default="2,3,4,6")
args = p.parse_args()

sha = lambda p: hashlib.sha256(open(p,'rb').read()).hexdigest()
print("source hashes:", sha(args.cif)[:16], sha(args.ops)[:16])

# ---- mmCIF parse (minimal but real: header-driven; validate chain inventory) ----
hdr = None; cols = None
ca = None; rows = []
ent = {}
with open(args.cif) as f:
    for ln in f:
        if ln.startswith("_atom_site."):
            pass
        if ln.startswith(("ATOM","HETATM")):
            f_ = ln.split()
            rows.append(f_)
# use header order
hdr_lines = [l.split()[0] for l in open(args.cif) if l.strip().startswith("_atom_site.")]
ix = hdr_lines.index("_atom_site.Cartn_x"); iy = hdr_lines.index("_atom_site.Cartn_y"); iz = hdr_lines.index("_atom_site.Cartn_z")
ich = hdr_lines.index("_atom_site.label_asym_id"); ires_seq = hdr_lines.index("_atom_site.label_seq_id")
iatom = hdr_lines.index("_atom_site.label_atom_id"); ires = hdr_lines.index("_atom_site.label_comp_id")
ient = hdr_lines.index("_atom_site.label_entity_id");

import numpy as np
xyz = np.array([[float(r[ix]), float(r[iy]), float(r[iz])] for r in rows])
chains = np.array([r[ich] for r in rows])
resid = np.array([(r[ires], int(r[ires_seq]) if r[ires_seq] not in (".","?") else -1) for r in rows], dtype=object)
atypes = np.array([r[iatom] for r in rows])
ents = np.array([r[ient] for r in rows])

# inventory from deposit: per-CHAIN classification is authoritative (entities vary)
n_protein_atoms = int((np.isin(chains, ["A","B","C","D","E"])).sum())
zn_atoms = int((chains == "F").sum())
residue_uniq_per_chain = {}
for c in np.unique(chains[chains != "F"]):
    rs = resid[(chains == c)]
    residue_uniq_per_chain[c] = len({(str(x[0]), int(x[1])) for x in rs})
# protein residues per chain (chain E = NC fragment)
chain_class = {c: ("CA" if c in "ABCD" else ("NC" if c == "E" else "Zn") if c == "F" else "OTHER") for c in set(chains.flat if hasattr(chains,'flat') else chains.tolist())}

ops = json.load(open(args.ops))
n_ops = len(ops)
print(f"deposit: {n_protein_atoms} protein/NC atoms, {zn_atoms} Zn, chains={sorted(set(chains[ents!='3']))}, ops={n_ops}")

# assembly inventory (explicit, independently checked)
mat = np.array([o["matrix"] for o in ops]).reshape(n_ops,3,3)
vec = np.array([o["vector"] for o in ops]).reshape(n_ops,3)
CA_chains = [c for c in chain_class if chain_class[c]=="CA"]
inst = []            # CA instances (chain,op)
nc_inst = []         # NC fragment instances
zn_inst = []         # Zn positions
mask_protein = ents != "3"
sel_nc = (chains == "E")
sel_zn = (chains == "F")
mask_protein = (chains != "F") & (chains != "E")
for oi in range(n_ops):
    for c in CA_chains:
        m = (chains == c) & mask_protein
        inst.append((oi, c, (mat[oi] @ xyz[m].T).T + vec[oi]))
    nc_inst.append((mat[oi] @ xyz[sel_nc].T).T + vec[oi])
    zn_inst.append((mat[oi] @ xyz[sel_zn].T).T + vec[oi])
inv = {"CA_instances": len(inst), "NC_instances": len(nc_inst), "Zn_instances": len(zn_inst),
       "nc_atoms_per_instance": int(len(xyz[sel_nc])), "zn_atoms_per_instance": int(len(xyz[sel_zn])),
       "expected": {"CA": 240, "NC": 60, "Zn": 60}}
assert inv["CA_instances"] == 240 and inv["NC_instances"] == 60 and inv["Zn_instances"] == 60, inv
assert inv["nc_atoms_per_instance"] > 0 and inv["zn_atoms_per_instance"] == 1, "vacuous inventory check forbidden"
print("inventory:", inv["CA_instances"], inv["NC_instances"], inv["Zn_instances"], "OK (240/60/60)")

# ---- per-instance features (sign-invariant orientation + residue-level contacts) ----
def inst_feats(inst):
    feats = []
    for oi, c, P in inst:
        cen = P.mean(axis=0)
        X = P - cen
        v0 = np.linalg.svd(X, full_matrices=False)[2][0]
        feats.append((oi, c, cen, np.abs(v0)))
    return feats

caps = None
feats = inst_feats(inst)
cen = np.array([f[2] for f in feats])
cap_ctr = cen.mean(axis=0)
radial = np.linalg.norm(cen - cap_ctr, axis=1)
# orientation feature: |cos(axis_PC1, radial_dir)|  (sign-invariant)
rd = (cen - cap_ctr); rd = rd / np.linalg.norm(rd, axis=1, keepdims=True)
ori_abs = []
for i, f in enumerate(feats):
    v0 = np.linalg.svd(inst[i][2] - f[2], full_matrices=False)[2][0]
    ori_abs.append([abs(v0 @ rd[i]), abs(v0 @ np.array([1,0,0])), abs(v0 @ np.array([0,0,1]))])
ori_abs = np.array(ori_abs)

# residue-level interface contacts: all-atom 5 A cross-instance atom pairs -> unique residue-residue pairs
allP = np.vstack([t[2] for t in inst]); iid = np.concatenate([np.full(len(t[2]), i) for i,t in enumerate(inst)])
tree = cKDTree(allP)
pairs = tree.query_pairs(5.0, output_type='ndarray')
pi, pj = iid[pairs[:,0]], iid[pairs[:,1]]
cross = pi != pj
pair_set = {(min(a,b), max(a,b)) for a,b in zip(pi[cross], pj[cross])}  # instance-level interface pairs
iface_deg = np.zeros(240, dtype=int)
for a,b in pair_set:
    iface_deg[a] += 1; iface_deg[b] += 1

# residue-contact count: pairs of residues with distance<5 across different instances
# atom -> (instance, chain-local res id)
resid_flat = [(str(x[0]), int(x[1])) for x in resid]
res_keys = []
ci = 0
for oi, c in [(t[0], t[1]) for t in inst]:
    m_idx = np.where(chains == c)[0]
    for mj in m_idx:
        res_keys.append((ci, resid_flat[mj][0], resid_flat[mj][1]))
    ci += 1
res_keys = [tuple(x) for x in res_keys]
# count unique cross-instance RESIDUE pairs (instance, resname, resseq on both ends)
rr = {(res_keys[u], res_keys[v]) for u,v in zip(pairs[cross][:,0], pairs[cross][:,1])
      if res_keys[u] != res_keys[v]}
res_pair_cnt = np.zeros(240, dtype=int)
for a,b in rr:
    res_pair_cnt[a[0]] += 1; res_pair_cnt[b[0]] += 1

X = np.column_stack([cen - cap_ctr, radial, ori_abs, iface_deg.astype(float), res_pair_cnt.astype(float)])
np.save(f"{args.out}/arc_features_level1_v2.npy", X)
feat_def = "xyz_rel_capsid_centroid(3) + radial_distance(1) + |cos|(PC1,radial)(1)+|cos|(PC1,X)(1)+|cos|(PC1,Y)(1) + interface_instances(1) + cross-residue-contact-count(1)"
mu, sd_ = X.mean(axis=0), X.std(axis=0); sd_[sd_==0]=1
XZ = (X-mu)/sd_

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

results_k = {}
all_partitions = {}
ks = [int(k) for k in args.ks.split(",")]
for k in ks:
    partitions = []
    cents = []
    for s in range(args.seed, args.seed + args.seeds):
        km = KMeans(n_clusters=k, n_init=1, random_state=s, max_iter=300).fit(XZ)
        partitions.append(km.labels_)
    partitions = np.array(partitions)
    # stability: fraction of same-pair co-assignments consistent across seeds
    coassign = []
    for i in range(args.seeds):
        for j in range(i+1, args.seeds):
            a, b = partitions[i], partitions[j]
            coassign.append(np.mean((a[:,None]==a[None,:]) == (b[:,None]==b[None,:])))
    stab = float(np.mean(coassign))
    sil = float(silhouette_score(XZ, partitions[0]))
    # symmetry agreement: fraction of (same chain label, different operator) instance pairs co-assigned
    labels_per_chain = {}
    for idx, (oi, c) in enumerate([(t[0], t[1]) for t in inst]):
        labels_per_chain.setdefault(c, []).append((oi, partitions[0][idx]))
    within_agree = {}
    for c, lst in labels_per_chain.items():
        lab = np.array([l for _, l in lst])
        both = []
        for a in range(60):
            for b in range(a+1, 60):
                both.append(1.0 if lab[a]==lab[b] else 0.0)
        within_agree[c] = float(np.mean(both))
    results_k[k] = {"silhouette": sil, "pairwise_seed_stability": stab,
                    "cluster_sizes": [int((partitions[0]==c).sum()) for c in range(k)],
                    "symmetry_within_chain_agreement": within_agree}
    all_partitions[k] = partitions

receipt = {
 "module":"M7ArcClusterBridge Level 1 v2", "evidence_class":"CALCULATED_EXPLORATORY",
 "corrections_from_review":"residue-level interfaces; sign-invariant orientation; inventory independently asserted (240 CA/60 NC/60 Zn); sklearn KMeans+silhouette; k-sweep; symmetry agreement; vendored inputs",
 "inputs":{"cif_sha256": sha(args.cif), "ops_sha256": sha(args.ops), "cif_path": args.cif, "ops_path": args.ops},
 "feature_definition": feat_def, "inventory": inv,
 "k_sweep": results_k,
 "feature_matrix_sha256": hashlib.sha256(open(f"{args.out}/arc_features_level1_v2.npy","rb").read()).hexdigest(),
 "caveats": ["geometric clusters != structural modules != biological functions; symmetry agreement measures self-consistency under chain labels, NOT icosahedral correctness"],
}
json.dump({k: v.tolist() for k,v in all_partitions.items()}, open(f"{args.out}/arc_level1_v2_partitions.json","w"))
json.dump(receipt, open(f"{args.out}/arc_level1_v2_receipt.json","w"), indent=1)
for k in ks:
    print(f"k={k}: sil={results_k[k]['silhouette']:.3f} stab={results_k[k]['pairwise_seed_stability']:.3f} sizes={results_k[k]['cluster_sizes']} chain-agree={ {c:round(v,2) for c,v in results_k[k]['symmetry_within_chain_agreement'].items()} }")
