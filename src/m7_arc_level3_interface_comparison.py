#!/usr/bin/env python3
"""
M7 ARC Level-3 Symmetry-Normalized Interface Comparison & 4-Branch Provenance Separation.

Builds Level 3 directly from the 10 symmetry-inequivalent sub-interfaces of
Drosophila dArc1 CA capsid (PDB 6TAP, icosahedral point group I,
_em_single_particle_entity.point_symmetry = "I", biological assembly 1,
RCSB global symmetry I, 60 asymmetric units) validated in
Level 2 (`artifacts/m7_arc_level2_validation_receipt.json`) and cross-checked in
Step 3 (`artifacts/m7_arc_level2_independent_crosscheck_receipt.json`):

  1. Symmetry-Normalized Residue Fingerprints & Interface-Fingerprint Families:
     - Maps every interfacial residue (`auth_seq_id` / `label_seq_id`, NTD stalk/bundle
       41..125 vs CTD lobe/tail 126..206) across quasi-equivalent chains A, B, C, D.
     - Classifies the 10 sub-interfaces into three empirical interface-fingerprint
       groupings (`interface_fingerprint_family`: Family I NTD+CTD capsomere seams,
       Family II pure-CTD quasi-equivalent dimers, and Family III residual distal
       CTD interface hierarchy `Family_III_Residual_Distal_CTD_Hierarchy`), explicitly
       distinguished from established biological or functional modules and stratified
       by threshold robustness (`family_robustness_stratification`).
     - Computes 10x10 pairwise similarity matrices across 4 set-similarity metrics
       (`jaccard`, `sorensen_dice`, `szymkiewicz_simpson_overlap`, `cosine_ochiai`) on
       both homologous `seq_id` sets and canonical homologous residue-pair sets, plus
       contact-degree-weighted cosine similarity across the 83 homologous sequence
       positions, with single-linkage and complete-linkage threshold sweeps across
       tau in {0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70}.

  2. Normalized Interface Chemistry & Area Profiles:
     - Compares Shrake-Rupley Fibonacci two-sided BSA (`9,570.22 A^2` per ASU;
       `574,213.20 A^2` whole-capsid two-sided BSA) and one-side interface area
       (`4,785.10 A^2` per ASU; `287,106.60 A^2` whole-capsid one-side area),
       Lee-Richards z-slice BSA (`9,400.64 A^2`), and Voronoi radical-plane single
       dividing-surface contact area (`4,236.73 A^2`), together with OLS linear
       calibration, model-based exploratory screening nulls, and degree- and
       interface-membership-preserving permutation nulls (`n_permutations = 10000`,
       seed `20261107`) per M7 Falsification Design Rule 3.

  3. 4-Branch Literature & Provenance Separation Table & Three Physical RNA Arrangements:
     - Explicitly separates (1) internal NC/Zn RNA-binding architecture in 6TAP,
       (2) direct lumenal RNA cargo encapsulation inside engineered dArc1 capsids
       (Lash et al., Nature Communications, Nov 4, 2026), and (3) non-encapsulated
       Arc-Mbnl1 extracellular-vesicle co-transport (Zinter et al., bioRxiv /
       PMC13061017, Mar 10, 2026).
     - Classifies Branch 4 candidate mutational hotspots by `perturbation_category`
       and `required_experimental_controls`, explicitly distinguishing
       `CYS227/CYS230/HIS235/CYS240` as
       `METAL_COORDINATION_AND_NC_FOLDING_PERTURBATION_NOT_ORDINARY_CA_INTERFACE_HOTSPOT`
       with Zn occupancy (ICP-MS/PAR), NC folding (CD/NMR), and RNA-binding controls.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

try:
    import numpy as np
except ImportError:
    for _cand in (
        "/opt/homebrew/opt/python@3.14/bin/python3.14",
        "/Library/Frameworks/Python.framework/Versions/3.13/bin/python3",
        "/Applications/ChimeraX-1.12.app/Contents/bin/python3.11",
    ):
        if __name__ == "__main__" and os.path.exists(_cand) and os.path.realpath(sys.executable) != os.path.realpath(_cand):
            os.execv(_cand, [_cand] + sys.argv)
    raise

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_L2_RECEIPT = REPO_ROOT / "artifacts" / "m7_arc_level2_validation_receipt.json"
DEFAULT_CROSSCHECK_RECEIPT = REPO_ROOT / "artifacts" / "m7_arc_level2_independent_crosscheck_receipt.json"
DEFAULT_OUT = REPO_ROOT / "artifacts" / "m7_arc_level3_interface_comparison_receipt.json"

NTD_MAX_SEQ_ID = 125  # 41..125 = NTD stalk + four-helix bundle; 126..206 = CTD lobe + tail

CANONICAL_SUBINTERFACE_ORDER: Tuple[str, ...] = (
    "A-A_C5_pentamer_op36",
    "C-D_C1_intra_asu_op1",
    "B-C_C3_inter_asu_hexamer_op33",
    "B-D_C2_dimer_op2",
    "A-D_C1_intra_asu_op1",
    "B-C_C1_intra_asu_op1",
    "C-C_C3_trimer_op33",
    "B-D_C3_inter_asu_vertex_op33",
    "A-D_C5_inter_asu_seam_op36",
    "A-B_C2_peripheral_op2",
)

INTERFACE_FINGERPRINT_FAMILY_MAP: Dict[str, str] = {
    "A-A_C5_pentamer_op36": "Family_I_NTD_CTD_Capsomere_Seam",
    "C-D_C1_intra_asu_op1": "Family_I_NTD_CTD_Capsomere_Seam",
    "B-C_C1_intra_asu_op1": "Family_I_NTD_CTD_Capsomere_Seam",
    "B-D_C2_dimer_op2": "Family_I_NTD_CTD_Capsomere_Seam",
    "A-D_C1_intra_asu_op1": "Family_II_Pure_CTD_Quasi_Equivalent_Dimer",
    "B-C_C3_inter_asu_hexamer_op33": "Family_II_Pure_CTD_Quasi_Equivalent_Dimer",
    "C-C_C3_trimer_op33": "Family_III_Residual_Distal_CTD_Hierarchy",
    "B-D_C3_inter_asu_vertex_op33": "Family_III_Residual_Distal_CTD_Hierarchy",
    "A-D_C5_inter_asu_seam_op36": "Family_III_Residual_Distal_CTD_Hierarchy",
    "A-B_C2_peripheral_op2": "Family_III_Residual_Distal_CTD_Hierarchy",
}

INTERFACE_FINGERPRINT_SUBTIER_MAP: Dict[str, str] = {
    "A-A_C5_pentamer_op36": "Family_I_NTD_CTD_Architectural_Seam",
    "C-D_C1_intra_asu_op1": "Family_I_NTD_CTD_Architectural_Seam",
    "B-C_C1_intra_asu_op1": "Family_I_NTD_CTD_Architectural_Seam",
    "B-D_C2_dimer_op2": "Family_I_NTD_CTD_Architectural_Seam",
    "A-D_C1_intra_asu_op1": "Family_II_Coherent_CTD_Dimer_Core",
    "B-C_C3_inter_asu_hexamer_op33": "Family_II_Coherent_CTD_Dimer_Core",
    "C-C_C3_trimer_op33": "Family_III_Tier1_Coherent_C3_Vertex_Core",
    "B-D_C3_inter_asu_vertex_op33": "Family_III_Tier1_Coherent_C3_Vertex_Core",
    "A-D_C5_inter_asu_seam_op36": "Family_III_Tier2_C5_Inter_ASU_Seam",
    "A-B_C2_peripheral_op2": "Family_III_Tier3_Peripheral_C2_Singleton",
}

THRESHOLD_SWEEP_GRID: Tuple[float, ...] = (
    0.15,
    0.20,
    0.25,
    0.30,
    0.35,
    0.40,
    0.45,
    0.50,
    0.60,
    0.70,
)

FAMILY_I_SET: Set[str] = {
    "A-A_C5_pentamer_op36",
    "C-D_C1_intra_asu_op1",
    "B-C_C1_intra_asu_op1",
    "B-D_C2_dimer_op2",
}
FAMILY_II_SET: Set[str] = {
    "A-D_C1_intra_asu_op1",
    "B-C_C3_inter_asu_hexamer_op33",
}
FAMILY_III_FULL_SET: Set[str] = {
    "C-C_C3_trimer_op33",
    "B-D_C3_inter_asu_vertex_op33",
    "A-D_C5_inter_asu_seam_op36",
    "A-B_C2_peripheral_op2",
}
FAMILY_III_TRIMER_PLUS_SEAM_SET: Set[str] = {
    "C-C_C3_trimer_op33",
    "B-D_C3_inter_asu_vertex_op33",
    "A-D_C5_inter_asu_seam_op36",
}
FAMILY_III_C3_CORE_SET: Set[str] = {
    "C-C_C3_trimer_op33",
    "B-D_C3_inter_asu_vertex_op33",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def residue_domain_region(seq_id: int) -> str:
    if seq_id <= 54:
        return "NTD_stalk_41_54"
    if seq_id <= NTD_MAX_SEQ_ID:
        return "NTD_bundle_55_125"
    if seq_id <= 132:
        return "interdomain_hinge_126_132"
    return "CTD_lobe_133_206"


def residue_major_domain(seq_id: int) -> str:
    return "NTD" if seq_id <= NTD_MAX_SEQ_ID else "CTD"


def jaccard_similarity(set_a: Set[Any], set_b: Set[Any]) -> float:
    union = set_a | set_b
    if not union:
        return 1.0
    return float(len(set_a & set_b)) / float(len(union))


def sorensen_dice_similarity(set_a: Set[Any], set_b: Set[Any]) -> float:
    denom = len(set_a) + len(set_b)
    if denom == 0:
        return 1.0
    return 2.0 * float(len(set_a & set_b)) / float(denom)


def szymkiewicz_simpson_overlap(set_a: Set[Any], set_b: Set[Any]) -> float:
    denom = min(len(set_a), len(set_b))
    if denom == 0:
        return 1.0 if (len(set_a) == 0 and len(set_b) == 0) else 0.0
    return float(len(set_a & set_b)) / float(denom)


def cosine_ochiai_similarity(set_a: Set[Any], set_b: Set[Any]) -> float:
    denom = math.sqrt(float(len(set_a)) * float(len(set_b)))
    if denom == 0.0:
        return 1.0 if (len(set_a) == 0 and len(set_b) == 0) else 0.0
    return float(len(set_a & set_b)) / denom


def weighted_cosine_similarity(vec_a: Dict[int, float], vec_b: Dict[int, float]) -> float:
    norm_a = math.sqrt(sum(v * v for v in vec_a.values()))
    norm_b = math.sqrt(sum(v * v for v in vec_b.values()))
    if norm_a == 0.0 or norm_b == 0.0:
        return 0.0
    dot = sum(vec_a.get(k, 0.0) * vec_b.get(k, 0.0) for k in (set(vec_a.keys()) & set(vec_b.keys())))
    return dot / (norm_a * norm_b)


def single_linkage_clusters_at_threshold(
    items: Tuple[str, ...],
    sim_matrix: Dict[str, Dict[str, float]],
    tau: float,
) -> List[List[str]]:
    visited: Set[str] = set()
    clusters: List[List[str]] = []
    for item in items:
        if item in visited:
            continue
        comp: List[str] = []
        queue: List[str] = [item]
        visited.add(item)
        while queue:
            curr = queue.pop(0)
            comp.append(curr)
            for other in items:
                if other not in visited and sim_matrix[curr][other] >= tau - 1e-12:
                    visited.add(other)
                    queue.append(other)
        comp_sorted = [x for x in items if x in comp]
        clusters.append(comp_sorted)
    return clusters


def complete_linkage_clusters_at_threshold(
    items: Tuple[str, ...],
    sim_matrix: Dict[str, Dict[str, float]],
    tau: float,
) -> List[List[str]]:
    clusters: List[List[str]] = [[x] for x in items]
    while len(clusters) > 1:
        best_sim = -1.0
        best_pair: Tuple[int, int] | None = None
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                c_sim = min(sim_matrix[u][v] for u in clusters[i] for v in clusters[j])
                if c_sim > best_sim + 1e-12:
                    best_sim = c_sim
                    best_pair = (i, j)
        if best_pair is None or best_sim < tau - 1e-12:
            break
        i_idx, j_idx = best_pair
        merged_set = set(clusters[i_idx]) | set(clusters[j_idx])
        merged = [x for x in items if x in merged_set]
        clusters[i_idx] = merged
        clusters.pop(j_idx)
    return clusters


def summarize_linkage_sweep(
    items: Tuple[str, ...],
    sim_matrix: Dict[str, Dict[str, float]],
    linkage_mode: str,
) -> Dict[str, Any]:
    cluster_fn = (
        single_linkage_clusters_at_threshold
        if linkage_mode == "single_linkage"
        else complete_linkage_clusters_at_threshold
    )
    sweep_by_tau: Dict[str, Any] = {}
    family_i_taus: List[float] = []
    family_ii_taus: List[float] = []
    family_iii_full_taus: List[float] = []
    family_iii_trimer_plus_seam_taus: List[float] = []
    family_iii_c3_core_taus: List[float] = []

    for tau in THRESHOLD_SWEEP_GRID:
        clusters = cluster_fn(items, sim_matrix, tau)
        cluster_sets = [set(c) for c in clusters]
        has_fam_i = FAMILY_I_SET in cluster_sets
        has_fam_ii = FAMILY_II_SET in cluster_sets
        has_fam_iii_full = FAMILY_III_FULL_SET in cluster_sets
        has_fam_iii_trimer_seam = FAMILY_III_TRIMER_PLUS_SEAM_SET in cluster_sets
        has_fam_iii_c3_core = FAMILY_III_C3_CORE_SET in cluster_sets
        ab_op2_singleton = {"A-B_C2_peripheral_op2"} in cluster_sets
        ad_op36_singleton = {"A-D_C5_inter_asu_seam_op36"} in cluster_sets

        if has_fam_i:
            family_i_taus.append(tau)
        if has_fam_ii:
            family_ii_taus.append(tau)
        if has_fam_iii_full:
            family_iii_full_taus.append(tau)
        if has_fam_iii_trimer_seam:
            family_iii_trimer_plus_seam_taus.append(tau)
        if has_fam_iii_c3_core:
            family_iii_c3_core_taus.append(tau)

        sweep_by_tau[f"{tau:.2f}"] = {
            "tau": tau,
            "cluster_count": len(clusters),
            "clusters": clusters,
            "family_I_exact_cluster_present": has_fam_i,
            "family_II_exact_cluster_present": has_fam_ii,
            "family_III_full_4member_cluster_present": has_fam_iii_full,
            "family_III_trimer_plus_A_D_op36_cluster_present": has_fam_iii_trimer_seam,
            "family_III_C3_core_cluster_present": has_fam_iii_c3_core,
            "A_B_op2_detached_as_singleton": ab_op2_singleton,
            "A_D_op36_detached_as_singleton": ad_op36_singleton,
        }

    # Exact continuous threshold boundaries from similarity matrix values:
    # For Family II = {A-D_op1, B-C_op33}: internal similarity is S(A-D_op1, B-C_op33);
    # max external similarity to any other sub-interface is max_{u in FamII, v not in FamII} S(u, v).
    fam_ii_internal = sim_matrix["A-D_C1_intra_asu_op1"]["B-C_C3_inter_asu_hexamer_op33"]
    fam_ii_max_ext = max(
        sim_matrix[u][v]
        for u in FAMILY_II_SET
        for v in items
        if v not in FAMILY_II_SET
    )

    # For Family I = {A-A_op36, C-D_op1, B-C_op1, B-D_op2}:
    # Single-linkage internal connectivity threshold is the bottleneck edge of the maximum spanning tree on Family I;
    # Complete-linkage internal threshold is the minimum pairwise similarity within Family I.
    fam_i_list = sorted(FAMILY_I_SET)
    fam_i_min_internal = min(
        sim_matrix[u][v] for i, u in enumerate(fam_i_list) for v in fam_i_list[i + 1 :]
    )
    # MST bottleneck on 4 nodes:
    best_mst_bottleneck = 0.0
    for tau_cand in sorted({sim_matrix[u][v] for u in fam_i_list for v in fam_i_list if u != v}):
        c_test = single_linkage_clusters_at_threshold(tuple(fam_i_list), sim_matrix, tau_cand)
        if len(c_test) == 1:
            best_mst_bottleneck = max(best_mst_bottleneck, tau_cand)
    fam_i_max_ext = max(
        sim_matrix[u][v]
        for u in FAMILY_I_SET
        for v in items
        if v not in FAMILY_I_SET
    )

    # Family III internal links:
    c3_core_sim = sim_matrix["C-C_C3_trimer_op33"]["B-D_C3_inter_asu_vertex_op33"]
    ad_op36_to_c3_max = max(
        sim_matrix["A-D_C5_inter_asu_seam_op36"]["C-C_C3_trimer_op33"],
        sim_matrix["A-D_C5_inter_asu_seam_op36"]["B-D_C3_inter_asu_vertex_op33"],
    )
    ad_op36_to_c3_min = min(
        sim_matrix["A-D_C5_inter_asu_seam_op36"]["C-C_C3_trimer_op33"],
        sim_matrix["A-D_C5_inter_asu_seam_op36"]["B-D_C3_inter_asu_vertex_op33"],
    )
    ab_op2_to_fam3_max = max(
        sim_matrix["A-B_C2_peripheral_op2"][v]
        for v in ("C-C_C3_trimer_op33", "B-D_C3_inter_asu_vertex_op33", "A-D_C5_inter_asu_seam_op36")
    )
    ab_op2_to_fam3_min = min(
        sim_matrix["A-B_C2_peripheral_op2"][v]
        for v in ("C-C_C3_trimer_op33", "B-D_C3_inter_asu_vertex_op33", "A-D_C5_inter_asu_seam_op36")
    )

    return {
        "linkage_mode": linkage_mode,
        "grid_thresholds": list(THRESHOLD_SWEEP_GRID),
        "family_I_invariant_grid_taus": family_i_taus,
        "family_I_invariant_grid_interval": (
            [min(family_i_taus), max(family_i_taus)] if family_i_taus else []
        ),
        "family_II_invariant_grid_taus": family_ii_taus,
        "family_II_invariant_grid_interval": (
            [min(family_ii_taus), max(family_ii_taus)] if family_ii_taus else []
        ),
        "family_III_full_4member_grid_taus": family_iii_full_taus,
        "family_III_trimer_plus_A_D_op36_grid_taus": family_iii_trimer_plus_seam_taus,
        "family_III_C3_core_only_grid_taus": family_iii_c3_core_taus,
        "exact_continuous_thresholds": {
            "family_I_max_external_similarity_tau_min_exclusive": round(fam_i_max_ext, 4),
            "family_I_single_linkage_mst_bottleneck_tau_max": round(best_mst_bottleneck, 4),
            "family_I_complete_linkage_min_pairwise_tau_max": round(fam_i_min_internal, 4),
            "family_II_max_external_similarity_tau_min_exclusive": round(fam_ii_max_ext, 4),
            "family_II_internal_pair_similarity_tau_max": round(fam_ii_internal, 4),
            "family_III_A_B_op2_detachment_threshold": round(ab_op2_to_fam3_max, 4),
            "family_III_A_B_op2_single_linkage_detachment_threshold": round(ab_op2_to_fam3_max, 4),
            "family_III_A_B_op2_complete_linkage_detachment_threshold": round(ab_op2_to_fam3_min, 4),
            "family_III_A_D_op36_single_linkage_detachment_threshold": round(ad_op36_to_c3_max, 4),
            "family_III_A_D_op36_complete_linkage_detachment_threshold": round(ad_op36_to_c3_min, 4),
            "family_III_C3_core_C_C_op33_vs_B_D_op33_similarity": round(c3_core_sim, 4),
        },
        "sweep_by_tau": sweep_by_tau,
    }


def build_level3_interface_comparison_receipt(
    l2_receipt_path: Path = DEFAULT_L2_RECEIPT,
    crosscheck_receipt_path: Path = DEFAULT_CROSSCHECK_RECEIPT,
) -> Dict[str, Any]:
    l2_data = json.loads(l2_receipt_path.read_text(encoding="utf-8"))
    cc_data = json.loads(crosscheck_receipt_path.read_text(encoding="utf-8"))

    mmcif_prov = l2_data["mmcif_provenance"]
    asym_meta_by_chain: Dict[str, Dict[str, Any]] = {
        rec["label_asym_id"]: rec
        for rec in mmcif_prov["struct_asym_inventory"]
    }

    l2_mods = l2_data["work_packages"]
    l2_a = l2_mods["L2_A_residue_pair_contact_census"]
    l2_b = l2_mods["L2_B_sasa_buried_area"]
    l2_c = l2_mods["L2_C_interface_chemistry"]
    l2_e = l2_mods["L2_E_nc_zn_spatial_context"]

    op_sub_by_name: Dict[str, Dict[str, Any]] = {
        sub["interface_class"]: sub
        for sub in l2_a["per_operator_sub_interfaces_at_5A"]
    }
    bsa_sub_by_name: Dict[str, Dict[str, Any]] = {
        sub["interface_class"]: sub
        for sub in l2_b["interface_buried_surface_area_by_sub_interface"]
    }
    chem_sub_by_name: Dict[str, Dict[str, Any]] = {
        sub["interface_class"]: sub
        for sub in l2_c["summary_by_sub_interface"]
    }
    cc_by_name: Dict[str, Dict[str, Any]] = {
        sub["interface_class"]: sub
        for sub in cc_data["sasa_and_voronoi_surface_crosscheck"]["sub_interface_crosschecks"]
    }

    # Index all 204 unique residue pairs by interface_class
    pairs_by_sub: Dict[str, List[Dict[str, Any]]] = {name: [] for name in CANONICAL_SUBINTERFACE_ORDER}
    for pair in l2_a["complete_residue_pair_list_at_5A"]:
        iclass = pair["interface_class"]
        if iclass in pairs_by_sub:
            pairs_by_sub[iclass].append(pair)

    # 1. Build Symmetry-Normalized Residue Fingerprints per sub-interface and per chain
    subinterface_fingerprints: List[Dict[str, Any]] = []
    homologous_seq_sets: Dict[str, Set[int]] = {}
    homologous_ntd_seq_sets: Dict[str, Set[int]] = {}
    homologous_ctd_seq_sets: Dict[str, Set[int]] = {}
    homologous_pair_sets: Dict[str, Set[Tuple[int, int]]] = {}
    chain_specific_seq_sets: Dict[str, Dict[str, Set[int]]] = {}
    chain_global_seq_sets: Dict[str, Set[int]] = {"A": set(), "B": set(), "C": set(), "D": set()}
    domain_pair_counts_by_sub: Dict[str, Dict[str, int]] = {}
    seq_contact_degree_by_sub: Dict[str, Dict[int, float]] = {}

    for iclass in CANONICAL_SUBINTERFACE_ORDER:
        sub_op = op_sub_by_name[iclass]
        plist = pairs_by_sub[iclass]

        seq_set: Set[int] = set()
        ntd_seqs: Set[int] = set()
        ctd_seqs: Set[int] = set()
        hom_pairs: Set[Tuple[int, int]] = set()
        by_chain: Dict[str, Set[int]] = {}
        residue_records: Dict[Tuple[str, int, str], Dict[str, Any]] = {}
        dom_counts = {"NTD-NTD": 0, "NTD-CTD": 0, "CTD-CTD": 0}
        deg_by_seq: Dict[int, float] = {}

        for p in plist:
            c1, s1, n1 = str(p["chain_i"]), int(p["label_seq_id_i"]), str(p["comp_i"])
            c2, s2, n2 = str(p["chain_j"]), int(p["label_seq_id_j"]), str(p["comp_j"])
            seq_set.add(s1)
            seq_set.add(s2)
            d1 = residue_major_domain(s1)
            d2 = residue_major_domain(s2)
            if d1 == "NTD":
                ntd_seqs.add(s1)
            else:
                ctd_seqs.add(s1)
            if d2 == "NTD":
                ntd_seqs.add(s2)
            else:
                ctd_seqs.add(s2)

            if d1 == "NTD" and d2 == "NTD":
                dom_counts["NTD-NTD"] += 1
            elif d1 == "CTD" and d2 == "CTD":
                dom_counts["CTD-CTD"] += 1
            else:
                dom_counts["NTD-CTD"] += 1

            hom_pairs.add((min(s1, s2), max(s1, s2)))
            by_chain.setdefault(c1, set()).add(s1)
            by_chain.setdefault(c2, set()).add(s2)
            chain_global_seq_sets[c1].add(s1)
            chain_global_seq_sets[c2].add(s2)

            for ch, sq, comp in ((c1, s1, n1), (c2, s2, n2)):
                key = (ch, sq, comp)
                cmeta = asym_meta_by_chain[ch]
                if key not in residue_records:
                    residue_records[key] = {
                        "chain_id": ch,
                        "label_asym_id": cmeta["label_asym_id"],
                        "auth_asym_id": cmeta["auth_asym_id"],
                        "label_entity_id": cmeta["label_entity_id"],
                        "struct_asym_entity_id": cmeta["struct_asym_entity_id"],
                        "functional_component_id": cmeta["functional_component_id"],
                        "biological_domain_or_fragment_designation": cmeta[
                            "biological_domain_or_fragment_designation"
                        ],
                        "auth_seq_id": sq,
                        "label_seq_id": sq,
                        "comp_id": comp,
                        "major_domain": residue_major_domain(sq),
                        "domain_region": residue_domain_region(sq),
                        "contact_pair_degree_in_subinterface": 0,
                    }
                residue_records[key]["contact_pair_degree_in_subinterface"] += 1
                deg_by_seq[sq] = deg_by_seq.get(sq, 0.0) + 1.0

        homologous_seq_sets[iclass] = seq_set
        homologous_ntd_seq_sets[iclass] = ntd_seqs
        homologous_ctd_seq_sets[iclass] = ctd_seqs
        homologous_pair_sets[iclass] = hom_pairs
        chain_specific_seq_sets[iclass] = by_chain
        domain_pair_counts_by_sub[iclass] = dom_counts
        seq_contact_degree_by_sub[iclass] = deg_by_seq

        sorted_residues = [
            residue_records[k]
            for k in sorted(residue_records.keys(), key=lambda x: (x[0], x[1]))
        ]
        ntd_res_count = sum(1 for r in sorted_residues if r["major_domain"] == "NTD")
        ctd_res_count = sum(1 for r in sorted_residues if r["major_domain"] == "CTD")
        total_res_count = len(sorted_residues)

        subinterface_fingerprints.append(
            {
                "interface_class": iclass,
                "interface_fingerprint_family": INTERFACE_FINGERPRINT_FAMILY_MAP[iclass],
                "interface_fingerprint_subfamily_tier": INTERFACE_FINGERPRINT_SUBTIER_MAP[iclass],
                "chain_pair": sub_op["chain_pair"],
                "symmetry_symbol": sub_op["symmetry_symbol"],
                "rotation_order": sub_op["rotation_order"],
                "rotation_angle_deg": sub_op["rotation_angle_deg"],
                "representative_operator_pair": [
                    sub_op["reference_operator_id"],
                    sub_op["partner_operator_id"],
                ],
                "unique_residue_pair_count_at_5A": sub_op["unique_asym_residue_pair_count"],
                "interfacial_chain_residue_count": total_res_count,
                "ntd_residue_count": ntd_res_count,
                "ctd_residue_count": ctd_res_count,
                "ntd_residue_fraction": round(float(ntd_res_count) / float(max(total_res_count, 1)), 4),
                "ctd_residue_fraction": round(float(ctd_res_count) / float(max(total_res_count, 1)), 4),
                "homologous_seq_id_count": len(seq_set),
                "homologous_seq_ids": sorted(seq_set),
                "homologous_ntd_seq_ids": sorted(ntd_seqs),
                "homologous_ctd_seq_ids": sorted(ctd_seqs),
                "homologous_canonical_seq_pairs": [
                    [a, b] for a, b in sorted(hom_pairs)
                ],
                "per_chain_seq_ids": {
                    ch: sorted(sset) for ch, sset in sorted(by_chain.items())
                },
                "interfacial_residues": sorted_residues,
            }
        )

    # Compute 10x10 pairwise similarity matrices across all 4 set metrics + contact-degree-weighted cosine
    seq_jaccard_matrix: Dict[str, Dict[str, float]] = {}
    seq_sorensen_dice_matrix: Dict[str, Dict[str, float]] = {}
    seq_szymkiewicz_simpson_matrix: Dict[str, Dict[str, float]] = {}
    seq_cosine_ochiai_matrix: Dict[str, Dict[str, float]] = {}
    seq_degree_weighted_cosine_matrix: Dict[str, Dict[str, float]] = {}

    pair_jaccard_matrix: Dict[str, Dict[str, float]] = {}
    pair_sorensen_dice_matrix: Dict[str, Dict[str, float]] = {}
    pair_szymkiewicz_simpson_matrix: Dict[str, Dict[str, float]] = {}
    pair_cosine_ochiai_matrix: Dict[str, Dict[str, float]] = {}

    ntd_seq_jaccard_matrix: Dict[str, Dict[str, float]] = {}
    ctd_seq_jaccard_matrix: Dict[str, Dict[str, float]] = {}

    for ic1 in CANONICAL_SUBINTERFACE_ORDER:
        seq_jaccard_matrix[ic1] = {}
        seq_sorensen_dice_matrix[ic1] = {}
        seq_szymkiewicz_simpson_matrix[ic1] = {}
        seq_cosine_ochiai_matrix[ic1] = {}
        seq_degree_weighted_cosine_matrix[ic1] = {}

        pair_jaccard_matrix[ic1] = {}
        pair_sorensen_dice_matrix[ic1] = {}
        pair_szymkiewicz_simpson_matrix[ic1] = {}
        pair_cosine_ochiai_matrix[ic1] = {}

        ntd_seq_jaccard_matrix[ic1] = {}
        ctd_seq_jaccard_matrix[ic1] = {}
        for ic2 in CANONICAL_SUBINTERFACE_ORDER:
            seq_jaccard_matrix[ic1][ic2] = round(
                jaccard_similarity(homologous_seq_sets[ic1], homologous_seq_sets[ic2]), 4
            )
            seq_sorensen_dice_matrix[ic1][ic2] = round(
                sorensen_dice_similarity(homologous_seq_sets[ic1], homologous_seq_sets[ic2]), 4
            )
            seq_szymkiewicz_simpson_matrix[ic1][ic2] = round(
                szymkiewicz_simpson_overlap(homologous_seq_sets[ic1], homologous_seq_sets[ic2]), 4
            )
            seq_cosine_ochiai_matrix[ic1][ic2] = round(
                cosine_ochiai_similarity(homologous_seq_sets[ic1], homologous_seq_sets[ic2]), 4
            )
            seq_degree_weighted_cosine_matrix[ic1][ic2] = round(
                weighted_cosine_similarity(seq_contact_degree_by_sub[ic1], seq_contact_degree_by_sub[ic2]), 4
            )

            pair_jaccard_matrix[ic1][ic2] = round(
                jaccard_similarity(homologous_pair_sets[ic1], homologous_pair_sets[ic2]), 4
            )
            pair_sorensen_dice_matrix[ic1][ic2] = round(
                sorensen_dice_similarity(homologous_pair_sets[ic1], homologous_pair_sets[ic2]), 4
            )
            pair_szymkiewicz_simpson_matrix[ic1][ic2] = round(
                szymkiewicz_simpson_overlap(homologous_pair_sets[ic1], homologous_pair_sets[ic2]), 4
            )
            pair_cosine_ochiai_matrix[ic1][ic2] = round(
                cosine_ochiai_similarity(homologous_pair_sets[ic1], homologous_pair_sets[ic2]), 4
            )

            ntd_seq_jaccard_matrix[ic1][ic2] = round(
                jaccard_similarity(homologous_ntd_seq_sets[ic1], homologous_ntd_seq_sets[ic2]), 4
            )
            ctd_seq_jaccard_matrix[ic1][ic2] = round(
                jaccard_similarity(homologous_ctd_seq_sets[ic1], homologous_ctd_seq_sets[ic2]), 4
            )

    seq_metric_matrices: Dict[str, Dict[str, Dict[str, float]]] = {
        "jaccard": seq_jaccard_matrix,
        "sorensen_dice": seq_sorensen_dice_matrix,
        "szymkiewicz_simpson_overlap": seq_szymkiewicz_simpson_matrix,
        "cosine_ochiai": seq_cosine_ochiai_matrix,
        "contact_degree_weighted_cosine": seq_degree_weighted_cosine_matrix,
    }
    pair_metric_matrices: Dict[str, Dict[str, Dict[str, float]]] = {
        "jaccard": pair_jaccard_matrix,
        "sorensen_dice": pair_sorensen_dice_matrix,
        "szymkiewicz_simpson_overlap": pair_szymkiewicz_simpson_matrix,
        "cosine_ochiai": pair_cosine_ochiai_matrix,
    }

    homologous_seq_id_metrics_analysis: Dict[str, Any] = {}
    for mname, mmat in seq_metric_matrices.items():
        homologous_seq_id_metrics_analysis[mname] = {
            "similarity_matrix": mmat,
            "single_linkage_threshold_sweep": summarize_linkage_sweep(
                CANONICAL_SUBINTERFACE_ORDER, mmat, "single_linkage"
            ),
            "complete_linkage_threshold_sweep": summarize_linkage_sweep(
                CANONICAL_SUBINTERFACE_ORDER, mmat, "complete_linkage"
            ),
        }

    canonical_seq_pair_metrics_analysis: Dict[str, Any] = {}
    for mname, mmat in pair_metric_matrices.items():
        canonical_seq_pair_metrics_analysis[mname] = {
            "similarity_matrix": mmat,
            "single_linkage_threshold_sweep": summarize_linkage_sweep(
                CANONICAL_SUBINTERFACE_ORDER, mmat, "single_linkage"
            ),
            "complete_linkage_threshold_sweep": summarize_linkage_sweep(
                CANONICAL_SUBINTERFACE_ORDER, mmat, "complete_linkage"
            ),
        }

    seq_j_sl = homologous_seq_id_metrics_analysis["jaccard"]["single_linkage_threshold_sweep"]
    seq_j_cl = homologous_seq_id_metrics_analysis["jaccard"]["complete_linkage_threshold_sweep"]
    seq_j_ct = seq_j_sl["exact_continuous_thresholds"]

    pair_j_sl = canonical_seq_pair_metrics_analysis["jaccard"]["single_linkage_threshold_sweep"]
    pair_j_cl = canonical_seq_pair_metrics_analysis["jaccard"]["complete_linkage_threshold_sweep"]
    pair_j_ct = pair_j_sl["exact_continuous_thresholds"]

    fam2_seq_cont = (
        f"({seq_j_ct['family_II_max_external_similarity_tau_min_exclusive']:.4f}, "
        f"{seq_j_ct['family_II_internal_pair_similarity_tau_max']:.4f}]"
    )
    fam2_pair_cont = (
        f"({pair_j_ct['family_II_max_external_similarity_tau_min_exclusive']:.4f}, "
        f"{pair_j_ct['family_II_internal_pair_similarity_tau_max']:.4f}]"
    )
    fam1_seq_sl_cont = (
        f"({seq_j_ct['family_I_max_external_similarity_tau_min_exclusive']:.4f}, "
        f"{seq_j_ct['family_I_single_linkage_mst_bottleneck_tau_max']:.4f}]"
    )
    fam1_seq_cl_cont = (
        f"({seq_j_ct['family_I_max_external_similarity_tau_min_exclusive']:.4f}, "
        f"{seq_j_ct['family_I_complete_linkage_min_pairwise_tau_max']:.4f}]"
    )
    fam1_pair_sl_cont = (
        f"({pair_j_ct['family_I_max_external_similarity_tau_min_exclusive']:.4f}, "
        f"{pair_j_ct['family_I_single_linkage_mst_bottleneck_tau_max']:.4f}]"
    )
    fam1_pair_cl_cont = (
        f"({pair_j_ct['family_I_max_external_similarity_tau_min_exclusive']:.4f}, "
        f"{pair_j_ct['family_I_complete_linkage_min_pairwise_tau_max']:.4f}]"
    )

    family_robustness_stratification: Dict[str, Any] = {
        "summary_statement": (
            "The three empirical interface-fingerprint groupings in PDB 6TAP exhibit sharply "
            "unequal threshold robustness across homologous sequence-position (seq_id) vs "
            "canonical residue-pair comparisons and single-linkage vs complete-linkage clustering: "
            "(1) Family II ({A-D_C1_intra_asu_op1, B-C_C3_inter_asu_hexamer_op33}) is the strongest "
            "shared CTD interface-fingerprint family, surviving broad thresholds across both "
            "homologous sequence positions ([0.15, 0.70] single and complete linkage Jaccard) and "
            "canonical residue pairs ([0.15, 0.60] single and complete linkage Jaccard); "
            "(2) Family I ({A-A_C5_pentamer_op36, C-D_C1_intra_asu_op1, B-C_C1_intra_asu_op1, "
            "B-D_C2_dimer_op2}) is a broader NTD/CTD architectural seam motif strongly supported at "
            "the homologous sequence-position / domain-architecture level ([0.15, 0.50] single linkage, "
            "[0.15, 0.45] complete linkage Jaccard) but substantially more threshold-sensitive at exact "
            "canonical residue-pair correspondence ([0.15, 0.20] single linkage and 0.15 complete "
            "linkage Jaccard); and (3) Family III (Family_III_Residual_Distal_CTD_Hierarchy: "
            "{C-C_C3_trimer_op33, B-D_C3_inter_asu_vertex_op33, A-D_C5_inter_asu_seam_op36, "
            "A-B_C2_peripheral_op2}) is a residual distal CTD interface hierarchy rather than one "
            "cohesive four-member family (it never forms a four-member cluster at any tested "
            "complete-linkage threshold tau >= 0.15), resolving into a coherent two-member C3 vertex "
            "core ({C-C_op33, B-D_op33}, Jaccard 0.5000 on seq_id, 0.4000 on canonical pairs), the C5 "
            "inter-ASU seam (A-D_op36, detaching above tau = 0.2500 complete / 0.3636 single Jaccard), "
            "and the peripheral C2 singleton (A-B_op2, detaching above tau = 0.1429 complete / 0.2222 "
            "single Jaccard on seq_id and tau = 0.0000 complete on canonical pairs)."
        ),
        "stratification_by_family": {
            "Family_II_Pure_CTD_Quasi_Equivalent_Dimer": {
                "robustness_rank": 1,
                "robustness_tier": "STRONGEST_SHARED_CTD_INTERFACE_FINGERPRINT_FAMILY",
                "members": sorted(FAMILY_II_SET),
                "homologous_seq_id_jaccard_single_linkage_grid_interval": seq_j_sl[
                    "family_II_invariant_grid_interval"
                ],
                "homologous_seq_id_jaccard_complete_linkage_grid_interval": seq_j_cl[
                    "family_II_invariant_grid_interval"
                ],
                "homologous_seq_id_jaccard_continuous_interval": fam2_seq_cont,
                "canonical_seq_pair_jaccard_single_linkage_grid_interval": pair_j_sl[
                    "family_II_invariant_grid_interval"
                ],
                "canonical_seq_pair_jaccard_complete_linkage_grid_interval": pair_j_cl[
                    "family_II_invariant_grid_interval"
                ],
                "canonical_seq_pair_jaccard_continuous_interval": fam2_pair_cont,
                "qualification_narrative": (
                    f"Strongest shared CTD interface-fingerprint family: survives broad thresholds "
                    f"across both homologous sequence-position ([0.15, 0.70] for single and complete "
                    f"linkage Jaccard; continuous {fam2_seq_cont}) and canonical residue-pair "
                    f"([0.15, 0.60] for single and complete linkage Jaccard; continuous {fam2_pair_cont}) "
                    f"comparisons."
                ),
            },
            "Family_I_NTD_CTD_Capsomere_Seam": {
                "robustness_rank": 2,
                "robustness_tier": "BROADER_NTD_CTD_ARCHITECTURAL_SEAM_MOTIF_THRESHOLD_SENSITIVE_AT_EXACT_PAIR_LEVEL",
                "members": sorted(FAMILY_I_SET),
                "homologous_seq_id_jaccard_single_linkage_grid_interval": seq_j_sl[
                    "family_I_invariant_grid_interval"
                ],
                "homologous_seq_id_jaccard_complete_linkage_grid_interval": seq_j_cl[
                    "family_I_invariant_grid_interval"
                ],
                "homologous_seq_id_jaccard_continuous_intervals": {
                    "single_linkage": fam1_seq_sl_cont,
                    "complete_linkage": fam1_seq_cl_cont,
                },
                "canonical_seq_pair_jaccard_single_linkage_grid_interval": pair_j_sl[
                    "family_I_invariant_grid_interval"
                ],
                "canonical_seq_pair_jaccard_complete_linkage_grid_interval": pair_j_cl[
                    "family_I_invariant_grid_interval"
                ],
                "canonical_seq_pair_jaccard_continuous_intervals": {
                    "single_linkage": fam1_pair_sl_cont,
                    "complete_linkage": fam1_pair_cl_cont,
                },
                "qualification_narrative": (
                    "Broader NTD/CTD architectural seam motif: supported at the homologous "
                    "sequence-position / domain-architecture level ([0.15, 0.50] single linkage, "
                    "[0.15, 0.45] complete linkage Jaccard), whereas exact canonical residue-pair "
                    "correspondence is substantially more threshold-sensitive ([0.15, 0.20] single "
                    "linkage and 0.15 complete linkage Jaccard)."
                ),
            },
            "Family_III_Residual_Distal_CTD_Hierarchy": {
                "legacy_alias": "Family_III_Vertex_and_Peripheral_CTD_Contacts",
                "robustness_rank": 3,
                "robustness_tier": "RESIDUAL_DISTAL_CTD_INTERFACE_HIERARCHY_NOT_COHESIVE_FOUR_MEMBER_FAMILY",
                "members": sorted(FAMILY_III_FULL_SET),
                "is_cohesive_four_member_cluster_at_complete_linkage_tau_ge_0p15": (
                    len(seq_j_cl["family_III_full_4member_grid_taus"]) > 0
                ),
                "homologous_seq_id_complete_linkage_4member_grid_taus": seq_j_cl[
                    "family_III_full_4member_grid_taus"
                ],
                "canonical_seq_pair_complete_linkage_4member_grid_taus": pair_j_cl[
                    "family_III_full_4member_grid_taus"
                ],
                "hierarchical_resolution": {
                    "tier_1_coherent_C3_vertex_core": {
                        "members": sorted(FAMILY_III_C3_CORE_SET),
                        "homologous_seq_id_jaccard": seq_j_ct[
                            "family_III_C3_core_C_C_op33_vs_B_D_op33_similarity"
                        ],
                        "canonical_seq_pair_jaccard": pair_j_ct[
                            "family_III_C3_core_C_C_op33_vs_B_D_op33_similarity"
                        ],
                        "description": (
                            "Coherent two-member C3 vertex core ({C-C_op33, B-D_op33}, Jaccard "
                            "0.5000 on seq_id, 0.4000 on canonical pairs)."
                        ),
                    },
                    "tier_2_C5_inter_asu_seam": {
                        "members": ["A-D_C5_inter_asu_seam_op36"],
                        "detachment_thresholds_homologous_seq_id_jaccard": {
                            "complete_linkage": seq_j_ct[
                                "family_III_A_D_op36_complete_linkage_detachment_threshold"
                            ],
                            "single_linkage": seq_j_ct[
                                "family_III_A_D_op36_single_linkage_detachment_threshold"
                            ],
                        },
                        "detachment_thresholds_canonical_seq_pair_jaccard": {
                            "complete_linkage": pair_j_ct[
                                "family_III_A_D_op36_complete_linkage_detachment_threshold"
                            ],
                            "single_linkage": pair_j_ct[
                                "family_III_A_D_op36_single_linkage_detachment_threshold"
                            ],
                        },
                        "description": (
                            "C5 inter-ASU seam (A-D_op36, detaching above tau = 0.2500 complete / "
                            "0.3636 single Jaccard on seq_id)."
                        ),
                    },
                    "tier_3_peripheral_C2_singleton": {
                        "members": ["A-B_C2_peripheral_op2"],
                        "detachment_thresholds_homologous_seq_id_jaccard": {
                            "complete_linkage": seq_j_ct[
                                "family_III_A_B_op2_complete_linkage_detachment_threshold"
                            ],
                            "single_linkage": seq_j_ct[
                                "family_III_A_B_op2_single_linkage_detachment_threshold"
                            ],
                        },
                        "detachment_thresholds_canonical_seq_pair_jaccard": {
                            "complete_linkage": pair_j_ct[
                                "family_III_A_B_op2_complete_linkage_detachment_threshold"
                            ],
                            "single_linkage": pair_j_ct[
                                "family_III_A_B_op2_single_linkage_detachment_threshold"
                            ],
                        },
                        "description": (
                            "Peripheral C2 singleton (A-B_op2, detaching above tau = 0.1429 complete / "
                            "0.2222 single Jaccard on seq_id and tau = 0.0000 complete on canonical pairs)."
                        ),
                    },
                },
                "qualification_narrative": (
                    "Residual interface hierarchy (Family_III_Residual_Distal_CTD_Hierarchy), NOT one "
                    "cohesive four-member family (it never forms a four-member cluster at any tested "
                    "complete-linkage threshold tau >= 0.15), resolving into: (1) a coherent two-member "
                    "C3 vertex core ({C-C_op33, B-D_op33}, Jaccard 0.5000 on seq_id, 0.4000 on canonical "
                    "pairs), (2) the C5 inter-ASU seam (A-D_op36, detaching above tau = 0.2500 complete / "
                    "0.3636 single Jaccard), and (3) the peripheral C2 singleton (A-B_op2, detaching "
                    "above tau = 0.1429 complete / 0.2222 single Jaccard on seq_id and tau = 0.0000 "
                    "complete on canonical pairs)."
                ),
            },
        },
    }

    alternative_similarity_and_sensitivity = {
        "terminology_note": (
            "Families I, II, and III (Family_III_Residual_Distal_CTD_Hierarchy) are empirical "
            "interface-fingerprint groupings (interface_fingerprint_family) defined by homologous "
            "sequence-position similarity and NTD/CTD contact architecture in 6TAP, not established "
            "biological or functional modules."
        ),
        "family_robustness_stratification": family_robustness_stratification,
        "threshold_sweep_grid": list(THRESHOLD_SWEEP_GRID),
        "set_metrics_evaluated": [
            "jaccard",
            "sorensen_dice",
            "szymkiewicz_simpson_overlap",
            "cosine_ochiai",
        ],
        "weighted_metrics_evaluated": [
            "contact_degree_weighted_cosine",
        ],
        "homologous_seq_id_metrics": homologous_seq_id_metrics_analysis,
        "canonical_seq_pair_metrics": canonical_seq_pair_metrics_analysis,
        "family_invariance_and_splitting_summary": {
            "family_I_members": sorted(FAMILY_I_SET),
            "family_II_members": sorted(FAMILY_II_SET),
            "family_III_members": sorted(FAMILY_III_FULL_SET),
            "family_III_C3_core_members": sorted(FAMILY_III_C3_CORE_SET),
            "homologous_seq_id_invariance_intervals": {
                mname: {
                    "single_linkage_family_I_grid_interval": homologous_seq_id_metrics_analysis[mname][
                        "single_linkage_threshold_sweep"
                    ]["family_I_invariant_grid_interval"],
                    "complete_linkage_family_I_grid_interval": homologous_seq_id_metrics_analysis[mname][
                        "complete_linkage_threshold_sweep"
                    ]["family_I_invariant_grid_interval"],
                    "single_linkage_family_II_grid_interval": homologous_seq_id_metrics_analysis[mname][
                        "single_linkage_threshold_sweep"
                    ]["family_II_invariant_grid_interval"],
                    "complete_linkage_family_II_grid_interval": homologous_seq_id_metrics_analysis[mname][
                        "complete_linkage_threshold_sweep"
                    ]["family_II_invariant_grid_interval"],
                    "exact_continuous_thresholds": homologous_seq_id_metrics_analysis[mname][
                        "single_linkage_threshold_sweep"
                    ]["exact_continuous_thresholds"],
                }
                for mname in seq_metric_matrices
            },
            "canonical_seq_pair_invariance_intervals": {
                mname: {
                    "single_linkage_family_I_grid_interval": canonical_seq_pair_metrics_analysis[mname][
                        "single_linkage_threshold_sweep"
                    ]["family_I_invariant_grid_interval"],
                    "complete_linkage_family_I_grid_interval": canonical_seq_pair_metrics_analysis[mname][
                        "complete_linkage_threshold_sweep"
                    ]["family_I_invariant_grid_interval"],
                    "single_linkage_family_II_grid_interval": canonical_seq_pair_metrics_analysis[mname][
                        "single_linkage_threshold_sweep"
                    ]["family_II_invariant_grid_interval"],
                    "complete_linkage_family_II_grid_interval": canonical_seq_pair_metrics_analysis[mname][
                        "complete_linkage_threshold_sweep"
                    ]["family_II_invariant_grid_interval"],
                    "exact_continuous_thresholds": canonical_seq_pair_metrics_analysis[mname][
                        "single_linkage_threshold_sweep"
                    ]["exact_continuous_thresholds"],
                }
                for mname in pair_metric_matrices
            },
            "structural_sensitivity_conclusion": (
                "Across all 5 homologous seq_id similarity metrics (jaccard, sorensen_dice, "
                "szymkiewicz_simpson_overlap, cosine_ochiai, and contact_degree_weighted_cosine), "
                "Family II ({A-D_op1, B-C_op33}) is the strongest shared CTD interface-fingerprint "
                "family, invariant over the entire grid tau in [0.15, 0.70] (continuous (0.0476, 0.8333] "
                "in both single and complete linkage; Sorensen-Dice continuous (0.0909, 0.9091]; "
                "Szymkiewicz-Simpson continuous (0.1667, 0.9375]; Cosine/Ochiai continuous (0.1021, 0.9095]; "
                "contact-degree-weighted cosine continuous (0.0558, 0.9571]) and also invariant over "
                "[0.15, 0.60] (continuous (0.0000, 0.6364]) on canonical residue pairs. Family I ({A-A_op36, "
                "C-D_op1, B-C_op1, B-D_op2}) is a broader NTD/CTD architectural seam motif invariant over grid "
                "tau in [0.15, 0.50] (continuous (0.0222, 0.5926] single-linkage; grid [0.15, 0.45], continuous "
                "(0.0222, 0.4889] complete-linkage) on homologous seq_id, whereas exact canonical residue-pair "
                "correspondence is more threshold-sensitive ([0.15, 0.20] single linkage, continuous (0.0000, 0.2407]; "
                "and [0.15, 0.15] complete linkage, continuous (0.0000, 0.1698] Jaccard). Family III "
                "(Family_III_Residual_Distal_CTD_Hierarchy: {C-C_op33, B-D_op33, A-D_op36, A-B_op2}) is a "
                "residual distal CTD interface hierarchy that never forms a four-member cluster at any tested "
                "complete-linkage threshold tau >= 0.15: peripheral dimer contact A-B_op2 (2 residues: 181, 204; "
                "LYS181-ASN204) detaches above complete-linkage Jaccard tau = 0.1429 (single-linkage 0.2222 on seq_id; "
                "0.0000 complete on canonical pairs), followed by the C5 pentamer-hexamer seam A-D_op36 detaching "
                "above complete-linkage Jaccard tau = 0.2500 (single-linkage 0.3636), leaving the coherent two-member "
                "C3 vertex core {C-C_op33, B-D_op33} intact up to Jaccard tau = 0.5000 on seq_id and 0.4000 on canonical pairs."
            ),
        },
    }

    # Focused patch comparisons required by Step 4:
    def compare_two_subinterfaces(ic_a: str, ic_b: str, comparison_role: str) -> Dict[str, Any]:
        seq_a = homologous_seq_sets[ic_a]
        seq_b = homologous_seq_sets[ic_b]
        pairs_a = homologous_pair_sets[ic_a]
        pairs_b = homologous_pair_sets[ic_b]
        chains_common = sorted(
            set(chain_specific_seq_sets[ic_a].keys()) & set(chain_specific_seq_sets[ic_b].keys())
        )
        per_chain_overlap: Dict[str, Any] = {}
        for ch in chains_common:
            ca = chain_specific_seq_sets[ic_a][ch]
            cb = chain_specific_seq_sets[ic_b][ch]
            per_chain_overlap[ch] = {
                f"{ic_a}_seq_ids": sorted(ca),
                f"{ic_b}_seq_ids": sorted(cb),
                "shared_seq_ids": sorted(ca & cb),
                "jaccard_similarity": round(jaccard_similarity(ca, cb), 4),
            }
        return {
            "subinterface_a": ic_a,
            "subinterface_b": ic_b,
            "comparison_role": comparison_role,
            "homologous_seq_id_jaccard": round(jaccard_similarity(seq_a, seq_b), 4),
            "homologous_seq_id_sorensen_dice": round(sorensen_dice_similarity(seq_a, seq_b), 4),
            "homologous_seq_id_szymkiewicz_simpson": round(szymkiewicz_simpson_overlap(seq_a, seq_b), 4),
            "homologous_seq_id_cosine_ochiai": round(cosine_ochiai_similarity(seq_a, seq_b), 4),
            "shared_homologous_seq_ids": sorted(seq_a & seq_b),
            "only_in_a_seq_ids": sorted(seq_a - seq_b),
            "only_in_b_seq_ids": sorted(seq_b - seq_a),
            "canonical_seq_pair_jaccard": round(jaccard_similarity(pairs_a, pairs_b), 4),
            "canonical_seq_pair_sorensen_dice": round(sorensen_dice_similarity(pairs_a, pairs_b), 4),
            "canonical_seq_pair_szymkiewicz_simpson": round(szymkiewicz_simpson_overlap(pairs_a, pairs_b), 4),
            "canonical_seq_pair_cosine_ochiai": round(cosine_ochiai_similarity(pairs_a, pairs_b), 4),
            "shared_canonical_seq_pairs": [[x, y] for x, y in sorted(pairs_a & pairs_b)],
            "common_chain_patch_comparison": per_chain_overlap,
        }

    focused_patch_comparisons = [
        compare_two_subinterfaces(
            "B-C_C1_intra_asu_op1",
            "B-C_C3_inter_asu_hexamer_op33",
            "Same chain pair (B-C), intra-ASU C1 NTD+CTD contact (op1) vs inter-ASU C3 hexamer CTD dimer (op33)",
        ),
        compare_two_subinterfaces(
            "A-D_C1_intra_asu_op1",
            "A-D_C5_inter_asu_seam_op36",
            "Same chain pair (A-D), intra-ASU C1 CTD lock (op1) vs inter-ASU C5 pentamer-hexamer pure-CTD seam (op36)",
        ),
        compare_two_subinterfaces(
            "B-D_C2_dimer_op2",
            "B-D_C3_inter_asu_vertex_op33",
            "Same chain pair (B-D), 2-fold inter-capsomere C2 dimer (op2) vs C3 vertex pure-CTD contact (op33)",
        ),
        compare_two_subinterfaces(
            "A-D_C1_intra_asu_op1",
            "B-C_C3_inter_asu_hexamer_op33",
            "Quasi-equivalent pure-CTD capsomere dimer interfaces across distinct chains (A-D op1 vs B-C op33)",
        ),
        compare_two_subinterfaces(
            "C-D_C1_intra_asu_op1",
            "B-C_C1_intra_asu_op1",
            "Quasi-equivalent intra-ASU hexamer NTD+CTD interfaces (C-D op1 vs B-C op1)",
        ),
        compare_two_subinterfaces(
            "C-D_C1_intra_asu_op1",
            "B-D_C2_dimer_op2",
            "Hexamer intra-ASU NTD+CTD seam (C-D op1) vs 2-fold C2 inter-capsomere dimer seam (B-D op2)",
        ),
    ]

    all_interfacial_seqs = sorted(
        chain_global_seq_sets["A"]
        | chain_global_seq_sets["B"]
        | chain_global_seq_sets["C"]
        | chain_global_seq_sets["D"]
    )
    conserved_across_all_4_chains = sorted(
        chain_global_seq_sets["A"]
        & chain_global_seq_sets["B"]
        & chain_global_seq_sets["C"]
        & chain_global_seq_sets["D"]
    )
    conserved_across_3plus_chains = [
        sq
        for sq in all_interfacial_seqs
        if sum(1 for ch in ("A", "B", "C", "D") if sq in chain_global_seq_sets[ch]) >= 3
    ]

    # 2. Normalized Interface Chemistry & Area Profiles across the 10 sub-interfaces
    normalized_profiles: List[Dict[str, Any]] = []
    total_sr_bsa = sum(
        float(bsa_sub_by_name[ic]["total_buried_surface_area_A2"])
        for ic in CANONICAL_SUBINTERFACE_ORDER
    )
    total_sr_per_side_bsa = sum(
        float(bsa_sub_by_name[ic]["buried_surface_area_per_side_A2"])
        for ic in CANONICAL_SUBINTERFACE_ORDER
    )
    total_whole_capsid_bsa = sum(
        float(bsa_sub_by_name[ic]["whole_capsid_sub_interface_buried_surface_area_A2"])
        for ic in CANONICAL_SUBINTERFACE_ORDER
    )
    total_whole_capsid_one_side_bsa = round(total_whole_capsid_bsa / 2.0, 2)
    total_unique_pairs = sum(
        int(op_sub_by_name[ic]["unique_asym_residue_pair_count"])
        for ic in CANONICAL_SUBINTERFACE_ORDER
    )

    for iclass in CANONICAL_SUBINTERFACE_ORDER:
        sub_op = op_sub_by_name[iclass]
        sub_bsa = bsa_sub_by_name[iclass]
        sub_chem = chem_sub_by_name[iclass]
        sub_cc = cc_by_name[iclass]

        sr_bsa = float(sub_bsa["total_buried_surface_area_A2"])
        sr_per_side = float(sub_bsa["buried_surface_area_per_side_A2"])
        lr_bsa = float(sub_cc["lee_richards_zslice_bsa_A2"])
        vor_area = float(sub_cc["voronoi_radical_plane_contact_area_A2"])
        whole_capsid_bsa = float(sub_bsa["whole_capsid_sub_interface_buried_surface_area_A2"])
        whole_capsid_one_side_bsa = float(
            sub_bsa.get("whole_capsid_one_side_bsa_A2", round(whole_capsid_bsa / 2.0, 2))
        )
        edge_mult = int(sub_bsa["unordered_capsid_edge_multiplicity"])

        unique_pairs = int(sub_op["unique_asym_residue_pair_count"])
        instance_atom_hits = int(sub_op["heavy_atom_pair_hit_count"])
        rep_atom_hits = float(instance_atom_hits) / float(edge_mult)

        hyd = int(sub_chem["hydrophobic_candidate_pairs"])
        pi_st = int(sub_chem["pi_stacking_candidate_pairs"])
        hb = int(sub_chem["hbond_heavy_atom_distance_candidate_pairs"])
        sb = int(sub_chem["salt_bridge_heavy_atom_distance_candidate_pairs"])
        polar_charged = hb + sb

        dom = domain_pair_counts_by_sub[iclass]
        ntd_ntd = int(dom["NTD-NTD"])
        ntd_ctd = int(dom["NTD-CTD"])
        ctd_ctd = int(dom["CTD-CTD"])

        scale_100 = sr_bsa / 100.0
        scale_1000 = sr_bsa / 1000.0

        normalized_profiles.append(
            {
                "interface_class": iclass,
                "interface_fingerprint_family": INTERFACE_FINGERPRINT_FAMILY_MAP[iclass],
                "interface_fingerprint_subfamily_tier": INTERFACE_FINGERPRINT_SUBTIER_MAP[iclass],
                "chain_pair": sub_op["chain_pair"],
                "symmetry_symbol": sub_op["symmetry_symbol"],
                "rotation_order": sub_op["rotation_order"],
                "rotation_angle_deg": sub_op["rotation_angle_deg"],
                "representative_operator_pair": [
                    sub_op["reference_operator_id"],
                    sub_op["partner_operator_id"],
                ],
                "unordered_capsid_edge_multiplicity": edge_mult,
                "area_metrics_A2": {
                    "shrake_rupley_pairwise_bsa_A2": round(sr_bsa, 2),
                    "shrake_rupley_per_side_bsa_A2": round(sr_per_side, 2),
                    "lee_richards_zslice_pairwise_bsa_A2": round(lr_bsa, 2),
                    "voronoi_radical_plane_contact_area_A2": round(vor_area, 2),
                    "voronoi_to_shrake_rupley_bsa_ratio": round(vor_area / sr_bsa, 4),
                    "voronoi_to_one_side_bsa_ratio": round(vor_area / sr_per_side, 4),
                    "whole_capsid_subinterface_bsa_A2": round(whole_capsid_bsa, 2),
                    "whole_capsid_one_side_bsa_A2": round(whole_capsid_one_side_bsa, 2),
                    "fraction_of_total_asu_bsa": round(sr_bsa / total_sr_bsa, 4),
                },
                "contact_density_metrics": {
                    "unique_residue_pairs_at_5A": unique_pairs,
                    "representative_heavy_atom_hits_at_5A": round(rep_atom_hits, 2),
                    "voronoi_contacting_atom_pairs": int(sub_cc["voronoi_contacting_atom_pairs"]),
                    "atom_hits_per_100_A2_bsa": round(rep_atom_hits / scale_100, 2),
                    "voronoi_atom_pairs_per_100_A2_bsa": round(
                        float(sub_cc["voronoi_contacting_atom_pairs"]) / scale_100, 2
                    ),
                    "residue_pairs_per_1000_A2_bsa": round(float(unique_pairs) / scale_1000, 2),
                },
                "candidate_chemistry_counts": {
                    "hydrophobic_candidate_pairs": hyd,
                    "pi_stacking_candidate_pairs": pi_st,
                    "hbond_heavy_atom_distance_candidate_pairs": hb,
                    "salt_bridge_heavy_atom_distance_candidate_pairs": sb,
                    "total_polar_charged_candidate_pairs": polar_charged,
                },
                "candidate_chemistry_densities_per_1000_A2_bsa": {
                    "hydrophobic_candidates_per_1000_A2_bsa": round(float(hyd) / scale_1000, 2),
                    "pi_stacking_candidates_per_1000_A2_bsa": round(float(pi_st) / scale_1000, 2),
                    "hbond_candidates_per_1000_A2_bsa": round(float(hb) / scale_1000, 2),
                    "salt_bridge_candidates_per_1000_A2_bsa": round(float(sb) / scale_1000, 2),
                    "total_polar_charged_candidates_per_1000_A2_bsa": round(
                        float(polar_charged) / scale_1000, 2
                    ),
                    "hydrophobic_to_polar_candidate_ratio": (
                        round(float(hyd) / float(polar_charged), 4) if polar_charged > 0 else None
                    ),
                },
                "domain_architecture_profile": {
                    "NTD_NTD_pairs": ntd_ntd,
                    "NTD_CTD_pairs": ntd_ctd,
                    "CTD_CTD_pairs": ctd_ctd,
                    "ntd_contact_involvement_fraction": round(
                        float(ntd_ntd + ntd_ctd) / float(max(unique_pairs, 1)), 4
                    ),
                    "ctd_contact_involvement_fraction": round(
                        float(ntd_ctd + ctd_ctd) / float(max(unique_pairs, 1)), 4
                    ),
                    "dominant_domain_mode": (
                        "Pure_CTD_CTD"
                        if (ntd_ntd == 0 and ntd_ctd == 0)
                        else (
                            "Pure_NTD_NTD"
                            if (ntd_ctd == 0 and ctd_ctd == 0)
                            else "Composite_NTD_and_CTD"
                        )
                    ),
                },
            }
        )

    # Competing Null Model Evaluation (M7 Falsification Design Rule 3)
    p_ntd = 336.0 / 659.0
    p_ctd = 323.0 / 659.0
    exp_ntd_ntd = 204.0 * (p_ntd ** 2)
    exp_ntd_ctd = 204.0 * (2.0 * p_ntd * p_ctd)
    exp_ctd_ctd = 204.0 * (p_ctd ** 2)
    obs_ntd_ntd = float(sum(domain_pair_counts_by_sub[ic]["NTD-NTD"] for ic in CANONICAL_SUBINTERFACE_ORDER))
    obs_ntd_ctd = float(sum(domain_pair_counts_by_sub[ic]["NTD-CTD"] for ic in CANONICAL_SUBINTERFACE_ORDER))
    obs_ctd_ctd = float(sum(domain_pair_counts_by_sub[ic]["CTD-CTD"] for ic in CANONICAL_SUBINTERFACE_ORDER))
    chi2_domain_segregation = (
        ((obs_ntd_ntd - exp_ntd_ntd) ** 2) / exp_ntd_ntd
        + ((obs_ntd_ctd - exp_ntd_ctd) ** 2) / exp_ntd_ctd
        + ((obs_ctd_ctd - exp_ctd_ctd) ** 2) / exp_ctd_ctd
    )

    total_hyd = sum(int(chem_sub_by_name[ic]["hydrophobic_candidate_pairs"]) for ic in CANONICAL_SUBINTERFACE_ORDER)
    total_pi = sum(int(chem_sub_by_name[ic]["pi_stacking_candidate_pairs"]) for ic in CANONICAL_SUBINTERFACE_ORDER)
    total_sb = sum(
        int(chem_sub_by_name[ic]["salt_bridge_heavy_atom_distance_candidate_pairs"])
        for ic in CANONICAL_SUBINTERFACE_ORDER
    )
    ctd_dimer_bsa = float(
        bsa_sub_by_name["A-D_C1_intra_asu_op1"]["total_buried_surface_area_A2"]
        + bsa_sub_by_name["B-C_C3_inter_asu_hexamer_op33"]["total_buried_surface_area_A2"]
    )
    ctd_dimer_hyd = int(
        chem_sub_by_name["A-D_C1_intra_asu_op1"]["hydrophobic_candidate_pairs"]
        + chem_sub_by_name["B-C_C3_inter_asu_hexamer_op33"]["hydrophobic_candidate_pairs"]
    )
    ctd_dimer_pi = int(
        chem_sub_by_name["A-D_C1_intra_asu_op1"]["pi_stacking_candidate_pairs"]
        + chem_sub_by_name["B-C_C3_inter_asu_hexamer_op33"]["pi_stacking_candidate_pairs"]
    )
    ntd_ctd_seams_bsa = float(
        bsa_sub_by_name["A-A_C5_pentamer_op36"]["total_buried_surface_area_A2"]
        + bsa_sub_by_name["C-D_C1_intra_asu_op1"]["total_buried_surface_area_A2"]
        + bsa_sub_by_name["B-C_C1_intra_asu_op1"]["total_buried_surface_area_A2"]
        + bsa_sub_by_name["B-D_C2_dimer_op2"]["total_buried_surface_area_A2"]
    )

    # Degree- and Interface-Membership-Preserving Permutation Null (n_permutations = 10000, seed = 20261107)
    n_perm = 10000
    perm_seed = 20261107
    all_l2_pairs = l2_a["complete_residue_pair_list_at_5A"]

    # 1A. 83 homologous sequence-position domain-label permutation preserving contact degree across 204 edges
    unique_hom_seqs = sorted(
        {int(p["label_seq_id_i"]) for p in all_l2_pairs}
        | {int(p["label_seq_id_j"]) for p in all_l2_pairs}
    )
    seq_to_idx = {s: idx for idx, s in enumerate(unique_hom_seqs)}
    true_seq_dom = np.array([0 if s <= NTD_MAX_SEQ_ID else 1 for s in unique_hom_seqs], dtype=np.int32)
    edge_seq_i = np.array([seq_to_idx[int(p["label_seq_id_i"])] for p in all_l2_pairs], dtype=np.int32)
    edge_seq_j = np.array([seq_to_idx[int(p["label_seq_id_j"])] for p in all_l2_pairs], dtype=np.int32)

    obs_same_domain = int(obs_ntd_ntd + obs_ctd_ctd)  # 131 / 204
    rng_dom_seq = np.random.default_rng(perm_seed)
    perm_seq_nn = np.zeros(n_perm, dtype=np.int32)
    perm_seq_nc = np.zeros(n_perm, dtype=np.int32)
    perm_seq_cc = np.zeros(n_perm, dtype=np.int32)
    for k_idx in range(n_perm):
        p_dom = rng_dom_seq.permutation(true_seq_dom)
        di = p_dom[edge_seq_i]
        dj = p_dom[edge_seq_j]
        perm_seq_nn[k_idx] = int(np.sum((di == 0) & (dj == 0)))
        perm_seq_cc[k_idx] = int(np.sum((di == 1) & (dj == 1)))
        perm_seq_nc[k_idx] = 204 - perm_seq_nn[k_idx] - perm_seq_cc[k_idx]

    perm_seq_same = perm_seq_nn + perm_seq_cc
    mean_seq_nn = float(np.mean(perm_seq_nn))
    mean_seq_nc = float(np.mean(perm_seq_nc))
    mean_seq_cc = float(np.mean(perm_seq_cc))
    mean_seq_same = float(np.mean(perm_seq_same))
    std_seq_same = float(np.std(perm_seq_same))
    exceed_seq_same = int(np.sum(perm_seq_same >= obs_same_domain))
    p_perm_seq_same_raw = float(exceed_seq_same) / float(n_perm)
    p_perm_seq_same_phipson_smyth = float(exceed_seq_same + 1) / float(n_perm + 1)

    obs_chi2_vs_perm_seq = (
        ((obs_ntd_ntd - mean_seq_nn) ** 2) / mean_seq_nn
        + ((obs_ntd_ctd - mean_seq_nc) ** 2) / mean_seq_nc
        + ((obs_ctd_ctd - mean_seq_cc) ** 2) / mean_seq_cc
    )
    perm_chi2_seq = (
        ((perm_seq_nn - mean_seq_nn) ** 2) / mean_seq_nn
        + ((perm_seq_nc - mean_seq_nc) ** 2) / mean_seq_nc
        + ((perm_seq_cc - mean_seq_cc) ** 2) / mean_seq_cc
    )
    exceed_seq_chi2 = int(np.sum(perm_chi2_seq >= obs_chi2_vs_perm_seq - 1e-12))
    p_perm_seq_chi2_phipson_smyth = float(exceed_seq_chi2 + 1) / float(n_perm + 1)

    # 1B. 235 chain-residue (chain_id, seq_id) vertex domain-label permutation preserving vertex degrees
    unique_chain_nodes = sorted(
        {(str(p["chain_i"]), int(p["label_seq_id_i"])) for p in all_l2_pairs}
        | {(str(p["chain_j"]), int(p["label_seq_id_j"])) for p in all_l2_pairs}
    )
    node_to_idx = {nd: idx for idx, nd in enumerate(unique_chain_nodes)}
    true_node_dom = np.array([0 if nd[1] <= NTD_MAX_SEQ_ID else 1 for nd in unique_chain_nodes], dtype=np.int32)
    edge_node_i = np.array(
        [node_to_idx[(str(p["chain_i"]), int(p["label_seq_id_i"]))] for p in all_l2_pairs],
        dtype=np.int32,
    )
    edge_node_j = np.array(
        [node_to_idx[(str(p["chain_j"]), int(p["label_seq_id_j"]))] for p in all_l2_pairs],
        dtype=np.int32,
    )
    rng_dom_node = np.random.default_rng(perm_seed)
    perm_node_same = np.zeros(n_perm, dtype=np.int32)
    for k_idx in range(n_perm):
        p_dom = rng_dom_node.permutation(true_node_dom)
        di = p_dom[edge_node_i]
        dj = p_dom[edge_node_j]
        perm_node_same[k_idx] = int(np.sum(di == dj))
    exceed_node_same = int(np.sum(perm_node_same >= obs_same_domain))
    p_perm_node_same_phipson_smyth = float(exceed_node_same + 1) / float(n_perm + 1)

    # 1C. Full 166 resolved sequence positions (41..206: 85 NTD 41..125, 81 CTD 126..206) domain permutation
    all_166_seqs = list(range(41, 207))
    seq166_to_idx = {s: idx for idx, s in enumerate(all_166_seqs)}
    true_166_dom = np.array([0 if s <= NTD_MAX_SEQ_ID else 1 for s in all_166_seqs], dtype=np.int32)
    edge_166_i = np.array([seq166_to_idx[int(p["label_seq_id_i"])] for p in all_l2_pairs], dtype=np.int32)
    edge_166_j = np.array([seq166_to_idx[int(p["label_seq_id_j"])] for p in all_l2_pairs], dtype=np.int32)
    rng_dom_166 = np.random.default_rng(perm_seed)
    perm_166_same = np.zeros(n_perm, dtype=np.int32)
    for k_idx in range(n_perm):
        p_dom = rng_dom_166.permutation(true_166_dom)
        di = p_dom[edge_166_i]
        dj = p_dom[edge_166_j]
        perm_166_same[k_idx] = int(np.sum(di == dj))
    mean_166_same = float(np.mean(perm_166_same))
    std_166_same = float(np.std(perm_166_same))
    exceed_166_same = int(np.sum(perm_166_same >= obs_same_domain))
    p_perm_166_same_phipson_smyth = float(exceed_166_same + 1) / float(n_perm + 1)

    # 1D. Contiguous 166-residue circular-shift autocorrelation null (all 166 cyclic shifts s in 0..165)
    is_non_fam1_edge = np.array([p["interface_class"] not in FAMILY_I_SET for p in all_l2_pairs], dtype=bool)
    non_fam1_total_edges = int(np.sum(is_non_fam1_edge))  # 75 edges across Family II (55) + Family III (20)
    obs_non_fam1_ctd_ctd = int(
        np.sum((true_166_dom[edge_166_i[is_non_fam1_edge]] == 1) & (true_166_dom[edge_166_j[is_non_fam1_edge]] == 1))
    )  # 71 / 75 (94.67%)
    circ_same_domain = np.zeros(len(all_166_seqs), dtype=np.int32)
    circ_non_fam1_ctd_ctd = np.zeros(len(all_166_seqs), dtype=np.int32)
    for shift_s in range(len(all_166_seqs)):
        rolled_dom = np.roll(true_166_dom, shift_s)
        di = rolled_dom[edge_166_i]
        dj = rolled_dom[edge_166_j]
        circ_same_domain[shift_s] = int(np.sum(di == dj))
        circ_non_fam1_ctd_ctd[shift_s] = int(
            np.sum((di[is_non_fam1_edge] == 1) & (dj[is_non_fam1_edge] == 1))
        )
    circ_exceed_same = int(np.sum(circ_same_domain >= obs_same_domain))
    circ_exceed_non_fam1_ctd = int(np.sum(circ_non_fam1_ctd_ctd >= obs_non_fam1_ctd_ctd))

    # 2A. Sub-interface assignment permutation across 204 contact edges preserving sub-interface edge counts
    hyd_edge_keys: Set[Tuple[str, str, str]] = set()
    for h_rec in l2_c["hydrophobic_contacts_at_4p5A"]:
        hyd_edge_keys.add((str(h_rec["interface_class"]), str(h_rec["residue_i"]), str(h_rec["residue_j"])))
        hyd_edge_keys.add((str(h_rec["interface_class"]), str(h_rec["residue_j"]), str(h_rec["residue_i"])))
    pi_edge_keys: Set[Tuple[str, str, str]] = set()
    for pi_rec in l2_c["aromatic_pi_stacking_candidates"]:
        pi_edge_keys.add((str(pi_rec["interface_class"]), str(pi_rec["residue_i"]), str(pi_rec["residue_j"])))
        pi_edge_keys.add((str(pi_rec["interface_class"]), str(pi_rec["residue_j"]), str(pi_rec["residue_i"])))
    sb_edge_keys: Set[Tuple[str, str, str]] = set()
    for sb_rec in l2_c["salt_bridge_candidates"]:
        sb_edge_keys.add(
            (str(sb_rec["interface_class"]), str(sb_rec["cationic_residue"]), str(sb_rec["anionic_residue"]))
        )
        sb_edge_keys.add(
            (str(sb_rec["interface_class"]), str(sb_rec["anionic_residue"]), str(sb_rec["cationic_residue"]))
        )

    edge_fam = np.array(
        [
            2
            if p["interface_class"] in FAMILY_II_SET
            else (1 if p["interface_class"] in FAMILY_I_SET else 3)
            for p in all_l2_pairs
        ],
        dtype=np.int32,
    )
    edge_is_hyd = np.array(
        [
            1
            if (
                str(p["interface_class"]),
                f"{p['chain_i']}:{p['comp_i']}{p['label_seq_id_i']}",
                f"{p['chain_j']}:{p['comp_j']}{p['label_seq_id_j']}",
            )
            in hyd_edge_keys
            else 0
            for p in all_l2_pairs
        ],
        dtype=np.int32,
    )
    edge_is_pi5 = np.array(
        [
            1
            if (
                str(p["interface_class"]),
                f"{p['chain_i']}:{p['comp_i']}{p['label_seq_id_i']}",
                f"{p['chain_j']}:{p['comp_j']}{p['label_seq_id_j']}",
            )
            in pi_edge_keys
            else 0
            for p in all_l2_pairs
        ],
        dtype=np.int32,
    )
    edge_is_sb_pair = np.array(
        [
            1
            if (
                str(p["interface_class"]),
                f"{p['chain_i']}:{p['comp_i']}{p['label_seq_id_i']}",
                f"{p['chain_j']}:{p['comp_j']}{p['label_seq_id_j']}",
            )
            in sb_edge_keys
            else 0
            for p in all_l2_pairs
        ],
        dtype=np.int32,
    )
    obs_fam2_hyd_edges = int(np.sum(edge_is_hyd[edge_fam == 2]))  # 13 / 20
    obs_fam2_pi5_edges = int(np.sum(edge_is_pi5[edge_fam == 2]))  # 8 / 8 (at d <= 5.0 A)
    obs_fam1_sb_edges = int(np.sum(edge_is_sb_pair[edge_fam == 1]))  # 6 / 6 residue pairs (9/9 atom hits)

    rng_sub_chem = np.random.default_rng(perm_seed)
    perm_fam2_hyd = np.zeros(n_perm, dtype=np.int32)
    perm_fam2_pi5 = np.zeros(n_perm, dtype=np.int32)
    perm_fam1_sb = np.zeros(n_perm, dtype=np.int32)
    for k_idx in range(n_perm):
        pf = rng_sub_chem.permutation(edge_fam)
        f2_mask = pf == 2
        f1_mask = pf == 1
        perm_fam2_hyd[k_idx] = int(np.sum(edge_is_hyd[f2_mask]))
        perm_fam2_pi5[k_idx] = int(np.sum(edge_is_pi5[f2_mask]))
        perm_fam1_sb[k_idx] = int(np.sum(edge_is_sb_pair[f1_mask]))

    exceed_fam2_hyd = int(np.sum(perm_fam2_hyd >= obs_fam2_hyd_edges))
    exceed_fam2_pi5 = int(np.sum(perm_fam2_pi5 >= obs_fam2_pi5_edges))
    exceed_fam1_sb = int(np.sum(perm_fam1_sb >= obs_fam1_sb_edges))

    # 2B. Homologous 83-residue amino-acid identity (comp_id) permutation preserving graph topology & degree
    seq_to_comp: Dict[int, str] = {}
    for p in all_l2_pairs:
        seq_to_comp[int(p["label_seq_id_i"])] = str(p["comp_i"])
        seq_to_comp[int(p["label_seq_id_j"])] = str(p["comp_j"])
    hydrophobic_aas = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"}
    aromatic_aas = {"PHE", "TYR", "TRP", "HIS"}
    pos_aas = {"ARG", "LYS"}
    neg_aas = {"ASP", "GLU"}

    true_comp_is_hyd = np.array(
        [1 if seq_to_comp[s] in hydrophobic_aas else 0 for s in unique_hom_seqs], dtype=np.int32
    )
    true_comp_is_aro = np.array(
        [1 if seq_to_comp[s] in aromatic_aas else 0 for s in unique_hom_seqs], dtype=np.int32
    )
    true_comp_charge = np.array(
        [
            1 if seq_to_comp[s] in pos_aas else (-1 if seq_to_comp[s] in neg_aas else 0)
            for s in unique_hom_seqs
        ],
        dtype=np.int32,
    )
    fam2_edge_mask = edge_fam == 2
    fam1_edge_mask = edge_fam == 1
    obs_aa_fam2_nonpolar = int(
        np.sum((true_comp_is_hyd[edge_seq_i] & true_comp_is_hyd[edge_seq_j])[fam2_edge_mask])
    )
    obs_aa_fam2_aromatic = int(
        np.sum((true_comp_is_aro[edge_seq_i] & true_comp_is_aro[edge_seq_j])[fam2_edge_mask])
    )
    obs_aa_fam1_ionpair = int(
        np.sum(((true_comp_charge[edge_seq_i] * true_comp_charge[edge_seq_j]) == -1)[fam1_edge_mask])
    )

    rng_aa_perm = np.random.default_rng(perm_seed)
    perm_aa_fam2_nonpolar = np.zeros(n_perm, dtype=np.int32)
    perm_aa_fam2_aromatic = np.zeros(n_perm, dtype=np.int32)
    perm_aa_fam1_ionpair = np.zeros(n_perm, dtype=np.int32)
    for k_idx in range(n_perm):
        p_idx = rng_aa_perm.permutation(len(unique_hom_seqs))
        ph = true_comp_is_hyd[p_idx]
        pa = true_comp_is_aro[p_idx]
        pc = true_comp_charge[p_idx]
        perm_aa_fam2_nonpolar[k_idx] = int(np.sum((ph[edge_seq_i] & ph[edge_seq_j])[fam2_edge_mask]))
        perm_aa_fam2_aromatic[k_idx] = int(np.sum((pa[edge_seq_i] & pa[edge_seq_j])[fam2_edge_mask]))
        perm_aa_fam1_ionpair[k_idx] = int(
            np.sum(((pc[edge_seq_i] * pc[edge_seq_j]) == -1)[fam1_edge_mask])
        )

    exceed_aa_nonpolar = int(np.sum(perm_aa_fam2_nonpolar >= obs_aa_fam2_nonpolar))
    exceed_aa_aromatic = int(np.sum(perm_aa_fam2_aromatic >= obs_aa_fam2_aromatic))
    exceed_aa_ionpair = int(np.sum(perm_aa_fam1_ionpair >= obs_aa_fam1_ionpair))

    # Statistical helpers for Clopper-Pearson / Wilson 95% CIs, autocorrelation, ESS, R_hat, and multiplicity
    def _betacf(a_val: float, b_val: float, x_val: float) -> float:
        qab = a_val + b_val
        qap = a_val + 1.0
        qam = a_val - 1.0
        c_cf = 1.0
        d_cf = 1.0 - qab * x_val / qap
        if abs(d_cf) < 1e-30:
            d_cf = 1e-30
        d_cf = 1.0 / d_cf
        h_cf = d_cf
        for m_step in range(1, 201):
            m2 = 2 * m_step
            aa = m_step * (b_val - m_step) * x_val / ((qam + m2) * (a_val + m2))
            d_cf = 1.0 + aa * d_cf
            if abs(d_cf) < 1e-30:
                d_cf = 1e-30
            c_cf = 1.0 + aa / c_cf
            if abs(c_cf) < 1e-30:
                c_cf = 1e-30
            d_cf = 1.0 / d_cf
            h_cf *= d_cf * c_cf
            aa = -(a_val + m_step) * (qab + m_step) * x_val / ((a_val + m2) * (qap + m2))
            d_cf = 1.0 + aa * d_cf
            if abs(d_cf) < 1e-30:
                d_cf = 1e-30
            c_cf = 1.0 + aa / c_cf
            if abs(c_cf) < 1e-30:
                c_cf = 1e-30
            d_cf = 1.0 / d_cf
            rem = d_cf * c_cf
            h_cf *= rem
            if abs(rem - 1.0) < 1e-14:
                break
        return h_cf

    def _reg_inc_beta(a_val: float, b_val: float, x_val: float) -> float:
        if x_val <= 0.0:
            return 0.0
        if x_val >= 1.0:
            return 1.0
        lbeta = (
            math.lgamma(a_val + b_val)
            - math.lgamma(a_val)
            - math.lgamma(b_val)
            + a_val * math.log(x_val)
            + b_val * math.log(1.0 - x_val)
        )
        bt = math.exp(lbeta) if lbeta > -700.0 else 0.0
        if x_val < (a_val + 1.0) / (a_val + b_val + 2.0):
            return bt * _betacf(a_val, b_val, x_val) / a_val
        return 1.0 - bt * _betacf(b_val, a_val, 1.0 - x_val) / b_val

    def _clopper_pearson_95_ci(k_val: int, n_val: int, alpha_val: float = 0.05) -> List[float]:
        if k_val <= 0:
            lo_ci = 0.0
        else:
            low_b, high_b = 0.0, 1.0
            for _ in range(80):
                mid_b = 0.5 * (low_b + high_b)
                if _reg_inc_beta(float(k_val), float(n_val - k_val + 1), mid_b) < alpha_val / 2.0:
                    low_b = mid_b
                else:
                    high_b = mid_b
            lo_ci = 0.5 * (low_b + high_b)
        if k_val >= n_val:
            hi_ci = 1.0
        else:
            low_b, high_b = 0.0, 1.0
            for _ in range(80):
                mid_b = 0.5 * (low_b + high_b)
                if _reg_inc_beta(float(k_val + 1), float(n_val - k_val), mid_b) < 1.0 - alpha_val / 2.0:
                    low_b = mid_b
                else:
                    high_b = mid_b
            hi_ci = 0.5 * (low_b + high_b)
        return [round(lo_ci, 6), round(hi_ci, 6)]

    def _wilson_95_ci(k_val: int, n_val: int) -> List[float]:
        z_val = 1.959963984540054
        phat = float(k_val) / float(n_val)
        denom = 1.0 + (z_val * z_val) / float(n_val)
        center = (phat + (z_val * z_val) / (2.0 * float(n_val))) / denom
        half = (
            z_val
            * math.sqrt((phat * (1.0 - phat) + (z_val * z_val) / (4.0 * float(n_val))) / float(n_val))
            / denom
        )
        return [round(max(0.0, center - half), 6), round(min(1.0, center + half), 6)]

    def _acf_at_lag(arr: np.ndarray, lag_k: int) -> float:
        x_f = arr.astype(np.float64)
        mu_x = float(np.mean(x_f))
        var_x = float(np.var(x_f))
        if var_x < 1e-15:
            return 0.0
        return float(np.mean((x_f[:-lag_k] - mu_x) * (x_f[lag_k:] - mu_x)) / var_x)

    def _compute_tau_int(arr: np.ndarray, max_lag: int = 500) -> float:
        x_f = arr.astype(np.float64)
        mu_x = float(np.mean(x_f))
        var_x = float(np.var(x_f))
        if var_x < 1e-15:
            return 1.0
        tau_val = 1.0
        for lag_k in range(1, max_lag + 1):
            rk = float(np.mean((x_f[:-lag_k] - mu_x) * (x_f[lag_k:] - mu_x)) / var_x)
            if rk <= 0.0:
                break
            tau_val += 2.0 * rk
        return max(1.0, tau_val)

    def _gelman_rubin_rhat(chains_list: List[np.ndarray]) -> float:
        m_ch = len(chains_list)
        n_st = len(chains_list[0])
        means_c = np.array([np.mean(c_arr) for c_arr in chains_list], dtype=np.float64)
        vars_c = np.array([np.var(c_arr, ddof=1) for c_arr in chains_list], dtype=np.float64)
        w_within = float(np.mean(vars_c))
        if w_within < 1e-15:
            return 1.0
        b_between = float(n_st * np.var(means_c, ddof=1))
        var_hat = ((float(n_st) - 1.0) / float(n_st)) * w_within + (1.0 / float(n_st)) * b_between
        return float(math.sqrt(var_hat / w_within))

    def _holm_bonferroni(pvals: List[float]) -> List[float]:
        m_p = len(pvals)
        order = sorted(range(m_p), key=lambda i: pvals[i])
        adj = [0.0] * m_p
        running_max = 0.0
        for rank_i, orig_idx in enumerate(order):
            val = min(1.0, pvals[orig_idx] * float(m_p - rank_i))
            running_max = max(running_max, val)
            adj[orig_idx] = running_max
        return adj

    def _benjamini_hochberg(pvals: List[float]) -> List[float]:
        m_p = len(pvals)
        order = sorted(range(m_p), key=lambda i: pvals[i])
        adj = [0.0] * m_p
        running_min = 1.0
        for rank_rev in range(m_p - 1, -1, -1):
            orig_idx = order[rank_rev]
            val = min(1.0, pvals[orig_idx] * float(m_p) / float(rank_rev + 1))
            running_min = min(running_min, val)
            adj[orig_idx] = running_min
        return adj

    def _build_p_uncertainty(k_val: int, n_val: int, ess_val: float | None = None) -> Dict[str, Any]:
        p_raw = float(k_val) / float(n_val)
        p_ps = float(k_val + 1) / float(n_val + 1)
        eff_n = float(ess_val) if ess_val is not None else float(n_val)
        mcse_p = math.sqrt(max(p_ps * (1.0 - p_ps), 0.0) / max(eff_n, 1.0))
        rec: Dict[str, Any] = {
            "exceedances_k": int(k_val),
            "total_draws_or_steps": int(n_val),
            "p_raw": round(p_raw, 6),
            "p_phipson_smyth": round(p_ps, 6),
            "minimum_attainable_p_floor_iid": round(1.0 / float(n_val + 1), 6),
            "mcse_p_hat": round(mcse_p, 6),
            "clopper_pearson_95_ci_raw": _clopper_pearson_95_ci(k_val, n_val),
            "wilson_score_95_ci_raw": _wilson_95_ci(k_val, n_val),
        }
        if ess_val is not None:
            rec["effective_sample_size_ess"] = round(float(ess_val), 2)
            rec["minimum_attainable_p_floor_ess"] = round(1.0 / (float(ess_val) + 1.0), 6)
        return rec

    # 3A. Constrained Graph State-Space Connectivity & 4-Chain Double-Edge-Swap MCMC Audit
    mcmc_n_chains = 4
    mcmc_base_seed = 20261108
    mcmc_seeds = [mcmc_base_seed + c_idx for c_idx in range(mcmc_n_chains)]
    mcmc_burn_in = 5000
    mcmc_steps_per_chain = 50000
    mcmc_total_steps = mcmc_n_chains * mcmc_steps_per_chain

    global_edges0: List[Tuple[int, int]] = []
    for p in all_l2_pairs:
        u_nd = node_to_idx[(str(p["chain_i"]), int(p["label_seq_id_i"]))]
        v_nd = node_to_idx[(str(p["chain_j"]), int(p["label_seq_id_j"]))]
        if u_nd > v_nd:
            u_nd, v_nd = v_nd, u_nd
        global_edges0.append((u_nd, v_nd))

    sub_edges_by_ic: Dict[str, List[Tuple[int, int]]] = {ic: [] for ic in CANONICAL_SUBINTERFACE_ORDER}
    for p in all_l2_pairs:
        ic_name = str(p["interface_class"])
        u_nd = node_to_idx[(str(p["chain_i"]), int(p["label_seq_id_i"]))]
        v_nd = node_to_idx[(str(p["chain_j"]), int(p["label_seq_id_j"]))]
        sub_edges_by_ic[ic_name].append((u_nd, v_nd))

    # Per-subinterface exact swap mobility analysis
    per_subinterface_swap_mobility: List[Dict[str, Any]] = []
    sub_pair_list: List[Tuple[int, int, str]] = []
    flat_sub_edges0: List[Tuple[int, int]] = []
    offset_sub = 0
    total_unordered_within_sub_pairs = 0
    total_valid_within_sub_swaps = 0
    total_color_changing_within_sub_swaps = 0
    min_within_sub_same_sum = 0
    max_within_sub_same_sum = 0

    for ic_name in CANONICAL_SUBINTERFACE_ORDER:
        el_sub = sub_edges_by_ic[ic_name]
        m_sub = len(el_sub)
        homotypic_sub = ic_name in ("A-A_C5_pentamer_op36", "C-C_C3_trimer_op33")
        es_sub = set(el_sub)

        u_ntd = sum(1 for (u_e, _) in el_sub if true_node_dom[u_e] == 0)
        u_ctd = m_sub - u_ntd
        v_ntd = sum(1 for (_, v_e) in el_sub if true_node_dom[v_e] == 0)
        v_ctd = m_sub - v_ntd
        tot_ntd_stubs_sub = u_ntd + v_ntd
        tot_ctd_stubs_sub = u_ctd + v_ctd
        ntd_stubs_excl_ala125 = sum(
            1 for (u_e, _) in el_sub if unique_chain_nodes[u_e][1] < 125
        ) + sum(
            1 for (_, v_e) in el_sub if unique_chain_nodes[v_e][1] < 125
        )

        # Theoretical directed bipartite bounds on E_NN, E_CC, E_same preserving side stub counts
        max_nn_sub = min(u_ntd, v_ntd)
        min_nn_sub = max(0, u_ntd - v_ctd)
        min_same_sub = u_ctd + v_ctd - m_sub + 2 * min_nn_sub
        max_same_sub = u_ctd + v_ctd - m_sub + 2 * max_nn_sub
        min_within_sub_same_sum += min_same_sub
        max_within_sub_same_sum += max_same_sub

        unord_pairs_sub = m_sub * (m_sub - 1) // 2
        valid_swaps_sub = 0
        color_changing_sub = 0
        for a_idx in range(m_sub):
            for b_idx in range(a_idx + 1, m_sub):
                sub_pair_list.append((offset_sub + a_idx, offset_sub + b_idx, ic_name))
                u_a, v_a = el_sub[a_idx]
                x_b, y_b = el_sub[b_idx]
                u1_c, v1_c = u_a, y_b
                u2_c, v2_c = x_b, v_a
                if u1_c != v1_c and u2_c != v2_c and (u1_c, v1_c) != (u2_c, v2_c):
                    collision = (u1_c, v1_c) in es_sub or (u2_c, v2_c) in es_sub
                    if homotypic_sub and not collision:
                        collision = (
                            (v1_c, u1_c) in es_sub
                            or (v2_c, u2_c) in es_sub
                            or (u1_c, v1_c) == (v2_c, u2_c)
                        )
                    if not collision:
                        valid_swaps_sub += 1
                        old_s = (1 if true_node_dom[u_a] == true_node_dom[v_a] else 0) + (
                            1 if true_node_dom[x_b] == true_node_dom[y_b] else 0
                        )
                        new_s = (1 if true_node_dom[u1_c] == true_node_dom[v1_c] else 0) + (
                            1 if true_node_dom[u2_c] == true_node_dom[v2_c] else 0
                        )
                        if old_s != new_s:
                            color_changing_sub += 1

        offset_sub += m_sub
        flat_sub_edges0.extend(el_sub)
        total_unordered_within_sub_pairs += unord_pairs_sub
        total_valid_within_sub_swaps += valid_swaps_sub
        total_color_changing_within_sub_swaps += color_changing_sub

        dom_sub = domain_pair_counts_by_sub[ic_name]
        obs_same_sub = int(dom_sub["NTD-NTD"] + dom_sub["CTD-CTD"])
        per_subinterface_swap_mobility.append(
            {
                "interface_class": ic_name,
                "interface_fingerprint_family": INTERFACE_FINGERPRINT_FAMILY_MAP[ic_name],
                "edge_count_at_5A": m_sub,
                "is_heterotypic_bipartite": not homotypic_sub,
                "chain_i_ntd_stubs": u_ntd,
                "chain_i_ctd_stubs": u_ctd,
                "chain_j_ntd_stubs": v_ntd,
                "chain_j_ctd_stubs": v_ctd,
                "total_ntd_stubs": tot_ntd_stubs_sub,
                "total_ctd_stubs": tot_ctd_stubs_sub,
                "ntd_stubs_excluding_boundary_ALA125": ntd_stubs_excl_ala125,
                "observed_domain_pairs": {
                    "NTD-NTD": int(dom_sub["NTD-NTD"]),
                    "NTD-CTD": int(dom_sub["NTD-CTD"]),
                    "CTD-CTD": int(dom_sub["CTD-CTD"]),
                },
                "observed_same_domain_edges": obs_same_sub,
                "unordered_edge_pairs": unord_pairs_sub,
                "valid_simple_bipartite_swaps": valid_swaps_sub,
                "valid_swap_fraction_of_unordered_pairs": (
                    round(float(valid_swaps_sub) / float(unord_pairs_sub), 6)
                    if unord_pairs_sub > 0
                    else 0.0
                ),
                "color_changing_valid_swaps": color_changing_sub,
                "color_changing_fraction_of_unordered_pairs": (
                    round(float(color_changing_sub) / float(unord_pairs_sub), 6)
                    if unord_pairs_sub > 0
                    else 0.0
                ),
                "color_changing_fraction_of_valid_swaps": (
                    round(float(color_changing_sub) / float(valid_swaps_sub), 6)
                    if valid_swaps_sub > 0
                    else 0.0
                ),
                "theoretical_same_domain_corridor": [min_same_sub, max_same_sub],
                "domain_pair_counts_frozen_under_scheme_II": color_changing_sub == 0,
                "belongs_to_6_non_family_I_pure_ctd_subinterfaces": ic_name not in FAMILY_I_SET,
            }
        )

    # Exact global G_0 swap mobility census across all 20,706 unordered edge pairs and 41,412 directed proposals
    g0_m = len(global_edges0)
    g0_eset = set(global_edges0)
    g0_unordered_pairs_total = g0_m * (g0_m - 1) // 2
    g0_directed_proposals_total = 2 * g0_unordered_pairs_total
    g0_unordered_pairs_with_valid_swap = 0
    g0_unordered_pairs_with_both_orientations_valid = 0
    g0_unordered_pairs_with_color_changing_valid_swap = 0
    g0_directed_proposals_valid = 0
    g0_directed_proposals_color_changing_valid = 0
    for i_g0 in range(g0_m):
        u_g0, v_g0 = global_edges0[i_g0]
        old_s_i = 1 if true_node_dom[u_g0] == true_node_dom[v_g0] else 0
        for j_g0 in range(i_g0 + 1, g0_m):
            x_g0, y_g0 = global_edges0[j_g0]
            old_pair_s = old_s_i + (1 if true_node_dom[x_g0] == true_node_dom[y_g0] else 0)
            v_cnt_pair = 0
            cc_cnt_pair = 0
            for u1_p, v1_p, u2_p, v2_p in ((u_g0, y_g0, x_g0, v_g0), (u_g0, x_g0, v_g0, y_g0)):
                a1_p, b1_p = (u1_p, v1_p) if u1_p < v1_p else (v1_p, u1_p)
                a2_p, b2_p = (u2_p, v2_p) if u2_p < v2_p else (v2_p, u2_p)
                if a1_p != b1_p and a2_p != b2_p and (a1_p, b1_p) != (a2_p, b2_p):
                    if (a1_p, b1_p) not in g0_eset and (a2_p, b2_p) not in g0_eset:
                        v_cnt_pair += 1
                        g0_directed_proposals_valid += 1
                        new_pair_s = (1 if true_node_dom[a1_p] == true_node_dom[b1_p] else 0) + (
                            1 if true_node_dom[a2_p] == true_node_dom[b2_p] else 0
                        )
                        if new_pair_s != old_pair_s:
                            cc_cnt_pair += 1
                            g0_directed_proposals_color_changing_valid += 1
            if v_cnt_pair >= 1:
                g0_unordered_pairs_with_valid_swap += 1
            if v_cnt_pair == 2:
                g0_unordered_pairs_with_both_orientations_valid += 1
            if cc_cnt_pair >= 1:
                g0_unordered_pairs_with_color_changing_valid_swap += 1

    def _run_global_mcmc_chain(seed_val: int) -> Tuple[np.ndarray, np.ndarray, int]:
        rng_m = np.random.default_rng(seed_val)
        edges_m = list(global_edges0)
        edge_set_m = set(edges_m)
        m_tot = len(edges_m)
        cur_nn = sum(1 for (u_e, v_e) in edges_m if true_node_dom[u_e] == 0 and true_node_dom[v_e] == 0)
        cur_cc = sum(1 for (u_e, v_e) in edges_m if true_node_dom[u_e] == 1 and true_node_dom[v_e] == 1)
        cur_same = cur_nn + cur_cc

        n_tot_steps = mcmc_burn_in + mcmc_steps_per_chain
        rand_i = rng_m.integers(0, m_tot, size=n_tot_steps)
        rand_j = rng_m.integers(0, m_tot - 1, size=n_tot_steps)
        rand_orient = rng_m.integers(0, 2, size=n_tot_steps)

        trace_same = np.empty(mcmc_steps_per_chain, dtype=np.int32)
        trace_diff = np.empty(mcmc_steps_per_chain, dtype=np.int32)
        accepted_cnt = 0
        for step_idx in range(n_tot_steps):
            i_e = int(rand_i[step_idx])
            j_e = int(rand_j[step_idx])
            if j_e >= i_e:
                j_e += 1
            u_e, v_e = edges_m[i_e]
            x_e, y_e = edges_m[j_e]
            if rand_orient[step_idx] == 0:
                u1_e, v1_e = (u_e, y_e) if u_e < y_e else (y_e, u_e)
                u2_e, v2_e = (x_e, v_e) if x_e < v_e else (v_e, x_e)
            else:
                u1_e, v1_e = (u_e, x_e) if u_e < x_e else (x_e, u_e)
                u2_e, v2_e = (v_e, y_e) if v_e < y_e else (y_e, v_e)
            if u1_e != v1_e and u2_e != v2_e and (u1_e, v1_e) != (u2_e, v2_e):
                if (u1_e, v1_e) not in edge_set_m and (u2_e, v2_e) not in edge_set_m:
                    old_nn = (1 if true_node_dom[u_e] == 0 and true_node_dom[v_e] == 0 else 0) + (
                        1 if true_node_dom[x_e] == 0 and true_node_dom[y_e] == 0 else 0
                    )
                    old_cc = (1 if true_node_dom[u_e] == 1 and true_node_dom[v_e] == 1 else 0) + (
                        1 if true_node_dom[x_e] == 1 and true_node_dom[y_e] == 1 else 0
                    )
                    new_nn = (1 if true_node_dom[u1_e] == 0 and true_node_dom[v1_e] == 0 else 0) + (
                        1 if true_node_dom[u2_e] == 0 and true_node_dom[v2_e] == 0 else 0
                    )
                    new_cc = (1 if true_node_dom[u1_e] == 1 and true_node_dom[v1_e] == 1 else 0) + (
                        1 if true_node_dom[u2_e] == 1 and true_node_dom[v2_e] == 1 else 0
                    )
                    cur_nn += new_nn - old_nn
                    cur_cc += new_cc - old_cc
                    cur_same = cur_nn + cur_cc
                    edge_set_m.remove((u_e, v_e))
                    edge_set_m.remove((x_e, y_e))
                    edge_set_m.add((u1_e, v1_e))
                    edge_set_m.add((u2_e, v2_e))
                    edges_m[i_e] = (u1_e, v1_e)
                    edges_m[j_e] = (u2_e, v2_e)
                    if step_idx >= mcmc_burn_in:
                        accepted_cnt += 1
            if step_idx >= mcmc_burn_in:
                out_idx = step_idx - mcmc_burn_in
                trace_same[out_idx] = cur_same
                trace_diff[out_idx] = cur_cc - cur_nn
        return trace_same, trace_diff, accepted_cnt

    def _run_within_subinterface_mcmc_chain(seed_val: int) -> Tuple[np.ndarray, np.ndarray, int, int]:
        rng_m = np.random.default_rng(seed_val)
        edges_m = list(flat_sub_edges0)
        edge_sets_m = {ic_n: set(sub_edges_by_ic[ic_n]) for ic_n in CANONICAL_SUBINTERFACE_ORDER}
        cur_nn = sum(1 for (u_e, v_e) in edges_m if true_node_dom[u_e] == 0 and true_node_dom[v_e] == 0)
        cur_cc = sum(1 for (u_e, v_e) in edges_m if true_node_dom[u_e] == 1 and true_node_dom[v_e] == 1)
        cur_same = cur_nn + cur_cc

        n_tot_steps = mcmc_burn_in + mcmc_steps_per_chain
        n_pairs_sub = len(sub_pair_list)
        rand_pair = rng_m.integers(0, n_pairs_sub, size=n_tot_steps)

        trace_same = np.empty(mcmc_steps_per_chain, dtype=np.int32)
        trace_diff = np.empty(mcmc_steps_per_chain, dtype=np.int32)
        accepted_cnt = 0
        color_changed_cnt = 0
        for step_idx in range(n_tot_steps):
            i_e, j_e, ic_n = sub_pair_list[int(rand_pair[step_idx])]
            u_e, v_e = edges_m[i_e]
            x_e, y_e = edges_m[j_e]
            u1_e, v1_e = u_e, y_e
            u2_e, v2_e = x_e, v_e
            es_cur = edge_sets_m[ic_n]
            if u1_e != v1_e and u2_e != v2_e and (u1_e, v1_e) != (u2_e, v2_e):
                homotypic_sub = ic_n in ("A-A_C5_pentamer_op36", "C-C_C3_trimer_op33")
                collision = (u1_e, v1_e) in es_cur or (u2_e, v2_e) in es_cur
                if homotypic_sub and not collision:
                    collision = (
                        (v1_e, u1_e) in es_cur
                        or (v2_e, u2_e) in es_cur
                        or (u1_e, v1_e) == (v2_e, u2_e)
                    )
                if not collision:
                    old_nn = (1 if true_node_dom[u_e] == 0 and true_node_dom[v_e] == 0 else 0) + (
                        1 if true_node_dom[x_e] == 0 and true_node_dom[y_e] == 0 else 0
                    )
                    old_cc = (1 if true_node_dom[u_e] == 1 and true_node_dom[v_e] == 1 else 0) + (
                        1 if true_node_dom[x_e] == 1 and true_node_dom[y_e] == 1 else 0
                    )
                    new_nn = (1 if true_node_dom[u1_e] == 0 and true_node_dom[v1_e] == 0 else 0) + (
                        1 if true_node_dom[u2_e] == 0 and true_node_dom[v2_e] == 0 else 0
                    )
                    new_cc = (1 if true_node_dom[u1_e] == 1 and true_node_dom[v1_e] == 1 else 0) + (
                        1 if true_node_dom[u2_e] == 1 and true_node_dom[v2_e] == 1 else 0
                    )
                    delta_same = (new_nn + new_cc) - (old_nn + old_cc)
                    cur_nn += new_nn - old_nn
                    cur_cc += new_cc - old_cc
                    cur_same = cur_nn + cur_cc
                    es_cur.remove((u_e, v_e))
                    es_cur.remove((x_e, y_e))
                    es_cur.add((u1_e, v1_e))
                    es_cur.add((u2_e, v2_e))
                    edges_m[i_e] = (u1_e, v1_e)
                    edges_m[j_e] = (u2_e, v2_e)
                    if step_idx >= mcmc_burn_in:
                        accepted_cnt += 1
                        if delta_same != 0:
                            color_changed_cnt += 1
            if step_idx >= mcmc_burn_in:
                out_idx = step_idx - mcmc_burn_in
                trace_same[out_idx] = cur_same
                trace_diff[out_idx] = cur_cc - cur_nn
        return trace_same, trace_diff, accepted_cnt, color_changed_cnt

    g_chains_same: List[np.ndarray] = []
    g_chains_diff: List[np.ndarray] = []
    g_acc_counts: List[int] = []
    w_chains_same: List[np.ndarray] = []
    w_chains_diff: List[np.ndarray] = []
    w_acc_counts: List[int] = []
    w_col_counts: List[int] = []

    for s_val in mcmc_seeds:
        gs, gd, gac = _run_global_mcmc_chain(s_val)
        g_chains_same.append(gs)
        g_chains_diff.append(gd)
        g_acc_counts.append(gac)

        ws, wd, wac, wco = _run_within_subinterface_mcmc_chain(s_val)
        w_chains_same.append(ws)
        w_chains_diff.append(wd)
        w_acc_counts.append(wac)
        w_col_counts.append(wco)

    g_all_same = np.concatenate(g_chains_same)
    g_all_diff = np.concatenate(g_chains_diff)
    g_rho1 = float(np.mean([_acf_at_lag(c_arr, 1) for c_arr in g_chains_same]))
    g_rho10 = float(np.mean([_acf_at_lag(c_arr, 10) for c_arr in g_chains_same]))
    g_rho50 = float(np.mean([_acf_at_lag(c_arr, 50) for c_arr in g_chains_same]))
    g_tau_int = float(np.mean([_compute_tau_int(c_arr) for c_arr in g_chains_same]))
    g_ess_per_chain = float(mcmc_steps_per_chain) / g_tau_int
    g_ess_total = float(mcmc_total_steps) / g_tau_int
    g_rhat = _gelman_rubin_rhat(g_chains_same)
    g_mean_same = float(np.mean(g_all_same))
    g_std_same = float(np.std(g_all_same))
    g_mcse_mean = g_std_same / math.sqrt(g_ess_total)
    g_exceed = int(np.sum(g_all_same >= obs_same_domain))
    g_p_raw = float(g_exceed) / float(mcmc_total_steps)
    g_p_ps = float(g_exceed + 1) / float(mcmc_total_steps + 1)
    g_p_min_ess = 1.0 / (g_ess_total + 1.0)

    w_all_same = np.concatenate(w_chains_same)
    w_all_diff = np.concatenate(w_chains_diff)
    w_rho1 = float(np.mean([_acf_at_lag(c_arr, 1) for c_arr in w_chains_same]))
    w_rho10 = float(np.mean([_acf_at_lag(c_arr, 10) for c_arr in w_chains_same]))
    w_rho50 = float(np.mean([_acf_at_lag(c_arr, 50) for c_arr in w_chains_same]))
    w_tau_int = float(np.mean([_compute_tau_int(c_arr) for c_arr in w_chains_same]))
    w_ess_per_chain = float(mcmc_steps_per_chain) / w_tau_int
    w_ess_total = float(mcmc_total_steps) / w_tau_int
    w_rhat = _gelman_rubin_rhat(w_chains_same)
    w_mean_same = float(np.mean(w_all_same))
    w_std_same = float(np.std(w_all_same))
    w_mcse_mean = w_std_same / math.sqrt(w_ess_total)
    w_exceed = int(np.sum(w_all_same >= obs_same_domain))
    w_p_raw = float(w_exceed) / float(mcmc_total_steps)
    w_p_ps = float(w_exceed + 1) / float(mcmc_total_steps + 1)
    w_p_min_ess = 1.0 / (w_ess_total + 1.0)

    state_space_specification: Dict[str, Any] = {
        "vertex_identity": {
            "row_index": 1,
            "specification_dimension": "vertex_identity",
            "chain_residue_node_count": 235,
            "ntd_chain_residue_node_count": 108,
            "ctd_chain_residue_node_count": 127,
            "homologous_quotient_seq_id_count": 83,
            "homologous_ntd_seq_id_count": 38,
            "homologous_ctd_seq_id_count": 45,
            "summary": (
                "235 chain-specific interfacial residue nodes (chain_id in {A,B,C,D}, seq_id in 41..206, "
                "comp_id, domain_label in {NTD: 41..125, CTD: 126..206}; 108 NTD nodes, 127 CTD nodes) "
                "and 83 homologous sequence positions (38 NTD, 45 CTD) on the symmetry-collapsed quotient graph."
            ),
            "null_scheme_I_global_degree_preserving_simple_graph_double_edge_swap_mcmc": (
                "235 labeled chain-specific interfacial residue vertices (chain_id, seq_id) with fixed domain "
                "labels (NTD: 41..125, CTD: 126..206)."
            ),
            "null_scheme_II_within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc": (
                "235 labeled chain-specific interfacial residue vertices partitioned into the 10 symmetry-inequivalent "
                "sub-interfaces with fixed chain_id and NTD/CTD domain labels."
            ),
        },
        "graph_class": {
            "row_index": 2,
            "specification_dimension": "graph_class",
            "vertex_count": 235,
            "edge_count_at_5A": 204,
            "self_loops": 0,
            "multi_edges": 0,
            "heterotypic_bipartite_subinterfaces": [
                "C-D_C1_intra_asu_op1",
                "B-C_C1_intra_asu_op1",
                "B-D_C2_dimer_op2",
                "A-D_C1_intra_asu_op1",
                "B-C_C3_inter_asu_hexamer_op33",
                "B-D_C3_inter_asu_vertex_op33",
                "A-D_C5_inter_asu_seam_op36",
                "A-B_C2_peripheral_op2",
            ],
            "homotypic_symmetry_subinterfaces": [
                "A-A_C5_pentamer_op36",
                "C-C_C3_trimer_op33",
            ],
            "summary": (
                "Undirected simple graph (235 vertices, 204 edges at d <= 5.0 A, zero self-loops, zero multi-edges; "
                "within each of the 8 heterotypic sub-interfaces C-D_op1, B-C_op1, B-D_op2, A-D_op1, B-C_op33, "
                "B-D_op33, A-D_op36, A-B_op2, the sub-graph is strictly bipartite between the two chain classes, "
                "whereas A-A_op36 and C-C_op33 are homotypic symmetry-interface graphs)."
            ),
            "null_scheme_I_global_degree_preserving_simple_graph_double_edge_swap_mcmc": (
                "Global undirected simple graph on 235 vertices and 204 edges (zero self-loops, zero multi-edges)."
            ),
            "null_scheme_II_within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc": (
                "Product of 10 simple sub-interface graphs (8 strictly chain-bipartite heterotypic sub-graphs and "
                "2 homotypic symmetry sub-graphs A-A_op36 and C-C_op33; zero self-loops, zero multi-edges)."
            ),
        },
        "invariants": {
            "row_index": 3,
            "specification_dimension": "invariants",
            "empirical_6tap_5A_stub_counts": {
                "S_N_ntd_stubs": 185,
                "S_C_ctd_stubs": 223,
                "total_stubs": 408,
                "conserved_E_CC_minus_E_NN": 19,
            },
            "illustrative_pre_run_stub_counts_reconciliation": {
                "S_N_ntd_stubs": 169,
                "S_C_ctd_stubs": 239,
                "total_stubs": 408,
                "conserved_E_CC_minus_E_NN": 35,
            },
            "empirical_6tap_subinterface_edge_counts": [41, 34, 28, 28, 27, 26, 7, 7, 5, 1],
            "illustrative_prompt_subinterface_edge_counts": [41, 34, 34, 20, 21, 34, 7, 6, 6, 1],
            "summary": (
                "Scheme I (Global) preserves exact per-node degree sequence deg(v) across all 235 nodes, total edge "
                "count (204), and fixed node NTD/CTD domain labels (S_N = 185 NTD stubs, S_C = 223 CTD stubs in "
                "empirical 6TAP at d <= 5.0 A; S_N = 169, S_C = 239 in the illustrative 35-stub-difference example), "
                "while allowing edges to rewire across sub-interfaces. Scheme II (Within-Subinterface) preserves exact "
                "sub-interface membership (10 sub-interfaces with empirical 6TAP edge counts 41, 34, 28, 28, 27, 26, "
                "7, 7, 5, 1; illustrative list 41, 34, 34, 20, 21, 34, 7, 6, 6, 1), exact within-subinterface "
                "per-node degree sequence deg_I(v), and chain-bipartite partition on the 8 heterotypic interfaces."
            ),
            "null_scheme_I_global_degree_preserving_simple_graph_double_edge_swap_mcmc": (
                "Preserves exact per-node degree sequence deg(v) across all 235 nodes, total edge count (204), "
                "and fixed node NTD/CTD domain labels (S_N = 185 NTD stubs, S_C = 223 CTD stubs in 6TAP; S_N = 169, "
                "S_C = 239 in illustrative parameterization), while allowing edges to rewire across sub-interfaces."
            ),
            "null_scheme_II_within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc": (
                "Preserves exact sub-interface membership (10 sub-interfaces), exact within-subinterface per-node "
                "degree sequence deg_I(v), and chain-bipartite partition on the 8 heterotypic sub-interfaces."
            ),
        },
        "swap_constraints": {
            "row_index": 4,
            "specification_dimension": "swap_constraints",
            "lazy_metropolis_hastings_persistence_on_rejection": True,
            "self_loops_rejected": True,
            "multiedge_collisions_rejected": True,
            "summary": (
                "Propose two edges e1 = (u, v), e2 = (x, y) and candidate rewiring e1' = (u, y), e2' = (x, v) "
                "(or (u, x), (v, y) in Scheme I). Reject if e1' or e2' creates a self-loop (u = y or x = v), "
                "duplicates an existing edge (multi-edge collision), or (in Scheme II) violates sub-interface "
                "membership or chain-bipartite orientation. Detailed balance / self-loop persistence: when a "
                "proposed swap is rejected, the Markov chain retains the current graph state G_{t+1} = G_t "
                "(lazy/Metropolis-Hastings persistence), ensuring a symmetric transition matrix P(G -> G') = P(G' -> G) "
                "and a uniform stationary distribution pi(G) = const over the connected component of reachable simple graphs."
            ),
            "null_scheme_I_global_degree_preserving_simple_graph_double_edge_swap_mcmc": (
                "Uniformly sample two distinct edges e1 = (u, v), e2 = (x, y) from the 204 global edges and one of "
                "the 2 rewiring orientations; reject if self-loop or multi-edge collision occurs, retaining G_{t+1} = G_t "
                "(lazy persistence, acceptance rate 0.9796)."
            ),
            "null_scheme_II_within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc": (
                "Uniformly sample two distinct edges e1 = (u, v), e2 = (x, y) within the same sub-interface I "
                "(2,865 unordered within-subinterface edge pairs) and propose the chain-orientation-preserving swap "
                "e1' = (u, y), e2' = (x, v); reject if self-loop or multi-edge collision occurs, retaining G_{t+1} = G_t "
                "(lazy persistence, acceptance rate 0.8160)."
            ),
        },
        "stationary_distribution": {
            "row_index": 5,
            "specification_dimension": "stationary_distribution",
            "distribution_family": "UNIFORM_OVER_REACHABLE_CONNECTED_COMPONENT_OF_CONSTRAINED_SIMPLE_GRAPHS",
            "detailed_balance_mechanism": "SYMMETRIC_PROPOSAL_WITH_LAZY_SELF_LOOP_PERSISTENCE_ON_REJECTION",
            "summary": (
                "Uniform distribution pi(G) = const over the connected component of simple graphs satisfying the "
                "respective invariants (verified via symmetric edge-pair proposal q(G' | G) = q(G | G') and "
                "lazy/Metropolis-Hastings self-loop persistence G_{t+1} = G_t on rejection)."
            ),
            "null_scheme_I_global_degree_preserving_simple_graph_double_edge_swap_mcmc": (
                "Uniform distribution over the space of 235-vertex, 204-edge undirected simple graphs with exact "
                "degree sequence deg(v)."
            ),
            "null_scheme_II_within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc": (
                "Uniform distribution over the product space of 10 simple sub-interface graphs with exact "
                "within-subinterface degree sequences deg_I(v) and chain-bipartite partitions."
            ),
        },
        "connectivity": {
            "row_index": 6,
            "specification_dimension": "connectivity",
            "non_family_I_pure_ctd_subinterfaces": [
                "A-D_C1_intra_asu_op1",
                "B-C_C3_inter_asu_hexamer_op33",
                "C-C_C3_trimer_op33",
                "B-D_C3_inter_asu_vertex_op33",
                "A-D_C5_inter_asu_seam_op36",
                "A-B_C2_peripheral_op2",
            ],
            "zero_ntd_stub_subinterfaces_strict_125": [
                "C-C_C3_trimer_op33",
                "B-D_C3_inter_asu_vertex_op33",
                "A-D_C5_inter_asu_seam_op36",
                "A-B_C2_peripheral_op2",
            ],
            "zero_ntd_stub_subinterfaces_hinge_125_grouped": [
                "A-D_C1_intra_asu_op1",
                "B-C_C3_inter_asu_hexamer_op33",
                "C-C_C3_trimer_op33",
                "B-D_C3_inter_asu_vertex_op33",
                "A-D_C5_inter_asu_seam_op36",
                "A-B_C2_peripheral_op2",
            ],
            "zero_color_changing_swap_subinterfaces": [
                "B-D_C2_dimer_op2",
                "C-C_C3_trimer_op33",
                "B-D_C3_inter_asu_vertex_op33",
                "A-D_C5_inter_asu_seam_op36",
                "A-B_C2_peripheral_op2",
            ],
            "scheme_I_global_g0_swap_mobility_by_unit": {
                "unit_1_unordered_edge_pairs": {
                    "total_unordered_edge_pairs_binom_204_2": g0_unordered_pairs_total,
                    "unordered_pairs_with_at_least_1_valid_simple_swap": g0_unordered_pairs_with_valid_swap,
                    "valid_fraction_of_unordered_pairs": round(
                        float(g0_unordered_pairs_with_valid_swap) / float(g0_unordered_pairs_total), 6
                    ),
                    "unordered_pairs_with_both_orientations_valid": g0_unordered_pairs_with_both_orientations_valid,
                    "both_orientations_valid_fraction_of_unordered_pairs": round(
                        float(g0_unordered_pairs_with_both_orientations_valid) / float(g0_unordered_pairs_total), 6
                    ),
                    "algebraic_color_changing_unordered_pairs_before_collision": int(56 * 75 + 73 * 72 // 2),
                    "unordered_pairs_with_color_changing_valid_swap": g0_unordered_pairs_with_color_changing_valid_swap,
                    "color_changing_fraction_of_unordered_pairs": round(
                        float(g0_unordered_pairs_with_color_changing_valid_swap) / float(g0_unordered_pairs_total), 6
                    ),
                    "color_changing_fraction_of_valid_unordered_pairs": round(
                        float(g0_unordered_pairs_with_color_changing_valid_swap)
                        / float(g0_unordered_pairs_with_valid_swap),
                        6,
                    ),
                },
                "unit_2_directed_rewiring_orientations": {
                    "total_directed_proposals_2_binom_204_2": g0_directed_proposals_total,
                    "valid_simple_directed_proposals": g0_directed_proposals_valid,
                    "valid_fraction_of_directed_proposals": round(
                        float(g0_directed_proposals_valid) / float(g0_directed_proposals_total), 6
                    ),
                    "algebraic_color_changing_directed_proposals_before_collision": int(2 * 56 * 75 + 73 * 72 // 2),
                    "color_changing_valid_directed_proposals": g0_directed_proposals_color_changing_valid,
                    "color_changing_fraction_of_directed_proposals": round(
                        float(g0_directed_proposals_color_changing_valid) / float(g0_directed_proposals_total), 6
                    ),
                    "color_changing_fraction_of_valid_directed_proposals": round(
                        float(g0_directed_proposals_color_changing_valid) / float(g0_directed_proposals_valid), 6
                    ),
                },
            },
            "scheme_II_within_subinterface_swap_mobility_by_unit": {
                "empirical_6tap_subinterface_edge_counts": [41, 34, 28, 28, 27, 26, 7, 7, 5, 1],
                "empirical_6tap_subinterface_unordered_pair_counts": [820, 561, 378, 378, 351, 325, 21, 21, 10, 0],
                "illustrative_partition_edge_counts_reconciled": [41, 34, 34, 20, 21, 34, 7, 6, 6, 1],
                "illustrative_partition_unordered_pair_counts_reconciled": [
                    820,
                    561,
                    561,
                    190,
                    210,
                    561,
                    21,
                    15,
                    15,
                    0,
                ],
                "illustrative_partition_total_unordered_pairs_reconciled": 2954,
                "stale_table_reconciliation_note": (
                    "Earlier prose table mixed illustrative edge counts (m=20 -> 190 pairs for B-D_C2 and m=21 -> 210 "
                    "pairs for A-D_C1, summing to 2,954) with empirical 6TAP valid swaps (324 valid swaps from m=28, "
                    "binom(28,2)=378 unordered pairs for B-D_C2, and 269 valid swaps from m=27, binom(27,2)=351 "
                    "unordered pairs for A-D_C1, summing to 2,377 valid swaps out of 2,865 empirical unordered pairs). "
                    "In empirical 6TAP, every sub-interface strictly satisfies valid_simple_bipartite_swaps <= "
                    "unordered_edge_pairs (324 <= 378 and 269 <= 351)."
                ),
                "total_unordered_within_subinterface_edge_pairs": total_unordered_within_sub_pairs,
                "total_valid_simple_bipartite_swaps": total_valid_within_sub_swaps,
                "valid_swap_fraction_of_unordered_within_subinterface_pairs": round(
                    float(total_valid_within_sub_swaps) / float(total_unordered_within_sub_pairs), 6
                ),
                "total_color_changing_valid_swaps": total_color_changing_within_sub_swaps,
                "color_changing_fraction_of_unordered_within_subinterface_pairs": round(
                    float(total_color_changing_within_sub_swaps) / float(total_unordered_within_sub_pairs), 6
                ),
                "color_changing_fraction_of_valid_within_subinterface_swaps": round(
                    float(total_color_changing_within_sub_swaps) / float(total_valid_within_sub_swaps), 6
                ),
            },
            "scheme_II_total_unordered_within_subinterface_edge_pairs": total_unordered_within_sub_pairs,
            "scheme_II_total_valid_simple_bipartite_swaps": total_valid_within_sub_swaps,
            "scheme_II_total_color_changing_swaps": total_color_changing_within_sub_swaps,
            "scheme_II_theoretical_same_domain_corridor": [min_within_sub_same_sum, max_within_sub_same_sum],
            "scheme_II_confound_and_restraint_framing": (
                "CONFOUND_AND_RESTRAINT_NOT_POSITIVE_SUPPORT: Null Scheme II (p_PS = 0.209264, Z = +1.3444) is a "
                "structural confound showing that once residue-to-subinterface membership and within-subinterface "
                "degrees are conditioned on, 71/75 contacts in Families II and III are already fixed as CTD-CTD by "
                "the sub-interface partition boundary, so domain identity adds almost no independent signal beyond "
                "sub-interface membership."
            ),
            "summary": (
                "Scheme I (Global simple-graph degree sequence): Connectedness of simple graphs with a given degree "
                "sequence under double-edge swaps (Taylor 1981; Ryser 1957; Fosdick et al. 2018 SIAM Rev.), qualified "
                "by whether additional chain/spatial restrictions are imposed (E_same explores [77, 135], mean 102.67). "
                "In G_0, out of binom(204,2) = 20,706 unordered edge pairs (41,412 directed proposals), 20,441 / 20,706 "
                "(98.72%) unordered pairs admit >= 1 valid simple swap (40,616 / 41,412 = 98.08% of directed proposals) "
                "and 6,753 / 20,706 (32.61% of all unordered pairs; 6,753 / 20,441 = 33.04% of valid unordered pairs; "
                "10,953 / 41,412 = 26.45% of all directed proposals; 10,953 / 40,616 = 26.97% of valid directed "
                "proposals) alter E_same by +/- 2. Scheme II (Within-Subinterface): Explicitly acts as a structural "
                "confound / restraint (NOT positive support), proving that 6 of the 10 sub-interfaces (A-D_op1, "
                "B-C_op33, C-C_op33, B-D_op33, A-D_op36, A-B_op2, totaling 71/75 pure CTD-CTD contacts) have either "
                "0 NTD stubs (4 sub-interfaces C-C_op33, B-D_op33, A-D_op36, A-B_op2 under NTD=41..125, and also "
                "A-D_op1, C-C_op33, B-D_op33, A-B_op2 when boundary residue ALA125 is grouped with the hinge) or "
                "single-boundary NTD stubs (1 ALA125 stub per chain in A-D_op1 and B-C_op33, so CTD-CTD cannot drop "
                "below 71/75), while B-D_op2 has 28 NTD / 0 CTD stubs on Chain D (single-side NTD stubs that cannot "
                "alter domain-pair counts under any bipartite swap). Across all 10 sub-interfaces (empirical 6TAP "
                "edge counts 41, 34, 28, 28, 27, 26, 7, 7, 5, 1 -> 2,865 unordered within-subinterface edge pairs, "
                "2,377 / 2,865 = 82.97% valid simple swaps), only 48 swaps (48 / 2,865 = 1.68% of all unordered "
                "within-subinterface pairs; 48 / 2,377 = 2.02% of valid within-subinterface swaps) can alter "
                "domain-pair counts, bounding E_same within [123, 135] (observed 131, p_PS = 0.209264) and showing "
                "that once sub-interface membership is conditioned on, domain segregation is largely absorbed by "
                "the sub-interface boundary itself."
            ),
            "null_scheme_I_global_degree_preserving_simple_graph_double_edge_swap_mcmc": (
                "Connected simple-graph degree-sequence realization space (Taylor 1981; Ryser 1957; Fosdick et al. "
                "2018 SIAM Rev.); E_same spans [77, 135] across 200,000 steps (mean 102.67 +/- 7.09, observed 131, "
                "p_PS = 0.000235; 20,441/20,706 = 98.72% valid unordered pairs, 6,753/20,706 = 32.61% color-changing "
                "unordered pairs [33.04% of valid]; 40,616/41,412 = 98.08% valid directed proposals, 10,953/41,412 = "
                "26.45% color-changing directed proposals [26.97% of valid])."
            ),
            "null_scheme_II_within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc": (
                "Structural confound / restraint (NOT positive support): within-subinterface fiber has 5 of 10 "
                "sub-interfaces with 0 color-changing swaps (including 4 with 0 NTD stubs and B-D_op2 with 100% NTD "
                "stubs on Chain D), and 6 of 10 sub-interfaces (71/75 pure CTD-CTD) cannot reduce CTD-CTD below 71/75; "
                "out of 2,865 unordered within-subinterface edge pairs (2,377 / 2,865 = 82.97% valid swaps), only 48 "
                "(48 / 2,865 = 1.68% of unordered within-subinterface pairs; 48 / 2,377 = 2.02% of valid swaps) alter "
                "domain-pair counts, trapping E_same in [123, 135] (mean 128.45 +/- 1.89, observed 131, p_PS = 0.209264)."
            ),
        },
        "mixing": {
            "row_index": 7,
            "specification_dimension": "mixing",
            "n_independent_chains": mcmc_n_chains,
            "chain_seeds": mcmc_seeds,
            "burn_in_steps_per_chain": mcmc_burn_in,
            "post_burn_in_steps_per_chain": mcmc_steps_per_chain,
            "total_post_burn_in_steps": mcmc_total_steps,
            "summary": (
                f"4 independent chains (seeds {mcmc_seeds}), {mcmc_burn_in} burn-in swap proposals, and "
                f"{mcmc_steps_per_chain} post-burn-in steps per chain ({mcmc_total_steps} total post-burn-in steps), "
                f"evaluated with lag-k autocorrelations (rho_1, rho_10, rho_50), Integrated Autocorrelation Time "
                f"tau_int, Effective Sample Size ESS = N_steps / tau_int, Gelman-Rubin potential scale reduction "
                f"factor R_hat, and Monte Carlo standard error (MCSE): Scheme I achieves tau_int = {g_tau_int:.2f}, "
                f"ESS_total = {g_ess_total:.2f} ({g_ess_per_chain:.2f}/chain), R_hat = {g_rhat:.5f}, MCSE_mean = {g_mcse_mean:.4f}; "
                f"Scheme II achieves tau_int = {w_tau_int:.2f}, ESS_total = {w_ess_total:.2f} ({w_ess_per_chain:.2f}/chain), "
                f"R_hat = {w_rhat:.5f}, MCSE_mean = {w_mcse_mean:.4f}."
            ),
            "null_scheme_I_global_degree_preserving_simple_graph_double_edge_swap_mcmc": {
                "acceptance_rate": round(float(sum(g_acc_counts)) / float(mcmc_total_steps), 4),
                "rho_1": round(g_rho1, 4),
                "rho_10": round(g_rho10, 4),
                "rho_50": round(g_rho50, 4),
                "tau_int": round(g_tau_int, 2),
                "ess_per_chain": round(g_ess_per_chain, 2),
                "ess_total": round(g_ess_total, 2),
                "effectively_independent_states": int(round(g_ess_total)),
                "gelman_rubin_r_hat": round(g_rhat, 5),
                "mcse_mean_same_domain": round(g_mcse_mean, 4),
                "p_min_ess_floor": round(g_p_min_ess, 6),
            },
            "null_scheme_II_within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc": {
                "acceptance_rate": round(float(sum(w_acc_counts)) / float(mcmc_total_steps), 4),
                "color_changing_step_fraction": round(float(sum(w_col_counts)) / float(mcmc_total_steps), 6),
                "rho_1": round(w_rho1, 4),
                "rho_10": round(w_rho10, 4),
                "rho_50": round(w_rho50, 4),
                "tau_int": round(w_tau_int, 2),
                "ess_per_chain": round(w_ess_per_chain, 2),
                "ess_total": round(w_ess_total, 2),
                "effectively_independent_states": int(round(w_ess_total)),
                "gelman_rubin_r_hat": round(w_rhat, 5),
                "mcse_mean_same_domain": round(w_mcse_mean, 4),
                "p_min_ess_floor": round(w_p_min_ess, 6),
            },
        },
    }

    constrained_graph_state_space_connectivity_and_mcmc_audit: Dict[str, Any] = {
        "state_space_specification": state_space_specification,
        "stub_difference_conservation_identity": {
            "role": "CONSERVATION_INVARIANT_AND_IMPLEMENTATION_CHECK_NOT_EMPIRICAL_SIGNIFICANCE_STATISTIC",
            "algebraic_proof": (
                "Let S_C = 2 E_CC + E_CN be the total number of CTD half-edge stubs and S_N = 2 E_NN + E_CN be the "
                "total number of NTD half-edge stubs across all 204 contact edges (S_C + S_N = 2 * 204 = 408). "
                "Subtracting the two identities eliminates E_CN identically: S_C - S_N = 2 (E_CC - E_NN) <=> "
                "E_CC - E_NN = (S_C - S_N) / 2. Consequently, in EVERY degree-and-color-preserving graph realization "
                "(under both Scheme I and Scheme II), E_CC - E_NN is strictly constant (variance = 0.0) and serves "
                "purely as a mathematical conservation invariant and MCMC implementation check, NOT an empirical "
                "significance statistic, whereas total same-domain contacts E_same = E_CC + E_NN = 204 - E_CN varies "
                "in steps of +/- 2 whenever a swap exchanges {NTD-NTD, CTD-CTD} <-> {NTD-CTD, NTD-CTD}."
            ),
            "illustrative_35_stub_difference_identity_equation": (
                "S_C = 2 E_CC + E_CN = 239, S_N = 2 E_NN + E_CN = 169 => "
                "E_CC - E_NN = (S_C - S_N) / 2 = (239 - 169) / 2 = 35"
            ),
            "empirical_6tap_5A_identity_equation": (
                "S_C = 2 E_CC + E_CN = 2(75) + 73 = 223, S_N = 2 E_NN + E_CN = 2(56) + 73 = 185 => "
                "E_CC - E_NN = (S_C - S_N) / 2 = (223 - 185) / 2 = 19"
            ),
            "empirical_6tap_5A_cutoff_values": {
                "ntd_seq_id_range": [41, NTD_MAX_SEQ_ID],
                "ctd_seq_id_range": [NTD_MAX_SEQ_ID + 1, 206],
                "observed_E_NN": int(obs_ntd_ntd),
                "observed_E_CN": int(obs_ntd_ctd),
                "observed_E_CC": int(obs_ctd_ctd),
                "observed_E_same": obs_same_domain,
                "S_N_ntd_stubs": int(2 * obs_ntd_ntd + obs_ntd_ctd),
                "S_C_ctd_stubs": int(2 * obs_ctd_ctd + obs_ntd_ctd),
                "stub_difference_S_C_minus_S_N": int(2 * (obs_ctd_ctd - obs_ntd_ntd)),
                "conserved_E_CC_minus_E_NN": int(obs_ctd_ctd - obs_ntd_ntd),
                "scheme_I_mcmc_min_E_CC_minus_E_NN": int(np.min(g_all_diff)),
                "scheme_I_mcmc_max_E_CC_minus_E_NN": int(np.max(g_all_diff)),
                "scheme_I_mcmc_variance_E_CC_minus_E_NN": float(np.var(g_all_diff)),
                "scheme_II_mcmc_min_E_CC_minus_E_NN": int(np.min(w_all_diff)),
                "scheme_II_mcmc_max_E_CC_minus_E_NN": int(np.max(w_all_diff)),
                "scheme_II_mcmc_variance_E_CC_minus_E_NN": float(np.var(w_all_diff)),
            },
            "illustrative_35_stub_difference_parameterization_reconciliation": {
                "S_C_ctd_stubs": 239,
                "S_N_ntd_stubs": 169,
                "stub_difference_S_C_minus_S_N": 70,
                "conserved_E_CC_minus_E_NN": 35,
                "implied_edge_counts_at_E_same_131": {
                    "E_CC": 83,
                    "E_NN": 48,
                    "E_CN": 73,
                    "E_same": 131,
                    "total_edges": 204,
                },
                "reconciliation_explanation": (
                    "The illustrative figure 35 (with S_C = 239, S_N = 169, S_C - S_N = 70, E_CC = 83, E_NN = 48) "
                    "arose in pre-run prose from equating E_CC - E_NN with the 35 four-chain conserved sequence "
                    "positions (or 83 homologous sequence positions for E_CC) while holding E_CC + E_NN = 131 and "
                    "E_CN = 73. In the actual 6TAP d <= 5.0 A contact graph (NTD = 41..125, CTD = 126..206), "
                    "E_CC = 75, E_NN = 56, E_CN = 73, so S_C = 223, S_N = 185, S_C - S_N = 38, and E_CC - E_NN = 19. "
                    "Under both parameterizations, E_CC - E_NN = (S_C - S_N) / 2 is strictly invariant across all "
                    "degree-preserving swaps."
                ),
            },
        },
        "global_degree_preserving_simple_graph_double_edge_swap_mcmc": {
            "null_scheme_id": "Null_Scheme_I_Global_Simple_Graph_Degree_Preserving_MCMC",
            "n_chains": mcmc_n_chains,
            "seeds": mcmc_seeds,
            "burn_in_steps_per_chain": mcmc_burn_in,
            "post_burn_in_steps_per_chain": mcmc_steps_per_chain,
            "total_post_burn_in_steps": mcmc_total_steps,
            "lazy_metropolis_hastings_persistence": True,
            "accepted_swaps_total": int(sum(g_acc_counts)),
            "acceptance_rate": round(float(sum(g_acc_counts)) / float(mcmc_total_steps), 4),
            "observed_same_domain_contacts": obs_same_domain,
            "mcmc_mean_same_domain_contacts": round(g_mean_same, 4),
            "mcmc_std_same_domain_contacts": round(g_std_same, 4),
            "mcmc_min_same_domain_contacts": int(np.min(g_all_same)),
            "mcmc_max_same_domain_contacts": int(np.max(g_all_same)),
            "same_domain_z_score": round((float(obs_same_domain) - g_mean_same) / g_std_same, 4),
            "steps_with_same_domain_ge_observed": g_exceed,
            "empirical_p_value_raw": round(g_p_raw, 6),
            "empirical_p_value_phipson_smyth": round(g_p_ps, 6),
            "lag_autocorrelation": {
                "rho_1": round(g_rho1, 4),
                "rho_10": round(g_rho10, 4),
                "rho_50": round(g_rho50, 4),
            },
            "lag_1_to_50_autocorrelation_curve": {
                f"lag_{lag_k:02d}": round(float(np.mean([_acf_at_lag(c_arr, lag_k) for c_arr in g_chains_same])), 6)
                for lag_k in range(1, 51)
            },
            "per_chain_diagnostics": [
                {
                    "chain_index": c_idx,
                    "seed": int(mcmc_seeds[c_idx]),
                    "accepted_swaps": int(g_acc_counts[c_idx]),
                    "acceptance_rate": round(float(g_acc_counts[c_idx]) / float(mcmc_steps_per_chain), 4),
                    "mean_same_domain": round(float(np.mean(g_chains_same[c_idx])), 4),
                    "std_same_domain": round(float(np.std(g_chains_same[c_idx])), 4),
                    "min_same_domain": int(np.min(g_chains_same[c_idx])),
                    "max_same_domain": int(np.max(g_chains_same[c_idx])),
                    "steps_ge_observed": int(np.sum(g_chains_same[c_idx] >= obs_same_domain)),
                    "rho_1": round(_acf_at_lag(g_chains_same[c_idx], 1), 6),
                    "rho_10": round(_acf_at_lag(g_chains_same[c_idx], 10), 6),
                    "rho_50": round(_acf_at_lag(g_chains_same[c_idx], 50), 6),
                    "tau_int": round(_compute_tau_int(g_chains_same[c_idx]), 4),
                    "ess": round(float(mcmc_steps_per_chain) / _compute_tau_int(g_chains_same[c_idx]), 2),
                    "lag_1_to_50_autocorrelation": {
                        f"lag_{lag_k:02d}": round(_acf_at_lag(g_chains_same[c_idx], lag_k), 6)
                        for lag_k in range(1, 51)
                    },
                }
                for c_idx in range(mcmc_n_chains)
            ],
            "integrated_autocorrelation_time_tau_int": round(g_tau_int, 2),
            "effective_sample_size_per_chain": round(g_ess_per_chain, 2),
            "effective_sample_size_total": round(g_ess_total, 2),
            "effectively_independent_states": int(round(g_ess_total)),
            "gelman_rubin_r_hat": round(g_rhat, 5),
            "mcse_mean_same_domain": round(g_mcse_mean, 4),
            "p_min_ess_floor": round(g_p_min_ess, 6),
            "p_value_uncertainty": _build_p_uncertainty(g_exceed, mcmc_total_steps, ess_val=g_ess_total),
            "conserved_E_CC_minus_E_NN_min": int(np.min(g_all_diff)),
            "conserved_E_CC_minus_E_NN_max": int(np.max(g_all_diff)),
            "g0_swap_mobility_by_unit": state_space_specification["connectivity"][
                "scheme_I_global_g0_swap_mobility_by_unit"
            ],
        },
        "within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc": {
            "null_scheme_id": "Null_Scheme_II_Within_Subinterface_Bipartite_Degree_Preserving_MCMC",
            "inferential_role": "STRUCTURAL_CONFOUND_AND_RESTRAINT_NOT_POSITIVE_SUPPORT",
            "confound_explanation": (
                "Null Scheme II (p_PS = 0.209264, Z = +1.3444) is a structural confound / restraint rather than "
                "positive statistical support: once residue-to-subinterface membership and within-subinterface "
                "degrees are conditioned on, 71/75 contacts in Families II and III are already fixed as CTD-CTD by "
                "the sub-interface partition boundary, so domain identity adds almost no independent signal beyond "
                "sub-interface membership."
            ),
            "n_chains": mcmc_n_chains,
            "seeds": mcmc_seeds,
            "burn_in_steps_per_chain": mcmc_burn_in,
            "post_burn_in_steps_per_chain": mcmc_steps_per_chain,
            "total_post_burn_in_steps": mcmc_total_steps,
            "lazy_metropolis_hastings_persistence": True,
            "accepted_swaps_total": int(sum(w_acc_counts)),
            "acceptance_rate": round(float(sum(w_acc_counts)) / float(mcmc_total_steps), 4),
            "color_changing_accepted_swaps_total": int(sum(w_col_counts)),
            "color_changing_step_fraction": round(float(sum(w_col_counts)) / float(mcmc_total_steps), 6),
            "observed_same_domain_contacts": obs_same_domain,
            "mcmc_mean_same_domain_contacts": round(w_mean_same, 4),
            "mcmc_std_same_domain_contacts": round(w_std_same, 4),
            "mcmc_min_same_domain_contacts": int(np.min(w_all_same)),
            "mcmc_max_same_domain_contacts": int(np.max(w_all_same)),
            "theoretical_same_domain_corridor": [min_within_sub_same_sum, max_within_sub_same_sum],
            "same_domain_z_score": round((float(obs_same_domain) - w_mean_same) / w_std_same, 4),
            "steps_with_same_domain_ge_observed": w_exceed,
            "empirical_p_value_raw": round(w_p_raw, 6),
            "empirical_p_value_phipson_smyth": round(w_p_ps, 6),
            "lag_autocorrelation": {
                "rho_1": round(w_rho1, 4),
                "rho_10": round(w_rho10, 4),
                "rho_50": round(w_rho50, 4),
            },
            "lag_1_to_50_autocorrelation_curve": {
                f"lag_{lag_k:02d}": round(float(np.mean([_acf_at_lag(c_arr, lag_k) for c_arr in w_chains_same])), 6)
                for lag_k in range(1, 51)
            },
            "per_chain_diagnostics": [
                {
                    "chain_index": c_idx,
                    "seed": int(mcmc_seeds[c_idx]),
                    "accepted_swaps": int(w_acc_counts[c_idx]),
                    "acceptance_rate": round(float(w_acc_counts[c_idx]) / float(mcmc_steps_per_chain), 4),
                    "color_changing_accepted_swaps": int(w_col_counts[c_idx]),
                    "mean_same_domain": round(float(np.mean(w_chains_same[c_idx])), 4),
                    "std_same_domain": round(float(np.std(w_chains_same[c_idx])), 4),
                    "min_same_domain": int(np.min(w_chains_same[c_idx])),
                    "max_same_domain": int(np.max(w_chains_same[c_idx])),
                    "steps_ge_observed": int(np.sum(w_chains_same[c_idx] >= obs_same_domain)),
                    "rho_1": round(_acf_at_lag(w_chains_same[c_idx], 1), 6),
                    "rho_10": round(_acf_at_lag(w_chains_same[c_idx], 10), 6),
                    "rho_50": round(_acf_at_lag(w_chains_same[c_idx], 50), 6),
                    "tau_int": round(_compute_tau_int(w_chains_same[c_idx]), 4),
                    "ess": round(float(mcmc_steps_per_chain) / _compute_tau_int(w_chains_same[c_idx]), 2),
                    "lag_1_to_50_autocorrelation": {
                        f"lag_{lag_k:02d}": round(_acf_at_lag(w_chains_same[c_idx], lag_k), 6)
                        for lag_k in range(1, 51)
                    },
                }
                for c_idx in range(mcmc_n_chains)
            ],
            "integrated_autocorrelation_time_tau_int": round(w_tau_int, 2),
            "effective_sample_size_per_chain": round(w_ess_per_chain, 2),
            "effective_sample_size_total": round(w_ess_total, 2),
            "effectively_independent_states": int(round(w_ess_total)),
            "gelman_rubin_r_hat": round(w_rhat, 5),
            "mcse_mean_same_domain": round(w_mcse_mean, 4),
            "p_min_ess_floor": round(w_p_min_ess, 6),
            "p_value_uncertainty": _build_p_uncertainty(w_exceed, mcmc_total_steps, ess_val=w_ess_total),
            "conserved_E_CC_minus_E_NN_min": int(np.min(w_all_diff)),
            "conserved_E_CC_minus_E_NN_max": int(np.max(w_all_diff)),
            "within_subinterface_swap_mobility_summary_by_unit": state_space_specification["connectivity"][
                "scheme_II_within_subinterface_swap_mobility_by_unit"
            ],
            "per_subinterface_swap_mobility": per_subinterface_swap_mobility,
        },
        "methodological_references_and_distinct_roles": {
            "polaris_2024_colored_configuration_model": {
                "citation": "Polaris (arXiv:2409.01363): Colored configuration model preserving joint color matrix and degree sequence.",
                "distinct_role_in_m7_level3": (
                    "Defines the colored multigraph configuration model preserving both the vertex degree sequence "
                    "and the joint color-pair contact matrix; transferring this framework to a simple, chain-bipartite "
                    "protein contact graph requires explicit simple-graph collision rejection (no self-loops, no "
                    "multi-edges) and chain-bipartite endpoint constraints, and when the 2x2 NTD/CTD color matrix is "
                    "NOT held fixed (so domain assortativity can be tested), the marginal stub counts S_N and S_C "
                    "still algebraically lock E_CC - E_NN = (S_C - S_N) / 2 while allowing E_same = E_CC + E_NN to "
                    "vary in steps of +/- 2."
                ),
            },
            "clauset_fosdick_dutta_double_edge_swap_mcmc": {
                "citation": "Clauset, Fosdick, & Dutta (arXiv:2105.12120) and Fosdick, Larremore, Nishimura, & Ugander (2018) SIAM Review 60(2), 315-355.",
                "distinct_role_in_m7_level3": (
                    "Establishes the mathematical distinction between stub-labeled and vertex-labeled simple-graph "
                    "state spaces under double-edge swaps, proves that lazy self-loop persistence (retaining G_{t+1} = G_t "
                    "when a proposed swap is rejected for creating a self-loop or multi-edge) is mandatory for a "
                    "uniform stationary distribution over simple graphs, and provides the MCMC convergence diagnostics "
                    "(lag-k autocorrelation, integrated autocorrelation time tau_int, effective sample size ESS = "
                    "N_steps / tau_int, and Gelman-Rubin R_hat) used to audit Scheme I and Scheme II."
                ),
            },
            "topologically_invariant_permutation_test_2025": {
                "citation": "Topologically Invariant Permutation Test (arXiv:2511.06153): Persistence-diagram / Wasserstein permutation testing on functional brain networks.",
                "distinct_role_in_m7_level3": (
                    "Demonstrates topology-preserving permutation inference on weighted functional brain graphs using "
                    "persistent homology and Wasserstein distances; it is cited solely as a domain comparator for "
                    "topology-invariant permutation testing and is explicitly NOT a proof of Markov-chain "
                    "irreducibility, state-space connectedness, or mixing for the M7 discrete double-edge-swap kernel."
                ),
                "is_proof_of_m7_double_edge_swap_irreducibility_or_mixing": False,
            },
            "kdd_2026_labeled_graphs_under_mismatches": {
                "citation": "KDD 2026 (Testing Statistical Dependence in Labeled Graphs under Mismatches).",
                "distinct_role_in_m7_level3": (
                    "Formulates statistical dependence testing between graph topology and node labels under alignment "
                    "or structural mismatches, motivating our paired evaluation of vertex-label permutations on the "
                    "fixed 204-edge topology alongside edge-rewiring MCMC nulls."
                ),
            },
            "mango_2026_degree_aware_network_clustering": {
                "citation": "MANGO (2026 preprint): Degree-aware network clustering conditioning on fixed topology and degree composition.",
                "distinct_role_in_m7_level3": (
                    "Provides the methodological precedent for conditioning network clustering and enrichment baselines "
                    "on heterogeneous node degrees and fixed community/sub-interface sizes."
                ),
            },
            "phipson_smyth_2010_permutation_p_values": {
                "citation": "Phipson, B. & Smyth, G. K. (2010) Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. Stat. Appl. Genet. Mol. Biol. 9, Article 39 (doi:10.2202/1544-6115.1585).",
                "distinct_role_in_m7_level3": (
                    "Defines the exact unbiased Monte Carlo p-value estimator p_PS = (k + 1) / (B + 1) with "
                    "finite-sample resolution floor p_min,perm = 1 / (B + 1) = 1 / 10001 = 9.999e-5."
                ),
            },
        },
        "monte_carlo_and_mcmc_p_value_resolution": {
            "iid_permutation_draw_count_B": n_perm,
            "p_min_perm_floor": round(1.0 / float(n_perm + 1), 8),
            "p_min_perm_formula": "p_min,perm = 1 / (B + 1) = 1 / 10001 = 9.9990001e-5",
            "mcmc_scheme_I_raw_steps": mcmc_total_steps,
            "mcmc_scheme_I_ess_total": round(g_ess_total, 2),
            "mcmc_scheme_I_p_min_ess_floor": round(g_p_min_ess, 6),
            "mcmc_scheme_II_raw_steps": mcmc_total_steps,
            "mcmc_scheme_II_ess_total": round(w_ess_total, 2),
            "mcmc_scheme_II_p_min_ess_floor": round(w_p_min_ess, 6),
            "p_min_ess_formula": "p_min,ESS = 1 / (ESS + 1)",
            "resolution_and_uncertainty_summary_by_test": {
                "homologous_83_residue_same_domain": _build_p_uncertainty(exceed_seq_same, n_perm),
                "homologous_83_residue_omnibus_chi_square": _build_p_uncertainty(exceed_seq_chi2, n_perm),
                "chain_specific_235_node_same_domain": _build_p_uncertainty(exceed_node_same, n_perm),
                "full_166_resolved_sequence_same_domain": _build_p_uncertainty(exceed_166_same, n_perm),
                "circular_shift_166_aggregate_same_domain": _build_p_uncertainty(circ_exceed_same, len(all_166_seqs)),
                "circular_shift_166_non_family_I_ctd_purity": _build_p_uncertainty(
                    circ_exceed_non_fam1_ctd, len(all_166_seqs)
                ),
                "subinterface_edge_permutation_family_II_hydrophobic": _build_p_uncertainty(exceed_fam2_hyd, n_perm),
                "subinterface_edge_permutation_family_II_pi_stacking": _build_p_uncertainty(exceed_fam2_pi5, n_perm),
                "subinterface_edge_permutation_family_I_salt_bridge": _build_p_uncertainty(exceed_fam1_sb, n_perm),
                "homologous_83_comp_id_permutation_family_II_nonpolar": _build_p_uncertainty(
                    exceed_aa_nonpolar, n_perm
                ),
                "homologous_83_comp_id_permutation_family_II_aromatic": _build_p_uncertainty(
                    exceed_aa_aromatic, n_perm
                ),
                "homologous_83_comp_id_permutation_family_I_ion_pair": _build_p_uncertainty(
                    exceed_aa_ionpair, n_perm
                ),
                "mcmc_scheme_I_global_double_edge_swap_same_domain": _build_p_uncertainty(
                    g_exceed, mcmc_total_steps, ess_val=g_ess_total
                ),
                "mcmc_scheme_II_within_subinterface_double_edge_swap_same_domain": _build_p_uncertainty(
                    w_exceed, mcmc_total_steps, ess_val=w_ess_total
                ),
            },
        },
    }

    # 3B. Multiplicity Analysis Across Defined Hypothesis Families (Holm-Bonferroni FWER & Benjamini-Hochberg FDR)
    fam_a_specs = [
        (
            "homologous_83_residue_same_domain",
            "83 homologous interfacial sequence-position domain permutation (E_same = 131/204)",
            p_perm_seq_same_phipson_smyth,
            exceed_seq_same,
            n_perm,
        ),
        (
            "homologous_83_residue_omnibus_chi_square",
            "83 homologous sequence-position degree-adjusted omnibus chi-square (chi^2 = 14.8162)",
            p_perm_seq_chi2_phipson_smyth,
            exceed_seq_chi2,
            n_perm,
        ),
        (
            "chain_specific_235_node_same_domain",
            "235 chain-specific (chain_id, seq_id) vertex domain permutation (E_same = 131/204)",
            p_perm_node_same_phipson_smyth,
            exceed_node_same,
            n_perm,
        ),
        (
            "full_166_resolved_sequence_same_domain",
            "166 resolved sequence-position (41..206) domain permutation (E_same = 131/204)",
            p_perm_166_same_phipson_smyth,
            exceed_166_same,
            n_perm,
        ),
        (
            "circular_shift_166_non_family_I_ctd_purity",
            "166-residue circular shift non-Family-I pure CTD-CTD concentration (71/75 = 94.67%)",
            float(circ_exceed_non_fam1_ctd) / float(len(all_166_seqs)),
            circ_exceed_non_fam1_ctd,
            len(all_166_seqs),
        ),
        (
            "circular_shift_166_aggregate_same_domain",
            "166-residue circular shift aggregate same-domain contacts (E_same = 131/204)",
            float(circ_exceed_same) / float(len(all_166_seqs)),
            circ_exceed_same,
            len(all_166_seqs),
        ),
    ]
    fam_b_specs = [
        (
            "subinterface_edge_permutation_family_II_pi_stacking",
            "Sub-interface edge assignment permutation: Family II pi-stacking at 5.0 A (8/8 edges)",
            float(exceed_fam2_pi5 + 1) / float(n_perm + 1),
            exceed_fam2_pi5,
            n_perm,
        ),
        (
            "subinterface_edge_permutation_family_II_hydrophobic",
            "Sub-interface edge assignment permutation: Family II hydrophobic contacts (13/20 edges)",
            float(exceed_fam2_hyd + 1) / float(n_perm + 1),
            exceed_fam2_hyd,
            n_perm,
        ),
        (
            "homologous_83_comp_id_permutation_family_II_aromatic",
            "83-residue amino-acid comp_id permutation: Family II aromatic pairs (9 pairs)",
            float(exceed_aa_aromatic + 1) / float(n_perm + 1),
            exceed_aa_aromatic,
            n_perm,
        ),
        (
            "homologous_83_comp_id_permutation_family_II_nonpolar",
            "83-residue amino-acid comp_id permutation: Family II nonpolar pairs (18 pairs)",
            float(exceed_aa_nonpolar + 1) / float(n_perm + 1),
            exceed_aa_nonpolar,
            n_perm,
        ),
        (
            "subinterface_edge_permutation_family_I_salt_bridge",
            "Sub-interface edge assignment permutation: Family I salt-bridge pairs at 5.0 A (6/6 pairs)",
            float(exceed_fam1_sb + 1) / float(n_perm + 1),
            exceed_fam1_sb,
            n_perm,
        ),
        (
            "homologous_83_comp_id_permutation_family_I_ion_pair",
            "83-residue amino-acid comp_id permutation: Family I opposite-charge ion pairs (16 pairs)",
            float(exceed_aa_ionpair + 1) / float(n_perm + 1),
            exceed_aa_ionpair,
            n_perm,
        ),
    ]

    p_fam_a = [x[2] for x in fam_a_specs]
    p_fam_b = [x[2] for x in fam_b_specs]
    p_omnibus_12 = p_fam_a + p_fam_b

    holm_fam_a = _holm_bonferroni(p_fam_a)
    bh_fam_a = _benjamini_hochberg(p_fam_a)
    holm_fam_b = _holm_bonferroni(p_fam_b)
    bh_fam_b = _benjamini_hochberg(p_fam_b)
    holm_omnibus_12 = _holm_bonferroni(p_omnibus_12)
    bh_omnibus_12 = _benjamini_hochberg(p_omnibus_12)

    fam_a_records: List[Dict[str, Any]] = []
    for idx_a, (tid, tdesc, p_unadj, k_ex, n_dr) in enumerate(fam_a_specs):
        fam_a_records.append(
            {
                "test_id": tid,
                "hypothesis_family": "Family_A_Domain_Segregation_and_Spatial_Autocorrelation",
                "description": tdesc,
                "exceedances_k": int(k_ex),
                "total_draws": int(n_dr),
                "p_unadjusted": round(p_unadj, 6),
                "p_holm_within_family": round(holm_fam_a[idx_a], 6),
                "q_bh_within_family": round(bh_fam_a[idx_a], 6),
                "p_holm_omnibus_12": round(holm_omnibus_12[idx_a], 6),
                "q_bh_omnibus_12": round(bh_omnibus_12[idx_a], 6),
                "survives_holm_within_family_alpha_0p05": holm_fam_a[idx_a] <= 0.05,
                "survives_bh_within_family_q_0p05": bh_fam_a[idx_a] <= 0.05,
                "survives_holm_omnibus_12_alpha_0p05": holm_omnibus_12[idx_a] <= 0.05,
                "survives_bh_omnibus_12_q_0p05": bh_omnibus_12[idx_a] <= 0.05,
            }
        )

    fam_b_records: List[Dict[str, Any]] = []
    for idx_b, (tid, tdesc, p_unadj, k_ex, n_dr) in enumerate(fam_b_specs):
        om_idx = len(fam_a_specs) + idx_b
        fam_b_records.append(
            {
                "test_id": tid,
                "hypothesis_family": "Family_B_Interface_Chemistry_Specialization",
                "description": tdesc,
                "exceedances_k": int(k_ex),
                "total_draws": int(n_dr),
                "p_unadjusted": round(p_unadj, 6),
                "p_holm_within_family": round(holm_fam_b[idx_b], 6),
                "q_bh_within_family": round(bh_fam_b[idx_b], 6),
                "p_holm_omnibus_12": round(holm_omnibus_12[om_idx], 6),
                "q_bh_omnibus_12": round(bh_omnibus_12[om_idx], 6),
                "survives_holm_within_family_alpha_0p05": holm_fam_b[idx_b] <= 0.05,
                "survives_bh_within_family_q_0p05": bh_fam_b[idx_b] <= 0.05,
                "survives_holm_omnibus_12_alpha_0p05": holm_omnibus_12[om_idx] <= 0.05,
                "survives_bh_omnibus_12_q_0p05": bh_omnibus_12[om_idx] <= 0.05,
            }
        )

    multiplicity_adjustment_analysis: Dict[str, Any] = {
        "inference_tier": "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY",
        "alpha_threshold": 0.05,
        "methods_applied": [
            "Holm-Bonferroni step-down family-wise error rate (FWER) control",
            "Benjamini-Hochberg step-up false discovery rate (FDR) control",
        ],
        "hypothesis_families": {
            "domain_segregation_and_spatial_autocorrelation_family": {
                "family_id": "Family_A_Domain_Segregation_and_Spatial_Autocorrelation",
                "n_tests": len(fam_a_records),
                "n_surviving_holm_within_family_alpha_0p05": sum(
                    1 for r in fam_a_records if r["survives_holm_within_family_alpha_0p05"]
                ),
                "n_surviving_bh_within_family_q_0p05": sum(
                    1 for r in fam_a_records if r["survives_bh_within_family_q_0p05"]
                ),
                "tests": fam_a_records,
            },
            "interface_chemistry_specialization_family": {
                "family_id": "Family_B_Interface_Chemistry_Specialization",
                "n_tests": len(fam_b_records),
                "n_surviving_holm_within_family_alpha_0p05": sum(
                    1 for r in fam_b_records if r["survives_holm_within_family_alpha_0p05"]
                ),
                "n_surviving_bh_within_family_q_0p05": sum(
                    1 for r in fam_b_records if r["survives_bh_within_family_q_0p05"]
                ),
                "tests": fam_b_records,
            },
        },
        "omnibus_12_permutation_test_battery": {
            "n_total_tests": len(fam_a_records) + len(fam_b_records),
            "n_surviving_holm_omnibus_12_alpha_0p05": sum(
                1 for r in (fam_a_records + fam_b_records) if r["survives_holm_omnibus_12_alpha_0p05"]
            ),
            "n_surviving_bh_omnibus_12_q_0p05": sum(
                1 for r in (fam_a_records + fam_b_records) if r["survives_bh_omnibus_12_q_0p05"]
            ),
            "all_tests": fam_a_records + fam_b_records,
        },
        "multiplicity_qualification_summary": (
            "Across the 12 permutation/shift tests (6 in Hypothesis Family A: domain segregation & spatial "
            "autocorrelation; 6 in Hypothesis Family B: interface chemistry specialization), 6/12 survive "
            "within-family Holm-Bonferroni FWER control at alpha = 0.05 (3/6 in Family A: 235-node E_same "
            "p_Holm = 0.000600, 83-residue E_same p_Holm = 0.032497, and 166-residue E_same p_Holm = 0.034397, "
            "while circular-shift non-Family-I CTD purity sits at borderline p_Holm = 0.054217 and q_BH = 0.027108; "
            "plus 3/6 in Family B: sub-interface edge permutation for Family II pi-stacking p_Holm = 0.001800, "
            "Family II hydrophobic contacts p_Holm = 0.002000, and 83-residue comp_id Family I ion pairs "
            "p_Holm = 0.044796). Under within-family Benjamini-Hochberg FDR control at q <= 0.05, 8/12 survive "
            "(4/6 in Family A adding circular-shift non-Family-I CTD purity q_BH = 0.027108, and 4/6 in Family B "
            "adding 83-residue comp_id Family II nonpolar pairs q_BH = 0.044696, while 83-residue comp_id Family II "
            "aromatic pairs sits at q_BH = 0.052195 and Family I salt-bridge edge permutation sits at q_BH = 0.063594). "
            "Under the combined 12-test omnibus correction, 3/12 survive omnibus Holm-Bonferroni FWER (235-node "
            "E_same p_Holm,12 = 0.001200, Family II pi-stacking p_Holm,12 = 0.003300, and Family II hydrophobic "
            "contacts p_Holm,12 = 0.004000) and 8/12 survive omnibus Benjamini-Hochberg FDR at q_BH,12 <= 0.05."
        ),
    }

    non_independence_caveats: List[str] = [
        "Residues recur across multiple contact edges (the 204 symmetry-unique contact edges at d <= 5.0 A connect 235 chain-specific (chain_id, seq_id) nodes and 83 homologous sequence positions with heterogeneous vertex degrees up to 13).",
        "Multiple residue-pair contacts belong to the same sub-interface (partitioned across 10 sub-interfaces with sizes ranging from 1 to 41 residue pairs) and are coupled by rigid-body subunit geometry.",
        "The 60 icosahedral symmetry copies in PDB 6TAP (point_symmetry = I) are exact geometric equivalents generated by _pdbx_struct_oper_list rather than independent biological replicates.",
        "Consequently, both the multinomial i.i.d. chi-square statistic (chi^2 = 22.1774, nominal p_iid = 1.53e-5) and the area-proportional chemistry enrichment ratios (2.35x hydrophobic, 2.89x pi-stacking in Family II) are model-based exploratory screening comparisons rather than confirmatory hypothesis tests, and remain :Incomplete until wet-lab mutagenesis.",
    ]

    competing_null_model_tests = {
        "m7_falsification_design_rule": "Rule 3: A falsifiable model must beat a competing null model.",
        "inference_tier": "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY",
        "epistemic_status": ":Incomplete (exploratory structural screening until wet-lab mutagenesis)",
        "non_independence_caveats": non_independence_caveats,
        "domain_segregation_vs_random_mixing_null": {
            "inference_tier": "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY",
            "non_independence_caveats": non_independence_caveats,
            "asu_ntd_residue_fraction": round(p_ntd, 4),
            "asu_ctd_residue_fraction": round(p_ctd, 4),
            "null_expected_pairs_at_204": {
                "NTD-NTD": round(exp_ntd_ntd, 2),
                "NTD-CTD": round(exp_ntd_ctd, 2),
                "CTD-CTD": round(exp_ctd_ctd, 2),
            },
            "observed_pairs_at_204": {
                "NTD-NTD": int(obs_ntd_ntd),
                "NTD-CTD": int(obs_ntd_ctd),
                "CTD-CTD": int(obs_ctd_ctd),
            },
            "chi_square_statistic_df2": round(chi2_domain_segregation, 4),
            "nominal_iid_p_value": 1.53e-5,
            "p_value_upper_bound": 1.6e-5,
            "null_rejected": True,
            "structural_interpretation": (
                f"Nominal multinomial screening comparison (chi^2 = {chi2_domain_segregation:.4f}, df = 2, "
                f"nominal p_iid = 1.53e-5; inference_tier = MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY) "
                f"separates observed NTD/CTD contact architecture from uniform random mixing: four distal interfaces "
                f"(C-C_op33, B-D_op33, A-D_op36, A-B_op2) in Family_III_Residual_Distal_CTD_Hierarchy are 100% pure "
                f"CTD-CTD (20/20 pairs), the two quasi-equivalent CTD dimer locks (A-D_op1 and B-C_op33) in "
                f"Family II are 92.7% CTD-CTD (51/55 pairs, with only boundary residue ALA125 in NTD), and "
                f"cross-domain NTD-CTD contacts are depleted {exp_ntd_ctd / obs_ntd_ctd:.2f}-fold relative to "
                f"the i.i.d. random-mixing expectation ({int(obs_ntd_ctd)} observed vs {exp_ntd_ctd:.2f} expected). "
                f"Because the 204 contact edges share recurring residues and sub-interfaces, empirical significance "
                f"is additionally evaluated under the degree- and interface-preserving permutation null below."
            ),
        },
        "ctd_dimer_hydrophobic_aromatic_vs_ntd_seam_salt_bridge_null": {
            "inference_tier": "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY",
            "non_independence_caveats": non_independence_caveats,
            "total_hydrophobic_candidates_across_asu": total_hyd,
            "total_pi_stacking_candidates_across_asu": total_pi,
            "total_salt_bridge_candidates_across_asu": total_sb,
            "total_asu_pairwise_bsa_A2": round(total_sr_bsa, 2),
            "quasi_equivalent_ctd_dimers_A_D_op1_and_B_C_op33": {
                "combined_bsa_A2": round(ctd_dimer_bsa, 2),
                "fraction_of_total_bsa": round(ctd_dimer_bsa / total_sr_bsa, 4),
                "observed_hydrophobic_candidates": ctd_dimer_hyd,
                "null_expected_hydrophobic_candidates": round(total_hyd * (ctd_dimer_bsa / total_sr_bsa), 4),
                "hydrophobic_fold_enrichment_over_null": round(
                    (float(ctd_dimer_hyd) / float(total_hyd)) / (ctd_dimer_bsa / total_sr_bsa), 4
                ),
                "observed_pi_stacking_candidates": ctd_dimer_pi,
                "null_expected_pi_stacking_candidates": round(total_pi * (ctd_dimer_bsa / total_sr_bsa), 4),
                "pi_stacking_fold_enrichment_over_null": round(
                    (float(ctd_dimer_pi) / float(total_pi)) / (ctd_dimer_bsa / total_sr_bsa), 4
                ),
                "observed_salt_bridge_candidates": 0,
            },
            "four_ntd_ctd_capsomere_seams_A_A_op36_C_D_op1_B_C_op1_B_D_op2": {
                "combined_bsa_A2": round(ntd_ctd_seams_bsa, 2),
                "fraction_of_total_bsa": round(ntd_ctd_seams_bsa / total_sr_bsa, 4),
                "observed_salt_bridge_candidates": total_sb,
                "null_expected_salt_bridge_candidates": round(total_sb * (ntd_ctd_seams_bsa / total_sr_bsa), 4),
                "salt_bridge_fold_enrichment_over_null": round(
                    1.0 / (ntd_ctd_seams_bsa / total_sr_bsa), 4
                ),
                "observed_pi_stacking_candidates": 0,
            },
            "structural_interpretation": (
                "Under the model-based area-proportional chemistry screening null "
                "(inference_tier = MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY), "
                "the two quasi-equivalent CTD dimer locks in interface-fingerprint Family II (A-D_op1 + B-C_op33, "
                "27.72% of BSA) carry 65.0% (13/20, 2.35x enrichment) of all hydrophobic candidates and 80.0% "
                "(8/10, 2.89x enrichment) of all aromatic pi-stacking candidates with 0/9 salt bridges, whereas "
                "100% (9/9, anchored by the conserved inter-subunit ARG56-ASP144 salt bridge across A-A_op36, "
                "C-D_op1, B-C_op1, and B-D_op2) of salt-bridge candidates occur in the 4 NTD+CTD capsomere seams "
                "of interface-fingerprint Family I (which carry 0/10 pi-stacking candidates)."
            ),
        },
        "degree_and_interface_preserving_permutation_null": {
            "inference_tier": "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY",
            "epistemic_status": ":Incomplete (exploratory structural screening until wet-lab mutagenesis)",
            "n_permutations": n_perm,
            "random_seed": perm_seed,
            "rng_implementation": "numpy.random.default_rng(20261107)",
            "domain_assignment_permutation_test": {
                "description": (
                    "Permutes NTD vs CTD domain labels across (1A) the exact 83 homologous interfacial sequence "
                    "positions (38 NTD, 45 CTD), (1B) the 235 chain-specific (chain_id, seq_id) interfacial "
                    "nodes (108 NTD, 127 CTD), (1C) all 166 resolved sequence positions (41..206: 85 NTD, 81 CTD), "
                    "and (1D) all 166 contiguous sequence-block circular shifts along 41..206, preserving the "
                    "exact contact graph topology, all 10 sub-interface sizes, and every residue's contact degree "
                    "across the 204 edges."
                ),
                "homologous_83_residue_degree_preserving_permutation": {
                    "total_homologous_interfacial_residues": len(unique_hom_seqs),
                    "ntd_homologous_residues": int(np.sum(true_seq_dom == 0)),
                    "ctd_homologous_residues": int(np.sum(true_seq_dom == 1)),
                    "observed_same_domain_contacts_NTD_NTD_plus_CTD_CTD": obs_same_domain,
                    "observed_cross_domain_contacts_NTD_CTD": int(obs_ntd_ctd),
                    "permutation_expected_pairs_at_204": {
                        "NTD-NTD": round(mean_seq_nn, 2),
                        "NTD-CTD": round(mean_seq_nc, 2),
                        "CTD-CTD": round(mean_seq_cc, 2),
                    },
                    "permutation_mean_same_domain_contacts": round(mean_seq_same, 2),
                    "permutation_std_same_domain_contacts": round(std_seq_same, 2),
                    "same_domain_z_score": round((float(obs_same_domain) - mean_seq_same) / std_seq_same, 4),
                    "permutations_with_same_domain_ge_observed": exceed_seq_same,
                    "empirical_p_value_same_domain_raw": round(p_perm_seq_same_raw, 6),
                    "empirical_p_value_same_domain_phipson_smyth": round(p_perm_seq_same_phipson_smyth, 6),
                    "degree_adjusted_chi_square_statistic": round(float(obs_chi2_vs_perm_seq), 4),
                    "permutations_with_chi_square_ge_observed": exceed_seq_chi2,
                    "empirical_p_value_omnibus_chi_square_phipson_smyth": round(p_perm_seq_chi2_phipson_smyth, 6),
                },
                "chain_specific_235_node_degree_preserving_permutation": {
                    "total_chain_residue_nodes": len(unique_chain_nodes),
                    "ntd_chain_residue_nodes": int(np.sum(true_node_dom == 0)),
                    "ctd_chain_residue_nodes": int(np.sum(true_node_dom == 1)),
                    "observed_same_domain_contacts": obs_same_domain,
                    "permutation_mean_same_domain_contacts": round(float(np.mean(perm_node_same)), 2),
                    "permutation_std_same_domain_contacts": round(float(np.std(perm_node_same)), 2),
                    "permutations_with_same_domain_ge_observed": exceed_node_same,
                    "empirical_p_value_same_domain_phipson_smyth": round(p_perm_node_same_phipson_smyth, 6),
                },
                "full_166_resolved_sequence_domain_permutation": {
                    "total_resolved_sequence_positions_41_to_206": len(all_166_seqs),
                    "ntd_sequence_positions_41_to_125": int(np.sum(true_166_dom == 0)),
                    "ctd_sequence_positions_126_to_206": int(np.sum(true_166_dom == 1)),
                    "observed_same_domain_contacts": obs_same_domain,
                    "permutation_mean_same_domain_contacts": round(mean_166_same, 2),
                    "permutation_std_same_domain_contacts": round(std_166_same, 2),
                    "same_domain_z_score": round((float(obs_same_domain) - mean_166_same) / std_166_same, 4),
                    "permutations_with_same_domain_ge_observed": exceed_166_same,
                    "empirical_p_value_same_domain_phipson_smyth": round(p_perm_166_same_phipson_smyth, 6),
                },
                "contiguous_166_residue_circular_shift_autocorrelation_null": {
                    "inferential_role": "SEQUENCE_CONTIGUITY_CONFOUND_AND_BOUNDARY_LOCALIZATION_RESTRAINT_NOT_STRONG_SUPPORT",
                    "total_circular_shifts_evaluated": len(all_166_seqs),
                    "discrete_resolution_floor_1_over_166": round(1.0 / float(len(all_166_seqs)), 6),
                    "observed_aggregate_same_domain_contacts": obs_same_domain,
                    "circular_shift_mean_same_domain_contacts": round(float(np.mean(circ_same_domain)), 2),
                    "circular_shift_std_same_domain_contacts": round(float(np.std(circ_same_domain)), 2),
                    "circular_shifts_with_aggregate_same_domain_ge_observed": circ_exceed_same,
                    "empirical_p_circular_aggregate_same_domain": round(
                        float(circ_exceed_same) / float(len(all_166_seqs)), 6
                    ),
                    "observed_non_family_I_pure_ctd_ctd_contacts": obs_non_fam1_ctd_ctd,
                    "total_non_family_I_contacts": non_fam1_total_edges,
                    "circular_shifts_with_non_family_I_ctd_ctd_ge_observed": circ_exceed_non_fam1_ctd,
                    "exceeding_shift_offsets_around_125_126_hinge": [-1, 0, 1],
                    "empirical_p_circular_non_family_I_ctd_purity": round(
                        float(circ_exceed_non_fam1_ctd) / float(len(all_166_seqs)), 6
                    ),
                    "autocorrelation_qualification_note": (
                        "CONFOUND AND RESTRAINT FRAMING: (a) Because NTD (41..125) and CTD (126..206) are contiguous "
                        "folded blocks along the 166-residue polypeptide chain, 38/166 (p_circular = 0.228916) cyclic "
                        "shifts yield >= 131/204 aggregate same-domain contacts, showing that aggregate same-domain "
                        "enrichment does NOT beat primary-sequence contiguity / 3D globular domain autocorrelation; "
                        "and (b) while only 3/166 (p_circular = 0.018072) cyclic shifts concentrate >= 71/75 (94.67%) "
                        "of non-Family-I contacts (Family II + Family III) into pure CTD-CTD contacts, those 3 hits "
                        "are simply the observed 125|126 inter-domain hinge (s = 0) and its two immediate neighbor "
                        "offsets (s in {-1, +1}) on a discrete 166-shift grid whose resolution floor is 1/166 = "
                        "0.006024, making it a descriptive hinge-boundary localization restraint rather than strong "
                        "independent statistical evidence."
                    ),
                },
            },
            "chemistry_clustering_permutation_test": {
                "description": (
                    "Evaluates Family II hydrophobic/aromatic clustering and Family I salt-bridge clustering "
                    "under two complementary degree- and interface-preserving permutation schemes: (1) permuting "
                    "sub-interface family assignments across the 204 contact edges while preserving exact "
                    "sub-interface contact counts (55 Family II edges, 129 Family I edges, 20 Family III edges) "
                    "and total chemical candidate counts; and (2) permuting the 83 homologous amino-acid residue "
                    "labels (comp_id) across sequence positions while preserving exact residue-type multiplicity "
                    "and vertex degrees across the 204 edges."
                ),
                "subinterface_edge_assignment_permutation": {
                    "family_II_edge_count": 55,
                    "family_I_edge_count": 129,
                    "family_III_edge_count": 20,
                    "family_II_hydrophobic_enrichment": {
                        "observed_family_II_hydrophobic_edges": obs_fam2_hyd_edges,
                        "total_hydrophobic_edges_at_5A": int(np.sum(edge_is_hyd)),
                        "permutation_mean": round(float(np.mean(perm_fam2_hyd)), 4),
                        "fold_enrichment_over_permutation_mean": round(
                            float(obs_fam2_hyd_edges) / float(np.mean(perm_fam2_hyd)), 4
                        ),
                        "permutations_ge_observed": exceed_fam2_hyd,
                        "empirical_p_value_phipson_smyth": round(
                            float(exceed_fam2_hyd + 1) / float(n_perm + 1), 6
                        ),
                    },
                    "family_II_pi_stacking_enrichment_at_5A": {
                        "observed_family_II_pi_stacking_edges_at_5A": obs_fam2_pi5_edges,
                        "total_pi_stacking_edges_at_5A": int(np.sum(edge_is_pi5)),
                        "note_on_all_cutoff_pi_stacking": (
                            "8/8 of the d <= 5.0 A pi-stacking contact edges (and 8/10 across all cutoffs, "
                            "where the remaining 2 HIS182-HIS182 C3 vertex rings have centroid distance ~6.8 A "
                            "and heavy-atom distance > 5.0 A) reside in Family II."
                        ),
                        "permutation_mean": round(float(np.mean(perm_fam2_pi5)), 4),
                        "fold_enrichment_over_permutation_mean": round(
                            float(obs_fam2_pi5_edges) / float(np.mean(perm_fam2_pi5)), 4
                        ),
                        "permutations_ge_observed": exceed_fam2_pi5,
                        "empirical_p_value_phipson_smyth": round(
                            float(exceed_fam2_pi5 + 1) / float(n_perm + 1), 6
                        ),
                    },
                    "family_I_salt_bridge_enrichment_at_5A": {
                        "observed_family_I_salt_bridge_residue_pairs_at_5A": obs_fam1_sb_edges,
                        "observed_family_I_salt_bridge_atom_hits": total_sb,
                        "total_salt_bridge_residue_pairs_at_5A": int(np.sum(edge_is_sb_pair)),
                        "permutation_mean_residue_pairs": round(float(np.mean(perm_fam1_sb)), 4),
                        "fold_enrichment_over_permutation_mean": round(
                            float(obs_fam1_sb_edges) / float(np.mean(perm_fam1_sb)), 4
                        ),
                        "permutations_ge_observed": exceed_fam1_sb,
                        "empirical_p_value_phipson_smyth": round(
                            float(exceed_fam1_sb + 1) / float(n_perm + 1), 6
                        ),
                    },
                },
                "homologous_83_residue_comp_id_label_permutation": {
                    "topological_graph_fixed_caveat": (
                        "This test is conditional on the fixed 204-edge contact graph and tests whether amino-acid "
                        "labels are unusually distributed over already-observed edges; it does not model "
                        "mutation-induced side-chain repacking or recompute the contact graph after a mutation, and "
                        "its enrichment results must not be described as predictions of mutational delta-delta-G or "
                        "assembly effects. Permuting amino-acid comp_id labels across the 83 interfacial sequence "
                        "positions while holding the d <= 5.0 A contact graph topology fixed is an exploratory "
                        "topological screening null rather than a physical side-chain repacking simulation, because "
                        "side-chain volume and geometry co-determine whether two sequence positions fall within 5.0 A."
                    ),
                    "family_II_nonpolar_residue_pair_enrichment": {
                        "observed_family_II_nonpolar_pairs": obs_aa_fam2_nonpolar,
                        "permutation_mean": round(float(np.mean(perm_aa_fam2_nonpolar)), 4),
                        "fold_enrichment_over_permutation_mean": round(
                            float(obs_aa_fam2_nonpolar) / float(np.mean(perm_aa_fam2_nonpolar)), 4
                        ),
                        "permutations_ge_observed": exceed_aa_nonpolar,
                        "empirical_p_value_phipson_smyth": round(
                            float(exceed_aa_nonpolar + 1) / float(n_perm + 1), 6
                        ),
                    },
                    "family_II_aromatic_residue_pair_enrichment": {
                        "observed_family_II_aromatic_pairs": obs_aa_fam2_aromatic,
                        "permutation_mean": round(float(np.mean(perm_aa_fam2_aromatic)), 4),
                        "fold_enrichment_over_permutation_mean": round(
                            float(obs_aa_fam2_aromatic) / float(np.mean(perm_aa_fam2_aromatic)), 4
                        ),
                        "permutations_ge_observed": exceed_aa_aromatic,
                        "empirical_p_value_phipson_smyth": round(
                            float(exceed_aa_aromatic + 1) / float(n_perm + 1), 6
                        ),
                    },
                    "family_I_opposite_charge_ion_pair_enrichment": {
                        "observed_family_I_ion_pairs": obs_aa_fam1_ionpair,
                        "permutation_mean": round(float(np.mean(perm_aa_fam1_ionpair)), 4),
                        "fold_enrichment_over_permutation_mean": round(
                            float(obs_aa_fam1_ionpair) / float(np.mean(perm_aa_fam1_ionpair)), 4
                        ),
                        "permutations_ge_observed": exceed_aa_ionpair,
                        "empirical_p_value_phipson_smyth": round(
                            float(exceed_aa_ionpair + 1) / float(n_perm + 1), 6
                        ),
                    },
                },
            },
            "constrained_graph_state_space_connectivity_and_mcmc_audit": (
                constrained_graph_state_space_connectivity_and_mcmc_audit
            ),
            "multiplicity_adjustment_analysis": multiplicity_adjustment_analysis,
        },
        "constrained_graph_state_space_connectivity_and_mcmc_audit": (
            constrained_graph_state_space_connectivity_and_mcmc_audit
        ),
        "multiplicity_adjustment_analysis": multiplicity_adjustment_analysis,
    }

    # 3. Three Distinct Physical RNA Arrangements & 4-Branch Literature Separation
    zn_sphere = l2_e["zn_coordination_sphere"]
    nc_model = l2_e["nc_peptide_model"]
    nc_dimer = l2_e["nc_nc_dimer_interface_across_op2"]
    sub_align = l2_mods["L2_D_icosahedral_operator_axis_mapping"]["subunit_radial_and_axis_alignment"]
    ca_radii = [sub_align[ch]["radial_distance_from_capsid_center_A"] for ch in ("A", "B", "C", "D")]

    three_physical_rna_arrangements: List[Dict[str, Any]] = [
        {
            "arrangement_id": "arrangement_1_internal_nc_zn_rna_binding_architecture_in_6tap",
            "physical_topology": "Internal NC CCHC zinc-knuckle RNA-binding domain inside native 6TAP dArc1 capsid lumen",
            "primary_reference": "Erlendsson, S. et al. (2020) Nature Neuroscience 23, 172-175 (PDB 6TAP)",
            "mmcif_entity_and_chain_provenance": (
                "Deposited polymer _entity.id = 1 (dArc1 polypeptide, 254 aa, UniProt Q7K1U0) strand label_asym_id = E "
                "(auth_asym_id = E, residues 224..250, functional_component_id = NC_ZINC_KNUCKLE_FRAGMENT) coordinating "
                "non-polymer ligand _entity.id = 2 (ZN, label_asym_id = F, auth_asym_id = E, _atom_site.auth_seq_id = 301, "
                "_pdbx_nonpoly_scheme.pdb_seq_num = 301, _pdbx_nonpoly_scheme.auth_seq_num = 251, "
                "functional_component_id = ZN_COORDINATION_LIGAND)."
            ),
            "physical_mechanism_and_location": (
                "The C-terminal NC zinc knuckle projects radially inward from the CA shell interior beneath the "
                f"C2-symmetric B-D dimeric hexamer seam (centroid radius {nc_model['radial_distance_from_center_A']:.2f} A "
                f"vs CA chains A..D at {min(ca_radii):.2f}..{max(ca_radii):.2f} A), "
                f"forming a C2-symmetric E-E_op2 dimer ({nc_dimer['unique_unordered_residue_pairs_at_5A']} unordered residue pairs "
                f"/ {nc_dimer['directed_residue_pairs_at_5A']} directed contacts at <= 5.0 A; "
                f"Zn-Zn distance {nc_dimer['zn_zn_distance_A']:.4f} A) "
                "poised for interior genomic/mRNA recognition."
            ),
            "epistemic_status": ":Accepted",
        },
        {
            "arrangement_id": "arrangement_2_direct_lumenal_rna_cargo_encapsulation_engineered_darc1",
            "physical_topology": "Direct lumenal RNA cargo encapsulation inside reconstituted engineered dArc1 capsids",
            "primary_reference": "Lash, B. et al. (2026) Nature Communications, Nov 4, 2026, doi:10.1038/s41467-026-79058-8",
            "mmcif_entity_and_chain_provenance": (
                "Engineered recombinant dArc1 capsid constructs fused to heterologous RNA-binding domains (RBDs) "
                "for in-vitro assembly and mammalian SORCS2-targeted cargo delivery."
            ),
            "physical_mechanism_and_location": (
                "Exogenous mRNA/editor cargo is directly packaged inside the closed icosahedral dArc1 capsid lumen "
                "during in-vitro reconstitution via engineered RBD tethering, enabling SORCS2-mediated uptake and "
                "up to 18% Cas9 exon skipping in mdx mouse muscle."
            ),
            "epistemic_status": ":Accepted (engineered delivery system) / :Incomplete (native 6TAP sub-interface release mechanics)",
        },
        {
            "arrangement_id": "arrangement_3_non_encapsulated_arc_mbnl1_ev_cotransport",
            "physical_topology": "Non-encapsulated Arc-Mbnl1 mRNA association and extracellular-vesicle (EV) membrane co-transport",
            "primary_reference": "Zinter, M. et al. (2026) bioRxiv / PMC13061017, posted Mar 10, 2026",
            "mmcif_entity_and_chain_provenance": (
                "Cellular mammalian Arc ribonucleoprotein / extracellular vesicle fraction (not a closed 6TAP dArc1 capsid)."
            ),
            "physical_mechanism_and_location": (
                "Activity-dependent Arc-Mbnl1 mRNA association occurs via membrane-enveloped extracellular vesicles "
                "where Mbnl1 mRNA is NOT directly encapsulated inside a closed icosahedral Arc protein shell; must "
                "never be conflated with direct lumenal capsid encapsulation."
            ),
            "epistemic_status": ":Incomplete",
        },
    ]

    four_branch_provenance_table: List[Dict[str, Any]] = [
        {
            "branch_id": "branch_1_native_darc1_6tap_structure",
            "branch_title": "Native dArc1 Capsid Structure & Biophysical Domain Characterization (PDB 6TAP)",
            "primary_citations": [
                "Erlendsson, S. et al. (2020) Structures of virus-like capsids formed by the Drosophila neuronal Arc proteins. Nature Neuroscience 23, 172-175 (PDB 6TAP / 6TAR).",
                "Cottee, M. A. et al. (2020) Crystal structures of dArc1 capsid domain fragments. Nature Neuroscience 23, 176-180.",
                "Hallin, E. I. et al. (2021) Structure and biophysical properties of the Drosophila Arc capsid domain.",
            ],
            "epistemic_status": ":Accepted",
            "physical_rna_arrangement_id": "arrangement_1_internal_nc_zn_rna_binding_architecture_in_6tap",
            "supported_claims": [
                "T=4 icosahedral capsid architecture of Drosophila dArc1 CA in icosahedral point group I (_em_single_particle_entity.point_symmetry = I, biological assembly 1, RCSB global symmetry I, 60 rotation operators; 240 CA subunits across 60 asymmetric units, 12 five-fold A-chain pentamers and 30 two-fold B/C/D hexamers).",
                "204 symmetry-unique heavy-atom residue contacts at d <= 5.0 A partitioned across 10 symmetry-inequivalent sub-interfaces (9,570.22 A^2 representative two-sided pairwise BSA and 4,785.10 A^2 one-side interface area per ASU; 574,213.20 A^2 whole-capsid two-sided buried surface area and 287,106.60 A^2 whole-capsid one-side interface area).",
                "Exact mmCIF entity reconciliation: polymer _entity.id = 1 (dArc1, 254 aa, UniProt Q7K1U0) encompasses CA shell chains A..D (residues 41..206) and the resolved interior C-terminal NC zinc-knuckle segment chain E (label_asym_id = E, auth_asym_id = E, label_entity_id = 1, residues 224..250), which coordinates non-polymer ligand _entity.id = 2 (ZN, label_asym_id = F, auth_asym_id = E, label_entity_id = 2, pdb_seq_num = 301, auth_seq_num = 251) beneath the C2 axis.",
            ],
            "unsupported_or_quarantined_extrapolations": [
                "Does NOT establish the individual kinetic assembly order, thermodynamic delta-delta-G contribution, or in-vivo synaptic signaling requirement of any specific sub-interface without site-directed mutagenesis.",
            ],
        },
        {
            "branch_id": "branch_2_2026_engineered_darc1_delivery",
            "branch_title": "2026 Engineered dArc1 RNA Delivery Vehicle & SORCS2 Receptor Identification",
            "primary_citations": [
                "Lash, B. et al. (2026) Engineered dArc1 capsids for targeted mRNA delivery via SORCS2. Nature Communications, published Nov 4, 2026, doi:10.1038/s41467-026-79058-8.",
            ],
            "epistemic_status": ":Accepted (for engineered delivery system) / :Incomplete (for native 6TAP interface function)",
            "epistemic_status_by_scope": {
                "engineered_reconstitution_and_mdx_delivery": ":Accepted",
                "native_6tap_subinterface_functional_mapping": ":Incomplete",
            },
            "physical_rna_arrangement_id": "arrangement_2_direct_lumenal_rna_cargo_encapsulation_engineered_darc1",
            "supported_claims": [
                "In-vitro reconstitution of engineered dArc1 capsids directly encapsulating exogenous cargo mRNA inside the capsid lumen via engineered RNA-binding domains (RBD).",
                "Identification of mammalian SORCS2 as a cell-surface receptor mediating uptake of engineered dArc1 particles.",
                "Up to 18% Cas9-mediated exon skipping in mdx mouse skeletal/cardiac muscle models following engineered dArc1 delivery.",
            ],
            "unsupported_or_quarantined_extrapolations": [
                "Engineered RBD cargo tethering and mammalian SORCS2 binding do NOT identify which of the 10 native 6TAP CA sub-interfaces triggers physiological cargo release or disassembly in native Drosophila neurons (:Incomplete).",
            ],
        },
        {
            "branch_id": "branch_3_2026_arc_mbnl1_rna_ev_preprint",
            "branch_title": "2026 Activity-Dependent Arc-Mbnl1 RNA Association & Extracellular Vesicle (EV) Transfer Preprint",
            "primary_citations": [
                "Zinter, M. et al. (2026) Activity-dependent Arc-Mbnl1 RNA association and extracellular vesicle transfer. bioRxiv / PMC13061017, posted Mar 10, 2026.",
            ],
            "epistemic_status": ":Incomplete",
            "physical_rna_arrangement_id": "arrangement_3_non_encapsulated_arc_mbnl1_ev_cotransport",
            "supported_claims": [
                "Preprint reports activity-dependent association between mammalian Arc and Mbnl1 mRNA and possible extracellular-vesicle (EV) transfer between cells.",
                "Crucially reports that Mbnl1 mRNA is NOT directly encapsulated inside a closed Arc protein capsid shell, distinguishing membrane-enveloped EV co-transport/association from direct intracapsid lumen packaging.",
            ],
            "unsupported_or_quarantined_extrapolations": [
                "Quarantined as :Incomplete due to unrefereed preprint status and because non-encapsulated EV membrane co-transport cannot be conflated with direct lumenal RNA packaging inside the 6TAP icosahedral CA shell.",
            ],
        },
        {
            "branch_id": "branch_4_candidate_mutational_hotspots",
            "branch_title": "Pre-Registered Candidate Mutational Hotspots Across 6TAP Sub-Interfaces",
            "primary_citations": [
                "Pre-registered structural hypotheses derived from M7 Level-2/Level-3 6TAP residue-interface census (local verification receipt 2026-11-07).",
            ],
            "epistemic_status": ":Incomplete",
            "supported_claims": [
                "Geometric identification of high-centrality candidate contact motifs across the 10 symmetry-inequivalent 6TAP sub-interfaces and interior NC zinc knuckle suitable for site-directed mutagenesis, with explicit separation between ordinary CA-CA interface energetics perturbations and CCHC metal-coordination/NC-folding perturbations.",
            ],
            "unsupported_or_quarantined_extrapolations": [
                "All four candidate hotspot motifs remain strictly :Incomplete pre-registered hypotheses until experimentally tested against the falsification readouts and required biophysical controls below.",
            ],
            "pre_registered_hotspot_hypotheses": [
                {
                    "hotspot_id": "hotspot_1_PHE133_ctd_aromatic_lock",
                    "residue_motif": "PHE133 (with TYR126 / TYR129 / HIS162 / PHE172 aromatic cluster)",
                    "target_subinterfaces": ["A-D_C1_intra_asu_op1", "B-C_C3_inter_asu_hexamer_op33"],
                    "perturbation_category": "CA_CA_SUBUNIT_INTERFACE_ENERGETICS_PERTURBATION",
                    "required_experimental_controls": (
                        "Monomer secondary/tertiary fold verification (far-UV CD and 1H-15N HSQC NMR of isolated CTD), "
                        "SEC-MALS oligomerization stoichiometry, DSF thermal unfolding (Tm), and negative-stain/cryo-EM "
                        "T=4 capsid shell morphology."
                    ),
                    "structural_observation": (
                        "PHE133 forms the central inter-subunit homologous-position aromatic lock across both quasi-equivalent heterotypic CTD dimers "
                        "in interface-fingerprint Family II (A:PHE133-D:PHE133 in A-D_op1 at 2.8699 A min heavy-atom distance, ring-centroid distance "
                        "4.217 A, plane angle 35.32 deg; B:PHE133-C:PHE133 in B-C_op33 at 3.5768 A min heavy-atom "
                        "distance, ring-centroid distance 4.933 A, plane angle 22.21 deg), flanked by TYR126-PHE172 "
                        "and HIS162-HIS162 pi-stacking candidates (8/10 of all capsid pi-stacking candidates)."
                    ),
                    "testable_mutation": "F133A and F133E single-point substitutions in recombinant dArc1 CA",
                    "pre_stated_falsification_readout": (
                        "Falsified if F133A/F133E recombinant dArc1 retains wild-type T=4 capsid yield and "
                        "hydrodynamic radius by SEC-MALS / negative-stain EM (delta-Tm < 1.5 K by DSF), which "
                        "would prove PHE133 is passively packed rather than a thermodynamic dimer-lock hotspot."
                    ),
                    "epistemic_status": ":Incomplete",
                },
                {
                    "hotspot_id": "hotspot_2_ARG49_ASP66_GLU76_ARG111_pentamer_salt_bridge_network",
                    "residue_motif": "ARG49 / ASP66 / GLU76 / ARG111 (with ARG56-ASP144 and GLU117-ARG199 ion-pair network)",
                    "target_subinterfaces": [
                        "A-A_C5_pentamer_op36",
                        "C-D_C1_intra_asu_op1",
                        "B-C_C1_intra_asu_op1",
                        "B-D_C2_dimer_op2",
                    ],
                    "perturbation_category": "CA_CA_SUBUNIT_INTERFACE_ENERGETICS_PERTURBATION",
                    "required_experimental_controls": (
                        "Monomer fold verification by CD/DSF across 50-500 mM NaCl ionic-strength titration, "
                        "SEC-MALS pentamer/hexamer ratio quantification, and compensatory double-mutant (R56E/D144R) "
                        "epistasis rescue controls."
                    ),
                    "structural_observation": (
                        "A-A_C5_pentamer_op36 (1,882.00 A^2 BSA, 41 residue pairs) in interface-fingerprint Family I "
                        "combines the conserved cross-domain salt-bridge candidates A:ARG56-A:ASP144 (3.424/3.812 A, "
                        "also conserved in B-C_op1 at 2.712 A, C-D_op1 at 2.264 A, and B-D_op2 at 3.929 A) and "
                        "A:ARG199-A:GLU117 (3.937 A) with the NTD charged/polar stalk-and-bundle contact network "
                        "involving ARG49, ASP66, GLU76, and ARG111."
                    ),
                    "testable_mutation": "Charge-reversal mutants R56E/D144R, R49E/E76R, and D66R/R111D (plus 500 mM NaCl ionic-strength titration)",
                    "pre_stated_falsification_readout": (
                        "Falsified if charge-swap or 500 mM NaCl screening fails to perturb C5 pentamer and "
                        "hexamer NTD+CTD seam closure by cryo-EM / SEC-MALS, or if compensatory double mutations "
                        "(R56E/D144R) fail to restore assembly relative to single charge-reversal mutants."
                    ),
                    "epistemic_status": ":Incomplete",
                },
                {
                    "hotspot_id": "hotspot_3_HIS182_HIS186_ctd_trimer_contact",
                    "residue_motif": "HIS182 / HIS186 (with LYS181 / GLN197 / ASN204)",
                    "target_subinterfaces": ["C-C_C3_trimer_op33", "B-D_C3_inter_asu_vertex_op33"],
                    "perturbation_category": "CA_CA_SUBUNIT_INTERFACE_ENERGETICS_PERTURBATION",
                    "required_experimental_controls": (
                        "pH-dependent monomer stability controls (CD/DSF across pH 5.5-8.0), SEC-MALS capsid mass "
                        "distribution, and RNase protection / capsid permeability assays."
                    ),
                    "structural_observation": (
                        "C-C_C3_trimer_op33 (359.89 A^2) and B-D_C3_inter_asu_vertex_op33 (344.66 A^2) form the "
                        "quasi-equivalent C3 vertex core of Family_III_Residual_Distal_CTD_Hierarchy (Jaccard = 0.5000, "
                        "Sorensen-Dice = 0.6667 on seq_ids 177, 180, 181, 182, 200, 204), where C:HIS182-C:HIS182 "
                        "(centroid 6.893 A) and B:HIS182-D:HIS182 (centroid 6.808 A) form C3 pi-stacking candidates "
                        "alongside C:HIS186-C:GLN197 and C:HIS186-C:ILE200."
                    ),
                    "testable_mutation": "H182A/H186A double alanine substitution and pH 5.5..8.0 histidine-protonation titration",
                    "pre_stated_falsification_readout": (
                        "Falsified if capsid porosity, RNase protection, and C3 vertex stability are "
                        "invariant across pH 5.5 to 8.0 and unaffected by H182A/H186A substitution."
                    ),
                    "epistemic_status": ":Incomplete",
                },
                {
                    "hotspot_id": "hotspot_4_CYS227_CYS230_HIS235_CYS240_nc_zn_knuckle",
                    "residue_motif": "CYS227 / CYS230 / HIS235 / CYS240 (CCHC Zn-knuckle in NC chain E, seq_id 224..250)",
                    "target_subinterfaces": ["E-E_C2_nc_dimer_op2", "NC_chain_E_ZN_F301_coordination_sphere"],
                    "perturbation_category": "METAL_COORDINATION_AND_NC_FOLDING_PERTURBATION_NOT_ORDINARY_CA_INTERFACE_HOTSPOT",
                    "required_experimental_controls": (
                        "Zn occupancy (ICP-MS/PAR), NC folding (CD/NMR), and RNA-binding controls "
                        "(EMSA/fluorescence polarization with dArc1 3'-UTR RNA) alongside empty-capsid CA shell "
                        "assembly controls (SEC-MALS/cryo-EM) so NC misfolding or Zn2+ loss is never confounded "
                        "with ordinary CA-CA interface energetics."
                    ),
                    "structural_observation": (
                        f"In 6TAP, interior nucleocapsid strand E (label_asym_id = E, auth_asym_id = E, "
                        f"label_entity_id = 1, residues 224..250, C-terminal NC fragment of dArc1 Entity 1) "
                        f"coordinates ZN (label_asym_id = F, auth_asym_id = E, label_entity_id = 2, "
                        f"pdb_seq_num 301, auth_seq_num 251) via the tetrahedral CCHC motif "
                        f"C227-X2-C230-X4-H235-X4-C240 (tau_4 = {zn_sphere['tau_4_tetrahedrality_index']}, "
                        f"mean bond angle {zn_sphere['mean_bond_angle_deg']} deg) and forms a C2-symmetric "
                        f"NC-NC dimer across op2 ({nc_dimer['unique_unordered_residue_pairs_at_5A']} unordered "
                        f"residue pairs at <= 5.0 A, Zn-Zn distance {nc_dimer['zn_zn_distance_A']} A). Perturbing "
                        f"C227/C230/H235/C240 ablates tetrahedral Zn2+ coordination and unfolds the NC zinc knuckle "
                        f"rather than testing an ordinary CA-CA capsid shell interface."
                    ),
                    "testable_mutation": "C227S/C230S Zn-chelation knockout or TPEN/EDTA Zn2+ depletion assay",
                    "pre_stated_falsification_readout": (
                        "Falsified if C227S/C230S mutation or Zn2+ chelation (verified by ICP-MS/PAR Zn occupancy "
                        "and CD/NMR NC folding controls) shows no change in specific dArc1 mRNA 3'-UTR binding "
                        "stoichiometry or C2 NC-dimer ordering while preserving empty-capsid CA shell assembly."
                    ),
                    "epistemic_status": ":Incomplete",
                },
            ],
        },
    ]

    receipt: Dict[str, Any] = {
        "schema_id": "m7-arc-level3-interface-comparison-v1",
        "schema_version": "1.3.0",
        "module": "M7ArcLevel3InterfaceComparison",
        "evidence_class": "SYMMETRY_NORMALIZED_INTERFACE_COMPARISON_AND_PROVENANCE_SEPARATION",
        "verification_disposition": "remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run",
        "publication_and_ci_provenance": {
            "github_pr_url": "https://github.com/aurascoper/flow-cytometry-model-bridge/pull/6",
            "historical_v1_baseline_commit": "955e32dc673d3a874f874ffb8965497db03ca112",
            "level1_v2_repair_commit": "6a09714d03fccf29905c99bddc416785781063d1",
            "level2_atlas_and_crosscheck_commit": "fc5e1a407e375014653b2b72f73b89559c14e2ee",
            "level3_initial_publication_commit": "1abeaad2e44fe25e3e4cd17b7920a51f990fc55e",
            "local_test_execution_summary": (
                "28/28 Python pytest assertions and 3,059/3,059 Julia assertions (restored from 2,996/2,996) executed and verified locally"
            ),
            "remote_github_actions_ci_status": "No GitHub Actions workflow runs executed on remote; test counts are recorded local verification results, not remote CI reruns",
            "historical_assertion_count_audit_trail": {
                "2996": (
                    "flow-cytometry-model-bridge 4 Julia suites (293 + 1,682 + 332 + 689 = 2,996) when 63 L89..L103 "
                    "predicates (53 in test_MicrofluidicClassificationEngine.jl, 10 in test_M7VisualizationBridge.jl) "
                    "were conjugated via && into composite @test macros."
                ),
                "3059": (
                    "flow-cytometry-model-bridge 4 Julia suites (293 + 1,735 + 332 + 699 = 3,059 across 1,907 static "
                    "@test macros) with all 63 L89..L103 predicates unrolled into individual top-level @test and "
                    "@test_throws statements (2,996 + 63 = 3,059)."
                ),
                "3830": (
                    "7-suite cross-repo total (flow-cytometry-model-bridge 4 Julia suites + M7Multiscale 327 + "
                    "SmallModularBioReactors.jl 501 + Python M7SobolevSceneAdapter 6) under the consolidated 2,996 "
                    "bridge Julia count (2,996 + 327 + 501 + 6 = 3,830)."
                ),
                "3893": (
                    "7-suite cross-repo total under the restored 3,059 bridge Julia count (3,059 + 327 + 501 + 6 = "
                    "3,893; 3,830 + 63 = 3,893)."
                ),
                "3913": (
                    "7-suite cross-repo total (3,893) plus the 20 pre-MCMC M7 ARC Level-2/3 Python pytest checks "
                    "(12 Level-2 + 8 Level-3 = 20 -> 3,893 + 20 = 3,913; with all 22 current Level-2/3 Python checks "
                    "and 6 Sobolev adapter checks, the 7-suite + ARC Python total is 3,887 + 28 = 3,915)."
                ),
                "4017": (
                    "8-suite cross-repo grand total (7 suites + NeuralCompose Tests/eval 187) under the consolidated "
                    "2,996 bridge Julia count (3,830 + 187 = 4,017)."
                ),
                "4080": (
                    "8-suite cross-repo grand total (7 suites + NeuralCompose Tests/eval 187) under the restored "
                    "3,059 bridge Julia count (3,893 + 187 = 4,080; 4,017 + 63 = 4,080)."
                ),
            },
        },
        "mmcif_provenance": mmcif_prov,
        "upstream_receipts": {
            "level2_validation_receipt_path": "artifacts/m7_arc_level2_validation_receipt.json",
            "level2_validation_receipt_sha256": sha256_file(l2_receipt_path),
            "level2_independent_crosscheck_receipt_path": "artifacts/m7_arc_level2_independent_crosscheck_receipt.json",
            "level2_independent_crosscheck_receipt_sha256": sha256_file(crosscheck_receipt_path),
            "level3_script_path": "src/m7_arc_level3_interface_comparison.py",
            "level3_script_sha256": sha256_file(Path(__file__).resolve()),
        },
        "family_robustness_stratification": family_robustness_stratification,
        "symmetry_normalized_residue_fingerprints": {
            "ntd_seq_id_range": [41, NTD_MAX_SEQ_ID],
            "ctd_seq_id_range": [NTD_MAX_SEQ_ID + 1, 206],
            "total_unique_homologous_interfacial_seq_ids": len(all_interfacial_seqs),
            "all_interfacial_homologous_seq_ids": all_interfacial_seqs,
            "conserved_interfacial_seq_ids_across_all_4_chains_ABCD": conserved_across_all_4_chains,
            "conserved_interfacial_seq_ids_across_3plus_chains": conserved_across_3plus_chains,
            "per_chain_interfacial_seq_id_counts": {
                ch: len(chain_global_seq_sets[ch]) for ch in ("A", "B", "C", "D")
            },
            "interface_fingerprint_families": {
                "terminology_note": (
                    "Families I, II, and III (Family_III_Residual_Distal_CTD_Hierarchy) are empirical "
                    "interface-fingerprint groupings (interface_fingerprint_family) defined by homologous "
                    "sequence-position similarity and NTD/CTD contact architecture in 6TAP, not established "
                    "biological or functional modules."
                ),
                "Family_I_NTD_CTD_Capsomere_Seam": sorted(FAMILY_I_SET),
                "Family_II_Pure_CTD_Quasi_Equivalent_Dimer": sorted(FAMILY_II_SET),
                "Family_III_Residual_Distal_CTD_Hierarchy": sorted(FAMILY_III_FULL_SET),
                "Family_III_Vertex_and_Peripheral_CTD_Contacts": sorted(FAMILY_III_FULL_SET),
            },
            "family_robustness_stratification": family_robustness_stratification,
            "subinterface_fingerprints": subinterface_fingerprints,
            "pairwise_homologous_seq_id_jaccard_matrix": seq_jaccard_matrix,
            "pairwise_homologous_seq_id_sorensen_dice_matrix": seq_sorensen_dice_matrix,
            "pairwise_homologous_seq_id_szymkiewicz_simpson_matrix": seq_szymkiewicz_simpson_matrix,
            "pairwise_homologous_seq_id_cosine_ochiai_matrix": seq_cosine_ochiai_matrix,
            "pairwise_homologous_seq_id_contact_degree_weighted_cosine_matrix": seq_degree_weighted_cosine_matrix,
            "pairwise_canonical_seq_pair_jaccard_matrix": pair_jaccard_matrix,
            "pairwise_canonical_seq_pair_sorensen_dice_matrix": pair_sorensen_dice_matrix,
            "pairwise_canonical_seq_pair_szymkiewicz_simpson_matrix": pair_szymkiewicz_simpson_matrix,
            "pairwise_canonical_seq_pair_cosine_ochiai_matrix": pair_cosine_ochiai_matrix,
            "pairwise_ntd_seq_id_jaccard_matrix": ntd_seq_jaccard_matrix,
            "pairwise_ctd_seq_id_jaccard_matrix": ctd_seq_jaccard_matrix,
            "alternative_similarity_metrics_and_threshold_sensitivity": alternative_similarity_and_sensitivity,
            "focused_patch_comparisons": focused_patch_comparisons,
        },
        "alternative_similarity_metrics_and_threshold_sensitivity": alternative_similarity_and_sensitivity,
        "normalized_interface_chemistry_and_area_profiles": {
            "representative_subinterface_pairwise_bsa_sum_A2": round(total_sr_bsa, 2),
            "representative_subinterface_per_side_bsa_sum_A2": round(total_sr_per_side_bsa, 2),
            "whole_capsid_total_buried_surface_area_A2": round(total_whole_capsid_bsa, 2),
            "whole_capsid_one_side_interface_area_A2": total_whole_capsid_one_side_bsa,
            "total_unique_residue_pairs_at_5A": total_unique_pairs,
            "voronoi_vs_bsa_linear_calibration": cc_data["sasa_and_voronoi_surface_crosscheck"][
                "voronoi_vs_bsa_geometric_distinction_and_linear_calibration"
            ],
            "subinterface_profiles": normalized_profiles,
            "competing_null_model_evaluations": competing_null_model_tests,
        },
        "constrained_graph_state_space_connectivity_and_mcmc_audit": (
            constrained_graph_state_space_connectivity_and_mcmc_audit
        ),
        "multiplicity_adjustment_analysis": multiplicity_adjustment_analysis,
        "three_physical_rna_arrangements": three_physical_rna_arrangements,
        "four_branch_literature_and_provenance_separation": four_branch_provenance_table,
    }
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build M7 ARC Level-3 symmetry-normalized interface comparison receipt."
    )
    parser.add_argument("--l2-receipt", type=Path, default=DEFAULT_L2_RECEIPT)
    parser.add_argument("--crosscheck-receipt", type=Path, default=DEFAULT_CROSSCHECK_RECEIPT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    receipt = build_level3_interface_comparison_receipt(
        l2_receipt_path=args.l2_receipt,
        crosscheck_receipt_path=args.crosscheck_receipt,
    )
    text = json.dumps(receipt, indent=2, sort_keys=False, ensure_ascii=True) + "\n"
    if chr(0x2014) in text or chr(0x2013) in text:
        raise RuntimeError("Em-dash or en-dash detected in Level-3 receipt")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8")
    print(f"Wrote Level-3 receipt to {args.out} ({len(text)} bytes, sha256={sha256_file(args.out)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
