#!/usr/bin/env python3
"""
Unit and verification tests for:
  - Step 1: Repaired Level-1 Interpretation Ledger (`B and D both lack homotypic interfaces`)
  - Step 2: Level-2 BSA normalization (`9570.22`, `4785.10`, `574213.20`) & candidate chemistry qualification
  - Step 3: Independent Structural & Surface Cross-Check (`m7_arc_level2_independent_crosscheck.py`)
  - Steps 4, 5 & 6: Level-3 Symmetry-Normalized Interface Comparison & 4-Branch Provenance Separation
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from m7_arc_level3_interface_comparison import (  # noqa: E402
    CANONICAL_SUBINTERFACE_ORDER,
    build_level3_interface_comparison_receipt,
)


class TestM7ArcLevel3InterfaceComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.l1_followup_path = REPO_ROOT / "artifacts" / "arc_level1_followup_receipt.json"
        cls.l2_receipt_path = REPO_ROOT / "artifacts" / "m7_arc_level2_validation_receipt.json"
        cls.cc_receipt_path = REPO_ROOT / "artifacts" / "m7_arc_level2_independent_crosscheck_receipt.json"
        cls.l3_receipt_path = REPO_ROOT / "artifacts" / "m7_arc_level3_interface_comparison_receipt.json"
        cls.l2_spec_path = REPO_ROOT / "docs" / "M7_ARC_LEVEL2_RESIDUE_SPEC.md"
        cls.l3_spec_path = REPO_ROOT / "docs" / "M7_ARC_LEVEL3_INTERFACE_COMPARISON_SPEC.md"

        cls.l1_followup = json.loads(cls.l1_followup_path.read_text(encoding="utf-8"))
        cls.l2_receipt = json.loads(cls.l2_receipt_path.read_text(encoding="utf-8"))
        cls.cc_receipt = json.loads(cls.cc_receipt_path.read_text(encoding="utf-8"))
        cls.l3_receipt = json.loads(cls.l3_receipt_path.read_text(encoding="utf-8"))

    def test_01_step1_level1_homotypic_contradiction_repaired(self) -> None:
        l1_text = self.l1_followup_path.read_text(encoding="utf-8")
        l2_spec_text = self.l2_spec_path.read_text(encoding="utf-8")

        forbidden = (
            "B is uniquely free of " + "homotypic interfaces",
            "A/C/D classes each include " + "same-chain partner contacts",
        )
        for phrase in forbidden:
            self.assertNotIn(phrase, l1_text)
            self.assertNotIn(phrase, l2_spec_text)

        required_statement = (
            "B and D both lack homotypic interfaces in the 6TAP residue-contact census "
            "(B-B = 0 and D-D = 0 across 3.5..6.0 A; only A-A = 41 and C-C = 7 are homotypic at 5.0 A)"
        )
        self.assertIn(required_statement, l1_text)
        self.assertIn(required_statement, l2_spec_text)
        self.assertEqual(
            self.l1_followup["earlier_acd_same_class_contact_rationale_status"],
            ":Superseded",
        )
        self.assertIn(
            "1,287.13 A^2",
            self.l1_followup["level2_re_evaluation_of_k2_B_vs_ACD_partition"]["heterotypic_partner_and_bsa_profile"],
        )

    def test_02_step2_bsa_normalization_and_candidate_chemistry(self) -> None:
        l2_b = self.l2_receipt["work_packages"]["L2_B_sasa_buried_area"]
        self.assertAlmostEqual(l2_b["representative_subinterface_pairwise_bsa_sum_A2"], 9570.22, places=2)
        self.assertAlmostEqual(l2_b["representative_subinterface_per_side_bsa_sum_A2"], 4785.10, places=2)
        self.assertAlmostEqual(l2_b["whole_capsid_total_buried_surface_area_A2"], 574213.20, places=2)
        self.assertAlmostEqual(l2_b["whole_capsid_one_side_interface_area_A2"], 287106.60, places=2)

        for sub in l2_b["interface_buried_surface_area_by_sub_interface"]:
            self.assertEqual(sub["unordered_capsid_edge_multiplicity"], 60)
            self.assertAlmostEqual(
                sub["whole_capsid_sub_interface_buried_surface_area_A2"],
                round(sub["total_buried_surface_area_A2"] * 60, 2),
                places=2,
            )
            self.assertAlmostEqual(
                sub["whole_capsid_one_side_bsa_A2"],
                round(sub["buried_surface_area_per_side_A2"] * 60, 2),
                places=2,
            )

        l2_c = self.l2_receipt["work_packages"]["L2_C_interface_chemistry"]
        self.assertIn("candidate_qualification_caveats", l2_c)
        for entry in l2_c["summary_by_sub_interface"]:
            for key in (
                "hydrophobic_candidate_pairs",
                "pi_stacking_candidate_pairs",
                "hbond_heavy_atom_distance_candidate_pairs",
                "salt_bridge_heavy_atom_distance_candidate_pairs",
            ):
                self.assertIn(key, entry)
        for pi_item in l2_c["aromatic_pi_stacking_candidates"]:
            self.assertIn("ring_plane_angle_deg", pi_item)
            self.assertIn("min_heavy_atom_distance_A", pi_item)

        eq = self.l2_receipt["epistemic_qualifications"]
        self.assertEqual(
            eq["verification_disposition"],
            "remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run",
        )
        self.assertEqual(
            eq["publication_and_ci_provenance"]["github_pr_url"],
            "https://github.com/aurascoper/flow-cytometry-model-bridge/pull/6",
        )
        self.assertEqual(
            eq["publication_and_ci_provenance"]["level3_initial_publication_commit"],
            "1abeaad2e44fe25e3e4cd17b7920a51f990fc55e",
        )
        self.assertEqual(len(eq["m7_falsification_design_rules"]), 4)

    def test_03_step3_independent_crosscheck_100_percent_identity(self) -> None:
        self.assertEqual(
            self.cc_receipt["verification_disposition"],
            "remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run",
        )
        self.assertEqual(
            self.cc_receipt["publication_and_ci_provenance"]["github_pr_url"],
            "https://github.com/aurascoper/flow-cytometry-model-bridge/pull/6",
        )
        self.assertEqual(
            self.cc_receipt["independence_scope"],
            "COMPUTATIONAL_IMPLEMENTATION_INDEPENDENT_SAME_DEPOSITED_6TAP_COORDINATES",
        )
        self.assertEqual(self.cc_receipt["point_symmetry"], "I")

        bm = self.cc_receipt["block_matrix_contact_crosscheck"]
        exact_5a = bm["exact_5A_residue_pair_set_comparison"]
        self.assertEqual(exact_5a["independent_unique_pair_count_at_5A"], 204)
        self.assertEqual(exact_5a["level2_atlas_unique_pair_count_at_5A"], 204)
        self.assertEqual(exact_5a["intersection_count"], 204)
        self.assertEqual(exact_5a["union_count"], 204)
        self.assertEqual(exact_5a["symmetric_difference_count"], 0)
        self.assertAlmostEqual(exact_5a["jaccard_similarity"], 1.0, places=6)
        self.assertTrue(exact_5a["hundred_percent_identity_verified"])
        self.assertEqual(exact_5a["atom_pair_hit_count_mismatches_across_204_pairs"], 0)

        expected_totals = {
            "3.5": 50,
            "4.0": 79,
            "4.5": 136,
            "5.0": 204,
            "5.5": 284,
            "6.0": 371,
        }
        for c_str, tot in expected_totals.items():
            entry = bm["multi_cutoff_census_verification"][c_str]
            self.assertEqual(entry["independent_unique_total"], tot)
            self.assertEqual(entry["level2_atlas_unique_total"], tot)
            self.assertEqual(entry["independent_canonical_pair_count"], tot)
            self.assertEqual(entry["level2_atlas_canonical_pair_count"], tot)
            self.assertEqual(entry["intersection_count"], tot)
            self.assertEqual(entry["union_count"], tot)
            self.assertEqual(entry["symmetric_difference_count"], 0)
            self.assertTrue(entry["all_10_chain_pairs_exact_match"])
            self.assertAlmostEqual(entry["jaccard_similarity_with_level2_atlas"], 1.0, places=6)

        surf = self.cc_receipt["sasa_and_voronoi_surface_crosscheck"]
        sums = surf["representative_subinterface_bsa_sums_A2"]
        self.assertAlmostEqual(sums["shrake_rupley_fibonacci_pairwise_bsa_sum_A2"], 9570.22, places=2)
        self.assertAlmostEqual(sums["lee_richards_zslice_pairwise_bsa_sum_A2"], 9400.64, places=2)
        self.assertAlmostEqual(sums["voronoi_radical_plane_contact_area_sum_A2"], 4236.73, places=2)
        self.assertGreater(sums["pearson_r_shrake_rupley_vs_lee_richards"], 0.999)
        self.assertGreater(sums["pearson_r_shrake_rupley_vs_voronoi"], 0.998)
        self.assertEqual(len(surf["sub_interface_crosschecks"]), 10)

        cal = surf["voronoi_vs_bsa_geometric_distinction_and_linear_calibration"]
        ols2 = cal["ols_vs_two_sided_shrake_rupley_pairwise_bsa"]
        ols1 = cal["ols_vs_one_side_shrake_rupley_bsa"]
        self.assertAlmostEqual(ols2["slope"], 0.444876, places=4)
        self.assertAlmostEqual(ols2["r_squared"], 0.998367, places=4)
        self.assertLess(ols2["rmse_A2"], 12.0)
        self.assertLess(ols2["max_abs_residual_A2"], 25.0)
        self.assertAlmostEqual(ols1["slope"], 0.889755, places=4)
        self.assertAlmostEqual(ols1["r_squared"], 0.998368, places=4)

        glob_asu = surf["global_asu_interfacial_residue_comparison"]
        self.assertEqual(glob_asu["distance_le_5p0A_residue_count"], 235)
        self.assertEqual(glob_asu["sasa_buried_gt_1p0A2_residue_count"], 282)
        self.assertEqual(glob_asu["sasa_buried_gt_0p1A2_residue_count"], 289)
        self.assertEqual(glob_asu["voronoi_contact_gt_0A2_residue_count"], 289)
        self.assertAlmostEqual(glob_asu["jaccard_sasa_1A2_vs_voronoi"], 0.9758, places=4)

        sub_by_cls = {s["interface_class"]: s for s in surf["sub_interface_crosschecks"]}
        self.assertEqual(
            sub_by_cls["C-C_C3_trimer_op33"]["interfacial_residue_counts"]["sasa_buried_gt_1p0A2_residues"],
            12,
        )
        self.assertAlmostEqual(
            sub_by_cls["C-C_C3_trimer_op33"]["interfacial_residue_set_agreement"]["jaccard_sasa_1A2_vs_voronoi"],
            1.0,
            places=4,
        )
        self.assertEqual(
            sub_by_cls["B-D_C3_inter_asu_vertex_op33"]["interfacial_residue_counts"]["sasa_buried_gt_1p0A2_residues"],
            12,
        )
        self.assertAlmostEqual(
            sub_by_cls["B-D_C3_inter_asu_vertex_op33"]["interfacial_residue_set_agreement"]["jaccard_sasa_1A2_vs_voronoi"],
            1.0,
            places=4,
        )

    def test_04_level3_symmetry_normalized_fingerprints_and_jaccard_matrices(self) -> None:
        self.assertEqual(
            self.l3_receipt["verification_disposition"],
            "remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run",
        )
        self.assertEqual(
            self.l3_receipt["publication_and_ci_provenance"]["github_pr_url"],
            "https://github.com/aurascoper/flow-cytometry-model-bridge/pull/6",
        )
        fp = self.l3_receipt["symmetry_normalized_residue_fingerprints"]
        self.assertEqual(fp["ntd_seq_id_range"], [41, 125])
        self.assertEqual(fp["ctd_seq_id_range"], [126, 206])
        self.assertEqual(len(fp["subinterface_fingerprints"]), 10)
        self.assertIn("interface_fingerprint_families", fp)
        self.assertIn("family_robustness_stratification", fp)
        strat = fp["family_robustness_stratification"]["stratification_by_family"]
        self.assertEqual(
            set(strat.keys()),
            {
                "Family_II_Pure_CTD_Quasi_Equivalent_Dimer",
                "Family_I_NTD_CTD_Capsomere_Seam",
                "Family_III_Residual_Distal_CTD_Hierarchy",
            },
        )
        self.assertEqual(strat["Family_II_Pure_CTD_Quasi_Equivalent_Dimer"]["robustness_rank"], 1)
        self.assertEqual(strat["Family_I_NTD_CTD_Capsomere_Seam"]["robustness_rank"], 2)
        self.assertEqual(strat["Family_III_Residual_Distal_CTD_Hierarchy"]["robustness_rank"], 3)
        self.assertEqual(
            strat["Family_II_Pure_CTD_Quasi_Equivalent_Dimer"][
                "homologous_seq_id_jaccard_continuous_interval"
            ],
            "(0.0476, 0.8333]",
        )
        self.assertEqual(
            strat["Family_II_Pure_CTD_Quasi_Equivalent_Dimer"][
                "canonical_seq_pair_jaccard_continuous_interval"
            ],
            "(0.0000, 0.6364]",
        )
        self.assertEqual(
            strat["Family_I_NTD_CTD_Capsomere_Seam"]["canonical_seq_pair_jaccard_continuous_intervals"][
                "single_linkage"
            ],
            "(0.0000, 0.2407]",
        )
        self.assertEqual(
            strat["Family_I_NTD_CTD_Capsomere_Seam"]["canonical_seq_pair_jaccard_continuous_intervals"][
                "complete_linkage"
            ],
            "(0.0000, 0.1698]",
        )
        self.assertFalse(
            strat["Family_III_Residual_Distal_CTD_Hierarchy"][
                "is_cohesive_four_member_cluster_at_complete_linkage_tau_ge_0p15"
            ]
        )
        self.assertEqual(
            len(strat["Family_III_Residual_Distal_CTD_Hierarchy"]["hierarchical_resolution"]),
            3,
        )
        tier2 = strat["Family_III_Residual_Distal_CTD_Hierarchy"]["hierarchical_resolution"][
            "tier_2_C5_inter_asu_seam"
        ]
        self.assertAlmostEqual(
            tier2["detachment_thresholds_canonical_seq_pair_jaccard"]["complete_linkage"], 0.0, places=4
        )
        self.assertAlmostEqual(
            tier2["detachment_thresholds_canonical_seq_pair_jaccard"]["single_linkage"], 0.0909, places=4
        )

        # Verify all 10 sub-interfaces sum to 204 unique residue pairs and carry interface_fingerprint_family
        total_pairs = sum(s["unique_residue_pair_count_at_5A"] for s in fp["subinterface_fingerprints"])
        self.assertEqual(total_pairs, 204)
        for s in fp["subinterface_fingerprints"]:
            self.assertIn(
                s["interface_fingerprint_family"],
                {
                    "Family_I_NTD_CTD_Capsomere_Seam",
                    "Family_II_Pure_CTD_Quasi_Equivalent_Dimer",
                    "Family_III_Residual_Distal_CTD_Hierarchy",
                },
            )
            for res in s["interfacial_residues"]:
                self.assertEqual(res["label_entity_id"], "1")
                self.assertEqual(res["struct_asym_entity_id"], "1")
                self.assertEqual(res["functional_component_id"], "CA_CAPSID_SHELL")
                self.assertIn("biological_domain_or_fragment_designation", res)

        seq_jacc = fp["pairwise_homologous_seq_id_jaccard_matrix"]
        pair_jacc = fp["pairwise_canonical_seq_pair_jaccard_matrix"]
        for ic1 in CANONICAL_SUBINTERFACE_ORDER:
            self.assertAlmostEqual(seq_jacc[ic1][ic1], 1.0, places=4)
            self.assertAlmostEqual(pair_jacc[ic1][ic1], 1.0, places=4)
            for ic2 in CANONICAL_SUBINTERFACE_ORDER:
                self.assertAlmostEqual(seq_jacc[ic1][ic2], seq_jacc[ic2][ic1], places=4)
                self.assertAlmostEqual(pair_jacc[ic1][ic2], pair_jacc[ic2][ic1], places=4)

        # Verify quasi-equivalent CTD dimers A-D_op1 and B-C_op33 have high Jaccard (0.8333 seq_id, 0.6364 pair)
        self.assertAlmostEqual(
            seq_jacc["A-D_C1_intra_asu_op1"]["B-C_C3_inter_asu_hexamer_op33"], 0.8333, places=4
        )
        self.assertAlmostEqual(
            pair_jacc["A-D_C1_intra_asu_op1"]["B-C_C3_inter_asu_hexamer_op33"], 0.6364, places=4
        )

        # Verify same-chain-pair distinct sub-interfaces use disjoint canonical residue pairs (Jaccard == 0.0)
        self.assertAlmostEqual(
            pair_jacc["B-C_C1_intra_asu_op1"]["B-C_C3_inter_asu_hexamer_op33"], 0.0, places=4
        )
        self.assertAlmostEqual(
            pair_jacc["A-D_C1_intra_asu_op1"]["A-D_C5_inter_asu_seam_op36"], 0.0, places=4
        )
        self.assertAlmostEqual(
            pair_jacc["B-D_C2_dimer_op2"]["B-D_C3_inter_asu_vertex_op33"], 0.0, places=4
        )

    def test_05_level3_normalized_chemistry_profiles_and_null_models(self) -> None:
        prof = self.l3_receipt["normalized_interface_chemistry_and_area_profiles"]
        self.assertAlmostEqual(prof["representative_subinterface_pairwise_bsa_sum_A2"], 9570.22, places=2)
        self.assertAlmostEqual(prof["representative_subinterface_per_side_bsa_sum_A2"], 4785.10, places=2)
        self.assertAlmostEqual(prof["whole_capsid_total_buried_surface_area_A2"], 574213.20, places=2)
        self.assertAlmostEqual(prof["whole_capsid_one_side_interface_area_A2"], 287106.60, places=2)
        self.assertEqual(prof["total_unique_residue_pairs_at_5A"], 204)
        self.assertEqual(len(prof["subinterface_profiles"]), 10)
        for sp in prof["subinterface_profiles"]:
            am = sp["area_metrics_A2"]
            self.assertAlmostEqual(
                am["whole_capsid_one_side_bsa_A2"],
                round(am["shrake_rupley_per_side_bsa_A2"] * 60, 2),
                places=2,
            )

        nulls = prof["competing_null_model_evaluations"]
        dom_null = nulls["domain_segregation_vs_random_mixing_null"]
        self.assertEqual(
            dom_null["inference_tier"],
            "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY",
        )
        self.assertEqual(len(dom_null["non_independence_caveats"]), 4)
        self.assertTrue(dom_null["null_rejected"])
        self.assertGreater(dom_null["chi_square_statistic_df2"], 20.0)

        chem_null = nulls["ctd_dimer_hydrophobic_aromatic_vs_ntd_seam_salt_bridge_null"]
        self.assertEqual(
            chem_null["inference_tier"],
            "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY",
        )
        self.assertEqual(len(chem_null["non_independence_caveats"]), 4)
        self.assertEqual(chem_null["total_hydrophobic_candidates_across_asu"], 20)
        self.assertEqual(chem_null["total_pi_stacking_candidates_across_asu"], 10)
        self.assertEqual(chem_null["total_salt_bridge_candidates_across_asu"], 9)
        self.assertGreater(
            chem_null["quasi_equivalent_ctd_dimers_A_D_op1_and_B_C_op33"]["hydrophobic_fold_enrichment_over_null"],
            2.0,
        )
        self.assertGreater(
            chem_null["quasi_equivalent_ctd_dimers_A_D_op1_and_B_C_op33"]["pi_stacking_fold_enrichment_over_null"],
            2.5,
        )

        perm_null = nulls["degree_and_interface_preserving_permutation_null"]
        self.assertEqual(
            perm_null["inference_tier"],
            "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY",
        )
        self.assertEqual(
            perm_null["epistemic_status"],
            ":Incomplete (exploratory structural screening until wet-lab mutagenesis)",
        )
        self.assertEqual(perm_null["n_permutations"], 10000)
        self.assertEqual(perm_null["random_seed"], 20261107)

        p1 = perm_null["domain_assignment_permutation_test"][
            "homologous_83_residue_degree_preserving_permutation"
        ]
        self.assertEqual(p1["observed_same_domain_contacts_NTD_NTD_plus_CTD_CTD"], 131)
        self.assertEqual(p1["permutations_with_same_domain_ge_observed"], 64)
        self.assertAlmostEqual(p1["empirical_p_value_same_domain_phipson_smyth"], 0.006499, places=5)
        self.assertEqual(p1["permutations_with_chi_square_ge_observed"], 455)
        self.assertAlmostEqual(
            p1["empirical_p_value_omnibus_chi_square_phipson_smyth"], 0.045595, places=5
        )

        p1_chain = perm_null["domain_assignment_permutation_test"][
            "chain_specific_235_node_degree_preserving_permutation"
        ]
        self.assertEqual(p1_chain["permutations_with_same_domain_ge_observed"], 0)
        self.assertAlmostEqual(p1_chain["empirical_p_value_same_domain_phipson_smyth"], 0.0001, places=5)

        p1_166 = perm_null["domain_assignment_permutation_test"][
            "full_166_resolved_sequence_domain_permutation"
        ]
        self.assertEqual(p1_166["permutations_with_same_domain_ge_observed"], 85)
        self.assertAlmostEqual(p1_166["empirical_p_value_same_domain_phipson_smyth"], 0.008599, places=5)

        p1_circ = perm_null["domain_assignment_permutation_test"][
            "contiguous_166_residue_circular_shift_autocorrelation_null"
        ]
        self.assertEqual(p1_circ["total_circular_shifts_evaluated"], 166)
        self.assertEqual(p1_circ["circular_shifts_with_aggregate_same_domain_ge_observed"], 38)
        self.assertAlmostEqual(
            p1_circ["empirical_p_circular_aggregate_same_domain"], 0.228916, places=5
        )
        self.assertEqual(p1_circ["circular_shifts_with_non_family_I_ctd_ctd_ge_observed"], 3)
        self.assertAlmostEqual(
            p1_circ["empirical_p_circular_non_family_I_ctd_purity"], 0.018072, places=5
        )

        p2_edge = perm_null["chemistry_clustering_permutation_test"][
            "subinterface_edge_assignment_permutation"
        ]
        self.assertEqual(
            p2_edge["family_II_hydrophobic_enrichment"]["permutations_ge_observed"],
            3,
        )
        self.assertEqual(
            p2_edge["family_II_pi_stacking_enrichment_at_5A"]["permutations_ge_observed"],
            2,
        )
        self.assertEqual(
            p2_edge["family_I_salt_bridge_enrichment_at_5A"]["permutations_ge_observed"],
            635,
        )

        p2_lbl = perm_null["chemistry_clustering_permutation_test"][
            "homologous_83_residue_comp_id_label_permutation"
        ]
        self.assertIn("topological_graph_fixed_caveat", p2_lbl)
        self.assertIn(
            "conditional on the fixed 204-edge contact graph and tests whether amino-acid labels are unusually distributed over already-observed edges",
            p2_lbl["topological_graph_fixed_caveat"],
        )
        self.assertEqual(
            p2_lbl["family_II_nonpolar_residue_pair_enrichment"]["permutations_ge_observed"],
            297,
        )
        self.assertEqual(
            p2_lbl["family_II_aromatic_residue_pair_enrichment"]["permutations_ge_observed"],
            434,
        )
        self.assertEqual(
            p2_lbl["family_I_opposite_charge_ion_pair_enrichment"]["permutations_ge_observed"],
            111,
        )

        # Commit 5: Explicit 7-row state_space_specification & 4-chain double-edge-swap MCMC audit
        mcmc = perm_null["constrained_graph_state_space_connectivity_and_mcmc_audit"]
        spec = mcmc["state_space_specification"]
        self.assertEqual(
            list(spec.keys()),
            [
                "vertex_identity",
                "graph_class",
                "invariants",
                "swap_constraints",
                "stationary_distribution",
                "connectivity",
                "mixing",
            ],
        )
        self.assertEqual(spec["vertex_identity"]["chain_residue_node_count"], 235)
        self.assertEqual(spec["vertex_identity"]["homologous_quotient_seq_id_count"], 83)
        self.assertEqual(spec["graph_class"]["edge_count_at_5A"], 204)
        self.assertEqual(
            spec["invariants"]["illustrative_pre_run_stub_counts_reconciliation"]["S_N_ntd_stubs"],
            169,
        )
        self.assertEqual(
            spec["invariants"]["illustrative_pre_run_stub_counts_reconciliation"]["S_C_ctd_stubs"],
            239,
        )
        self.assertTrue(spec["swap_constraints"]["lazy_metropolis_hastings_persistence_on_rejection"])
        self.assertIn(
            "lazy/Metropolis-Hastings persistence",
            spec["swap_constraints"]["summary"],
        )
        self.assertEqual(
            len(spec["connectivity"]["non_family_I_pure_ctd_subinterfaces"]),
            6,
        )
        self.assertEqual(
            len(spec["connectivity"]["zero_ntd_stub_subinterfaces_strict_125"]),
            4,
        )
        self.assertEqual(spec["connectivity"]["scheme_II_total_valid_simple_bipartite_swaps"], 2377)
        self.assertEqual(spec["connectivity"]["scheme_II_total_color_changing_swaps"], 48)
        self.assertEqual(spec["connectivity"]["scheme_II_theoretical_same_domain_corridor"], [123, 135])
        self.assertEqual(spec["mixing"]["n_independent_chains"], 4)
        self.assertEqual(spec["mixing"]["chain_seeds"], [20261108, 20261109, 20261110, 20261111])
        self.assertEqual(spec["mixing"]["burn_in_steps_per_chain"], 5000)
        self.assertEqual(spec["mixing"]["post_burn_in_steps_per_chain"], 50000)
        self.assertEqual(spec["mixing"]["total_post_burn_in_steps"], 200000)

        # Commit 5: Stub-difference conservation identity
        stub_id = mcmc["stub_difference_conservation_identity"]
        self.assertEqual(
            stub_id["role"],
            "CONSERVATION_INVARIANT_AND_IMPLEMENTATION_CHECK_NOT_EMPIRICAL_SIGNIFICANCE_STATISTIC",
        )
        self.assertEqual(
            stub_id["illustrative_35_stub_difference_parameterization_reconciliation"][
                "conserved_E_CC_minus_E_NN"
            ],
            35,
        )
        self.assertEqual(
            stub_id["empirical_6tap_5A_cutoff_values"]["conserved_E_CC_minus_E_NN"],
            19,
        )
        self.assertEqual(
            stub_id["empirical_6tap_5A_cutoff_values"]["scheme_I_mcmc_min_E_CC_minus_E_NN"],
            19,
        )
        self.assertEqual(
            stub_id["empirical_6tap_5A_cutoff_values"]["scheme_I_mcmc_max_E_CC_minus_E_NN"],
            19,
        )
        self.assertAlmostEqual(
            stub_id["empirical_6tap_5A_cutoff_values"]["scheme_I_mcmc_variance_E_CC_minus_E_NN"],
            0.0,
            places=8,
        )
        self.assertAlmostEqual(
            stub_id["empirical_6tap_5A_cutoff_values"]["scheme_II_mcmc_variance_E_CC_minus_E_NN"],
            0.0,
            places=8,
        )

        # Commit 5: Scheme I vs Scheme II MCMC results
        s1 = mcmc["global_degree_preserving_simple_graph_double_edge_swap_mcmc"]
        self.assertEqual(s1["steps_with_same_domain_ge_observed"], 46)
        self.assertAlmostEqual(s1["empirical_p_value_phipson_smyth"], 0.000235, places=6)
        self.assertAlmostEqual(
            s1["effective_sample_size_total"],
            954.15,
            places=2,
        )
        self.assertAlmostEqual(
            s1["gelman_rubin_r_hat"],
            1.00259,
            places=4,
        )
        self.assertAlmostEqual(
            s1["p_min_ess_floor"],
            0.001047,
            places=6,
        )

        s2 = mcmc["within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc"]
        self.assertEqual(s2["steps_with_same_domain_ge_observed"], 41852)
        self.assertAlmostEqual(s2["empirical_p_value_phipson_smyth"], 0.209264, places=6)
        self.assertEqual(s2["mcmc_min_same_domain_contacts"], 123)
        self.assertEqual(s2["mcmc_max_same_domain_contacts"], 135)
        self.assertAlmostEqual(
            s2["effective_sample_size_total"],
            1031.18,
            places=2,
        )
        self.assertAlmostEqual(
            s2["gelman_rubin_r_hat"],
            1.00159,
            places=4,
        )
        self.assertEqual(len(s2["per_subinterface_swap_mobility"]), 10)
        mob_list = s2["per_subinterface_swap_mobility"]
        self.assertEqual(
            [m["valid_simple_bipartite_swaps"] for m in mob_list],
            [695, 490, 303, 324, 269, 279, 10, 6, 1, 0],
        )
        self.assertEqual(
            [m["color_changing_valid_swaps"] for m in mob_list],
            [20, 12, 1, 0, 1, 14, 0, 0, 0, 0],
        )
        self.assertEqual(sum(m["chain_i_ntd_stubs"] for m in mob_list), 96)
        self.assertEqual(sum(m["chain_i_ctd_stubs"] for m in mob_list), 108)
        self.assertEqual(sum(m["chain_j_ntd_stubs"] for m in mob_list), 89)
        self.assertEqual(sum(m["chain_j_ctd_stubs"] for m in mob_list), 115)
        self.assertEqual(sum(m["theoretical_same_domain_corridor"][0] for m in mob_list), 123)
        self.assertEqual(sum(m["theoretical_same_domain_corridor"][1] for m in mob_list), 135)

        refs = mcmc["methodological_references_and_distinct_roles"]
        self.assertFalse(
            refs["topologically_invariant_permutation_test_2025"][
                "is_proof_of_m7_double_edge_swap_irreducibility_or_mixing"
            ]
        )

        pres = mcmc["monte_carlo_and_mcmc_p_value_resolution"]
        self.assertAlmostEqual(pres["p_min_perm_floor"], 1.0 / 10001.0, places=7)
        self.assertEqual(len(pres["resolution_and_uncertainty_summary_by_test"]), 14)

        # Commit 5: Multiplicity adjustment analysis
        mult = perm_null["multiplicity_adjustment_analysis"]
        self.assertEqual(mult["omnibus_12_permutation_test_battery"]["n_total_tests"], 12)
        self.assertEqual(
            mult["hypothesis_families"]["domain_segregation_and_spatial_autocorrelation_family"]["n_tests"],
            6,
        )
        self.assertEqual(
            mult["hypothesis_families"]["domain_segregation_and_spatial_autocorrelation_family"][
                "n_surviving_holm_within_family_alpha_0p05"
            ],
            3,
        )
        self.assertEqual(
            mult["hypothesis_families"]["domain_segregation_and_spatial_autocorrelation_family"][
                "n_surviving_bh_within_family_q_0p05"
            ],
            4,
        )
        self.assertEqual(
            mult["hypothesis_families"]["interface_chemistry_specialization_family"]["n_tests"],
            6,
        )
        self.assertEqual(
            mult["hypothesis_families"]["interface_chemistry_specialization_family"][
                "n_surviving_holm_within_family_alpha_0p05"
            ],
            3,
        )
        self.assertEqual(
            mult["hypothesis_families"]["interface_chemistry_specialization_family"][
                "n_surviving_bh_within_family_q_0p05"
            ],
            4,
        )
        self.assertEqual(
            mult["omnibus_12_permutation_test_battery"]["n_surviving_holm_omnibus_12_alpha_0p05"],
            3,
        )
        self.assertEqual(
            mult["omnibus_12_permutation_test_battery"]["n_surviving_bh_omnibus_12_q_0p05"],
            8,
        )

        ledger_path = REPO_ROOT / "VERIFICATION_LEDGER_2026-11-02.json"
        if ledger_path.exists():
            ledger_data = json.loads(ledger_path.read_text(encoding="utf-8"))
            l_rec = ledger_data["m7_arc_level2_residue_interface_atlas_reconciliation"]
            l_mcmc = l_rec["constrained_graph_state_space_connectivity_and_mcmc_audit_summary"]
            self.assertAlmostEqual(
                l_mcmc["scheme_I_global_mcmc"]["mean_same_domain"],
                s1["mcmc_mean_same_domain_contacts"],
                places=4,
            )
            self.assertEqual(
                l_mcmc["scheme_I_global_mcmc"]["k_ge_131"],
                s1["steps_with_same_domain_ge_observed"],
            )
            self.assertAlmostEqual(
                l_mcmc["scheme_I_global_mcmc"]["ess_total"],
                s1["effective_sample_size_total"],
                places=2,
            )
            self.assertAlmostEqual(
                l_mcmc["scheme_II_within_subinterface_mcmc"]["mean_same_domain"],
                s2["mcmc_mean_same_domain_contacts"],
                places=4,
            )
            self.assertEqual(
                l_mcmc["scheme_II_within_subinterface_mcmc"]["k_ge_131"],
                s2["steps_with_same_domain_ge_observed"],
            )
            l_mult = l_rec["multiplicity_adjustment_summary"]
            self.assertEqual(l_mult["family_A_holm_survived_count"], 3)
            self.assertEqual(l_mult["family_A_bh_fdr_survived_count"], 4)
            self.assertEqual(l_mult["family_B_holm_survived_count"], 3)
            self.assertEqual(l_mult["family_B_bh_fdr_survived_count"], 4)
            self.assertEqual(l_mult["omnibus_12_holm_survived_count"], 3)
            self.assertEqual(l_mult["omnibus_12_bh_fdr_survived_count"], 8)

    def test_06_four_branch_literature_and_provenance_separation(self) -> None:
        arrangements = self.l3_receipt["three_physical_rna_arrangements"]
        self.assertEqual(len(arrangements), 3)
        arr_ids = [a["arrangement_id"] for a in arrangements]
        self.assertEqual(
            arr_ids,
            [
                "arrangement_1_internal_nc_zn_rna_binding_architecture_in_6tap",
                "arrangement_2_direct_lumenal_rna_cargo_encapsulation_engineered_darc1",
                "arrangement_3_non_encapsulated_arc_mbnl1_ev_cotransport",
            ],
        )

        branches = self.l3_receipt["four_branch_literature_and_provenance_separation"]
        self.assertEqual(len(branches), 4)
        by_id = {b["branch_id"]: b for b in branches}

        b1 = by_id["branch_1_native_darc1_6tap_structure"]
        self.assertEqual(b1["epistemic_status"], ":Accepted")

        b2 = by_id["branch_2_2026_engineered_darc1_delivery"]
        self.assertIn("10.1038/s41467-026-79058-8", b2["primary_citations"][0])
        self.assertEqual(b2["epistemic_status_by_scope"]["engineered_reconstitution_and_mdx_delivery"], ":Accepted")
        self.assertEqual(b2["epistemic_status_by_scope"]["native_6tap_subinterface_functional_mapping"], ":Incomplete")
        self.assertTrue(any("SORCS2" in c for c in b2["supported_claims"]))
        self.assertTrue(any("18%" in c for c in b2["supported_claims"]))

        b3 = by_id["branch_3_2026_arc_mbnl1_rna_ev_preprint"]
        self.assertEqual(b3["epistemic_status"], ":Incomplete")
        self.assertIn("PMC13061017", b3["primary_citations"][0])
        self.assertTrue(any("NOT directly encapsulated" in c for c in b3["supported_claims"]))

        b4 = by_id["branch_4_candidate_mutational_hotspots"]
        self.assertEqual(b4["epistemic_status"], ":Incomplete")
        hotspots = b4["pre_registered_hotspot_hypotheses"]
        self.assertEqual(len(hotspots), 4)
        for h in hotspots:
            self.assertEqual(h["epistemic_status"], ":Incomplete")
            self.assertIn("Falsified if", h["pre_stated_falsification_readout"])
            self.assertIn("perturbation_category", h)
            self.assertIn("required_experimental_controls", h)
        self.assertNotIn("homotypic", hotspots[0]["structural_observation"].lower())
        self.assertIn("inter-subunit homologous-position aromatic lock", hotspots[0]["structural_observation"])

        h4 = hotspots[3]
        self.assertEqual(h4["hotspot_id"], "hotspot_4_CYS227_CYS230_HIS235_CYS240_nc_zn_knuckle")
        self.assertEqual(
            h4["perturbation_category"],
            "METAL_COORDINATION_AND_NC_FOLDING_PERTURBATION_NOT_ORDINARY_CA_INTERFACE_HOTSPOT",
        )
        self.assertIn("Zn occupancy (ICP-MS/PAR), NC folding (CD/NMR), and RNA-binding controls", h4["required_experimental_controls"])
        self.assertIn("12 unordered residue pairs", arrangements[0]["physical_mechanism_and_location"])
        self.assertIn("14.8186 A", arrangements[0]["physical_mechanism_and_location"])
        self.assertIn("163.70 A", arrangements[0]["physical_mechanism_and_location"])
        self.assertIn("_atom_site.auth_seq_id = 301", arrangements[0]["mmcif_entity_and_chain_provenance"])

    def test_07_point_symmetry_I_and_mmcif_entity_reconciliation(self) -> None:
        for receipt in (self.l2_receipt, self.l3_receipt):
            prov = receipt["mmcif_provenance"]
            self.assertEqual(prov["point_symmetry"], "I")
            self.assertIn("Cubic space-group notation I432 (No. 211, octahedral) is not applicable to 6TAP.", prov["symmetry_provenance_note"])
            inv = {r["label_asym_id"]: r for r in prov["struct_asym_inventory"]}
            self.assertEqual(set(inv.keys()), {"A", "B", "C", "D", "E", "F"})
            for ch in ("A", "B", "C", "D", "E"):
                self.assertEqual(inv[ch]["auth_asym_id"], ch)
                self.assertEqual(inv[ch]["label_entity_id"], "1")
                self.assertEqual(inv[ch]["struct_asym_entity_id"], "1")
            self.assertEqual(inv["E"]["functional_component_id"], "NC_ZINC_KNUCKLE_FRAGMENT")
            self.assertEqual(
                inv["E"]["biological_domain_or_fragment_designation"],
                "NC_zinc_knuckle_C_terminal_fragment_224_250_of_entity_1",
            )
            self.assertEqual(inv["F"]["auth_asym_id"], "E")
            self.assertEqual(inv["F"]["label_entity_id"], "2")
            self.assertEqual(inv["F"]["struct_asym_entity_id"], "2")
            self.assertEqual(inv["F"]["functional_component_id"], "ZN_COORDINATION_LIGAND")
            self.assertEqual(
                inv["F"]["biological_domain_or_fragment_designation"],
                "ZN_tetrahedral_CCHC_ligand_entity_2",
            )

    def test_08_alternative_similarity_metrics_and_threshold_sensitivity(self) -> None:
        alt = self.l3_receipt["alternative_similarity_metrics_and_threshold_sensitivity"]
        self.assertEqual(
            alt["threshold_sweep_grid"],
            [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70],
        )
        seq_metrics = alt["homologous_seq_id_metrics"]
        pair_metrics = alt["canonical_seq_pair_metrics"]
        for mname in (
            "jaccard",
            "sorensen_dice",
            "szymkiewicz_simpson_overlap",
            "cosine_ochiai",
            "contact_degree_weighted_cosine",
        ):
            self.assertIn(mname, seq_metrics)
            mat = seq_metrics[mname]["similarity_matrix"]
            self.assertEqual(len(mat), 10)
            for u in CANONICAL_SUBINTERFACE_ORDER:
                self.assertAlmostEqual(mat[u][u], 1.0, places=4)
                for v in CANONICAL_SUBINTERFACE_ORDER:
                    self.assertAlmostEqual(mat[u][v], mat[v][u], places=4)
            sl = seq_metrics[mname]["single_linkage_threshold_sweep"]
            cl = seq_metrics[mname]["complete_linkage_threshold_sweep"]
            self.assertEqual(len(sl["sweep_by_tau"]), 10)
            self.assertEqual(len(cl["sweep_by_tau"]), 10)
            self.assertTrue(len(sl["family_I_invariant_grid_interval"]) == 2)
            self.assertTrue(len(sl["family_II_invariant_grid_interval"]) == 2)
            self.assertIn("family_III_A_B_op2_complete_linkage_detachment_threshold", cl["exact_continuous_thresholds"])
            self.assertIn("family_III_A_B_op2_single_linkage_detachment_threshold", sl["exact_continuous_thresholds"])

        for mname in (
            "jaccard",
            "sorensen_dice",
            "szymkiewicz_simpson_overlap",
            "cosine_ochiai",
        ):
            self.assertIn(mname, pair_metrics)
            mat = pair_metrics[mname]["similarity_matrix"]
            self.assertEqual(len(mat), 10)

        # Verify Jaccard homologous seq_id exact interval and Family III splitting order
        j_sl = seq_metrics["jaccard"]["single_linkage_threshold_sweep"]
        j_cl = seq_metrics["jaccard"]["complete_linkage_threshold_sweep"]
        self.assertEqual(j_sl["family_I_invariant_grid_interval"], [0.15, 0.50])
        self.assertEqual(j_cl["family_I_invariant_grid_interval"], [0.15, 0.45])
        self.assertEqual(j_sl["family_II_invariant_grid_interval"], [0.15, 0.70])
        self.assertEqual(j_cl["family_II_invariant_grid_interval"], [0.15, 0.70])
        self.assertAlmostEqual(
            j_cl["exact_continuous_thresholds"]["family_III_A_B_op2_complete_linkage_detachment_threshold"],
            0.1429,
            places=4,
        )
        self.assertAlmostEqual(
            j_sl["exact_continuous_thresholds"]["family_III_A_B_op2_single_linkage_detachment_threshold"],
            0.2222,
            places=4,
        )
        self.assertIn(
            "181, 204",
            alt["family_invariance_and_splitting_summary"]["structural_sensitivity_conclusion"],
        )
        self.assertEqual(
            alt["family_robustness_stratification"],
            self.l3_receipt["symmetry_normalized_residue_fingerprints"]["family_robustness_stratification"],
        )
        # At tau = 0.25 under single-linkage Jaccard, A-B_op2 has detached while {C-C_op33, B-D_op33, A-D_op36} remains connected:
        self.assertTrue(j_sl["sweep_by_tau"]["0.25"]["A_B_op2_detached_as_singleton"])
        self.assertTrue(j_sl["sweep_by_tau"]["0.25"]["family_III_trimer_plus_A_D_op36_cluster_present"])
        # At tau = 0.40 under single-linkage Jaccard, A-D_op36 has also detached, leaving the C3 core {C-C_op33, B-D_op33}:
        self.assertTrue(j_sl["sweep_by_tau"]["0.40"]["A_D_op36_detached_as_singleton"])
        self.assertTrue(j_sl["sweep_by_tau"]["0.40"]["family_III_C3_core_cluster_present"])

    def test_09_reproducibility_and_zero_dashes(self) -> None:
        recomputed = build_level3_interface_comparison_receipt(
            l2_receipt_path=self.l2_receipt_path,
            crosscheck_receipt_path=self.cc_receipt_path,
        )
        self.assertEqual(recomputed, self.l3_receipt)

        l3_spec_text = self.l3_spec_path.read_text(encoding="utf-8")
        self.assertIn("`83` unique sequence positions", l3_spec_text)
        self.assertIn("`35` positions", l3_spec_text)
        self.assertIn("`9,400.64 A^2`", l3_spec_text)
        self.assertIn("`287,106.60`", l3_spec_text)
        self.assertIn("interface-fingerprint families", l3_spec_text)
        self.assertIn("EMD-10423", l3_spec_text)
        self.assertNotIn("EMD-10419", l3_spec_text)
        self.assertIn("`0.1429` / `0.2222`", l3_spec_text)
        self.assertIn("Family_III_Residual_Distal_CTD_Hierarchy", l3_spec_text)
        self.assertIn("family_robustness_stratification", l3_spec_text)
        self.assertIn("MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY", l3_spec_text)
        self.assertIn("degree_and_interface_preserving_permutation_null", l3_spec_text)
        self.assertIn("julia_assertion_count_reconciliation_3059_to_2996", l3_spec_text)
        self.assertIn("julia_assertion_count_restoration_2996_to_3059", l3_spec_text)
        self.assertIn("state_space_specification", l3_spec_text)
        self.assertIn(
            "CONSERVATION_INVARIANT_AND_IMPLEMENTATION_CHECK_NOT_EMPIRICAL_SIGNIFICANCE_STATISTIC",
            l3_spec_text,
        )
        self.assertIn("multiplicity_adjustment_analysis", l3_spec_text)
        self.assertIn(
            "| `A-A_C5_pentamer_op36` | I | `41` | `20 / 21` | `40 / 1` | `60 / 22` | `20 / 20 / 1` (`21`) | `820` | `695` (`84.76%`) | `20` (`2.44%` / `2.88%`) | `[19, 21]` | No |",
            l3_spec_text,
        )
        self.assertIn(
            "`48 / 2,865 = 1.68%` of all unordered within-subinterface pairs; `48 / 2,377 = 2.02%` of valid within-subinterface swaps",
            l3_spec_text,
        )
        self.assertIn("stale_table_reconciliation_note", l3_spec_text)
        self.assertIn(
            "| **Total Across 10 Sub-Interfaces** | **All** | **`204`** | **`96 / 108`** | **`89 / 115`** | **`185 / 223`** |",
            l3_spec_text,
        )
        self.assertIn(
            "remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run",
            l3_spec_text,
        )

        for path in (
            self.l1_followup_path,
            self.l2_receipt_path,
            self.cc_receipt_path,
            self.l3_receipt_path,
            self.l2_spec_path,
            self.l3_spec_path,
            REPO_ROOT / "src" / "m7_arc_level2_independent_crosscheck.py",
            REPO_ROOT / "src" / "m7_arc_level3_interface_comparison.py",
            REPO_ROOT / "code" / "test_m7_arc_level3_interface_comparison.py",
        ):
            txt = path.read_text(encoding="utf-8")
            self.assertNotIn(chr(0x2014), txt, f"Em-dash found in {path}")
            self.assertNotIn(chr(0x2013), txt, f"En-dash found in {path}")


if __name__ == "__main__":
    unittest.main()

