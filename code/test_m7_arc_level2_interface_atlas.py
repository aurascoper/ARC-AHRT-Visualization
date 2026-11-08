#!/usr/bin/env python3
"""Automated unit and integration tests for M7 / 6TAP Level-2 Residue Interface Atlas."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

for _cand in (
    "/opt/homebrew/opt/python@3.14/bin/python3.14",
    "/Library/Frameworks/Python.framework/Versions/3.13/bin/python3",
):
    if os.path.exists(_cand) and os.path.realpath(sys.executable) != os.path.realpath(_cand):
        try:
            import jsonschema  # noqa: F401
            import numpy  # noqa: F401
            break
        except ImportError:
            os.execv(_cand, [_cand] + sys.argv)

import jsonschema
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import m7_arc_level2_interface_atlas as atlas  # noqa: E402


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class TestM7ArcLevel2InterfaceAtlas(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.cif_path = REPO_ROOT / "structures" / "6TAP.cif"
        cls.ops_path = REPO_ROOT / "structures" / "6TAP_operators.json"
        cls.schema_path = REPO_ROOT / "schemas" / "m7_arc_level2_residue_interface_graph.schema.json"
        cls.spec_path = REPO_ROOT / "docs" / "M7_ARC_LEVEL2_RESIDUE_SPEC.md"
        cls.script_path = REPO_ROOT / "src" / "m7_arc_level2_interface_atlas.py"
        cls.receipt_path = REPO_ROOT / "artifacts" / "m7_arc_level2_validation_receipt.json"

        with cls.schema_path.open("r", encoding="utf-8") as f:
            cls.schema = json.load(f)
        with cls.receipt_path.open("r", encoding="utf-8") as f:
            cls.receipt = json.load(f)

    def test_01_json_schema_validation_and_negative_rejection(self) -> None:
        validator = jsonschema.Draft202012Validator(self.schema)
        errors = list(validator.iter_errors(self.receipt))
        self.assertEqual(errors, [], f"Receipt failed schema validation: {errors}")

        # Negative test 1: non-zero D-D self-contact count must be rejected by schema const=0
        bad_dd = copy.deepcopy(self.receipt)
        bad_dd["work_packages"]["L2_A_residue_pair_contact_census"]["census_by_cutoff"]["5.0"][
            "unique_asym_residue_pair_count"
        ]["D-D"] = 1
        self.assertFalse(validator.is_valid(bad_dd))

        # Negative test 2: non-zero B-B self-contact count must be rejected by schema const=0
        bad_bb = copy.deepcopy(self.receipt)
        bad_bb["work_packages"]["L2_A_residue_pair_contact_census"]["census_by_cutoff"]["5.0"][
            "unique_asym_residue_pair_count"
        ]["B-B"] = 1
        self.assertFalse(validator.is_valid(bad_bb))

        # Negative test 3: invalid functional_component_id must be rejected
        bad_ent = copy.deepcopy(self.receipt)
        bad_ent["mmcif_provenance"]["functional_components"][0]["functional_component_id"] = "INVALID_COMPONENT"
        self.assertFalse(validator.is_valid(bad_ent))

    def test_02_sha256_provenance_pins(self) -> None:
        hashes = self.receipt["artifact_hashes"]
        prov = self.receipt["mmcif_provenance"]
        self.assertEqual(prov["point_symmetry"], "I")
        self.assertNotIn("I432", prov["point_symmetry"])
        self.assertIn("icosahedral point-group symmetry I", prov["symmetry_provenance_note"])
        self.assertEqual(
            prov["cif_sha256"],
            "61397d5502c44f6fd83d887de353a698ed6f70fbd00604e2fc1c85fd1cd68e69",
        )
        self.assertEqual(hashes["cif_sha256"], sha256_file(self.cif_path))
        self.assertEqual(hashes["operators_sha256"], sha256_file(self.ops_path))
        self.assertEqual(hashes["schema_sha256"], sha256_file(self.schema_path))
        self.assertEqual(hashes["spec_sha256"], sha256_file(self.spec_path))
        self.assertEqual(hashes["generator_script_sha256"], sha256_file(self.script_path))

    def test_03_mmcif_entity_asym_operator_inventory_reconciliation(self) -> None:
        prov = self.receipt["mmcif_provenance"]
        self.assertEqual(prov["point_symmetry"], "I")
        self.assertEqual(len(prov["deposited_entities"]), 2)
        self.assertEqual(prov["deposited_entities"][0]["entity_id"], "1")
        self.assertEqual(prov["deposited_entities"][0]["entity_type"], "polymer")
        self.assertEqual(prov["deposited_entities"][1]["entity_id"], "2")
        self.assertEqual(prov["deposited_entities"][1]["entity_type"], "non-polymer")

        comp_ids = [row["functional_component_id"] for row in prov["functional_components"]]
        self.assertEqual(
            comp_ids,
            ["CA_CAPSID_SHELL", "NC_ZINC_KNUCKLE_FRAGMENT", "ZN_COORDINATION_LIGAND"],
        )
        comp_entities = [row["label_entity_id"] for row in prov["functional_components"]]
        self.assertEqual(comp_entities, ["1", "1", "2"])

        expected_asym_designations = {
            "A": ("A", "1", "1", "CA_CAPSID_SHELL", "CA_capsid_shell_domain_pentamer_A_41_205"),
            "B": ("B", "1", "1", "CA_CAPSID_SHELL", "CA_capsid_shell_domain_hexamer_B_43_205"),
            "C": ("C", "1", "1", "CA_CAPSID_SHELL", "CA_capsid_shell_domain_hexamer_C_43_206"),
            "D": ("D", "1", "1", "CA_CAPSID_SHELL", "CA_capsid_shell_domain_hexamer_D_41_206"),
            "E": ("E", "1", "1", "NC_ZINC_KNUCKLE_FRAGMENT", "NC_zinc_knuckle_C_terminal_fragment_224_250_of_entity_1"),
            "F": ("E", "2", "2", "ZN_COORDINATION_LIGAND", "ZN_tetrahedral_CCHC_ligand_entity_2"),
        }
        asym_by_id = {row["label_asym_id"]: row for row in prov["struct_asym_inventory"]}
        self.assertEqual(set(asym_by_id.keys()), {"A", "B", "C", "D", "E", "F"})
        for asym_id, (auth_id, lbl_ent, st_ent, fcomp, desig) in expected_asym_designations.items():
            row = asym_by_id[asym_id]
            self.assertEqual(row["auth_asym_id"], auth_id)
            self.assertEqual(row["label_entity_id"], lbl_ent)
            self.assertEqual(row["struct_asym_entity_id"], st_ent)
            self.assertEqual(row["functional_component_id"], fcomp)
            self.assertEqual(row["biological_domain_or_fragment_designation"], desig)

        inv = prov["inventory_reconciliation"]
        self.assertEqual(inv["deposited_ca_chains"], 4)
        self.assertEqual(inv["operators_count"], 60)
        self.assertEqual(inv["expanded_ca_instances"], 240)
        self.assertEqual(inv["deposited_nc_chains"], 1)
        self.assertEqual(inv["expanded_nc_instances"], 60)
        self.assertEqual(inv["deposited_zn_ions"], 1)
        self.assertEqual(inv["expanded_zn_instances"], 60)
        self.assertEqual(inv["deposited_water_or_other_ligands"], 0)
        self.assertEqual(inv["expanded_total_chain_instances"], 360)

        node_schema = self.receipt["node_schema"]
        self.assertEqual(node_schema["asu_residue_node_count_ca"], 659)
        self.assertEqual(node_schema["asu_residue_node_count_total"], 687)
        self.assertEqual(node_schema["expanded_residue_node_count_ca"], 39540)
        self.assertEqual(node_schema["expanded_residue_node_count_total"], 41220)

    def test_04_node_and_edge_schema_plus_seq_id_preservation(self) -> None:
        node_schema = self.receipt["node_schema"]
        self.assertEqual(
            node_schema["canonical_node_identity_tuple"],
            [
                "structure_id",
                "assembly_id",
                "entity_id",
                "asym_id",
                "operator_id",
                "label_seq_id",
                "auth_seq_id",
                "insertion_code",
            ],
        )
        edge_schema = self.receipt["edge_schema"]
        self.assertEqual(
            edge_schema["canonical_edge_identity_tuple"],
            [
                "residue_instance_i",
                "residue_instance_j",
                "cutoff_definition_id",
                "contact_metric_id",
            ],
        )
        self.assertEqual(
            edge_schema["required_edge_attributes"],
            [
                "minimum_heavy_atom_distance",
                "atom_pair_count",
                "residue_pair_contact",
                "buried_surface_area",
                "interaction_features",
                "interface_class",
                "operator_pair",
                "symmetry_axis_relation",
            ],
        )

        # Verify deposited 6TAP ZN node preserves label_asym_id="F", auth_asym_id="E",
        # label_seq_id=None, auth_seq_id=301, insertion_code="."
        asym_f = [
            row for row in self.receipt["mmcif_provenance"]["struct_asym_inventory"]
            if row["label_asym_id"] == "F"
        ][0]
        self.assertEqual(asym_f["auth_asym_id"], "E")
        self.assertEqual(asym_f["label_entity_id"], "2")
        self.assertEqual(asym_f["struct_asym_entity_id"], "2")
        self.assertEqual(asym_f["mmcif_label_entity_id"], "2")
        self.assertEqual(asym_f["functional_component_id"], "ZN_COORDINATION_LIGAND")
        self.assertEqual(
            asym_f["biological_domain_or_fragment_designation"],
            "ZN_tetrahedral_CCHC_ligand_entity_2",
        )
        self.assertEqual(asym_f["label_seq_id_range"], [None, None])
        self.assertEqual(asym_f["auth_seq_id_range"], [301, 301])

        sample_f = [n for n in node_schema["sample_nodes"] if n["asym_id"] == "F"][0]
        self.assertEqual(sample_f["node_id"], "6TAP.1.2.F.1...301..")
        self.assertIsNone(sample_f["label_seq_id"])
        self.assertEqual(sample_f["auth_seq_id"], 301)

    def test_05_synthetic_mmcif_label_vs_auth_seq_id_and_insertion_code(self) -> None:
        # Synthetic mmCIF where two distinct residues share auth_seq_id=100
        # ("100." and "100A") while having distinct label_seq_id (1 and 2).
        synthetic_cif = """data_TEST
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'Test peptide'
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1)' A
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'point symmetry operation' 1 0 0 0 0 1 0 0 0 0 1 0
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0.0 0.0 0.0 1.0 20.0 100 ALA X CA 1
ATOM 2 C CA . GLY A 1 2 A 3.8 0.0 0.0 1.0 20.0 100 GLY X CA 1
ATOM 3 C CA . VAL A 1 3 ? 7.6 0.0 0.0 1.0 20.0 101 VAL X CA 1
"""
        with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False, encoding="utf-8") as tmp:
            tmp.write(synthetic_cif)
            tmp_path = Path(tmp.name)
        try:
            parsed = atlas.parse_mmcif_file(tmp_path)
            residues = atlas.extract_residue_nodes_from_cif(parsed)["A"]
            self.assertEqual(len(residues), 3)
            self.assertEqual(residues[0]["label_seq_id"], 1)
            self.assertEqual(residues[0]["auth_seq_id"], 100)
            self.assertEqual(residues[0]["insertion_code"], ".")
            self.assertEqual(residues[0]["auth_asym_id"], "X")

            self.assertEqual(residues[1]["label_seq_id"], 2)
            self.assertEqual(residues[1]["auth_seq_id"], 100)
            self.assertEqual(residues[1]["insertion_code"], "A")
            self.assertNotEqual(residues[0]["node_id"], residues[1]["node_id"])
            self.assertEqual(residues[0]["node_id"], "TEST.1.1.A.1.1.100..")
            self.assertEqual(residues[1]["node_id"], "TEST.1.1.A.1.2.100.A")
        finally:
            tmp_path.unlink(missing_ok=True)

        # Fail-closed test on empty or malformed mmCIF
        with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False, encoding="utf-8") as bad_tmp:
            bad_tmp.write("")
            bad_path = Path(bad_tmp.name)
        try:
            with self.assertRaises(ValueError):
                atlas.parse_mmcif_file(bad_path)
        finally:
            bad_path.unlink(missing_ok=True)

    def test_06_l2a_contact_metrics_and_self_contact_invariants(self) -> None:
        l2a = self.receipt["work_packages"]["L2_A_residue_pair_contact_census"]
        cutoffs = ["3.5", "4.0", "4.5", "5.0", "5.5", "6.0"]
        self.assertEqual([f"{c:.1f}" for c in l2a["cutoffs_angstrom"]], cutoffs)

        prev_unique = 0
        prev_inst = 0
        prev_hits = 0
        for c_key in cutoffs:
            entry = l2a["census_by_cutoff"][c_key]
            u_tot = entry["totals"]["unique_asym_residue_pair_count_total"]
            i_tot = entry["totals"]["instance_residue_pair_observation_count_total"]
            h_tot = entry["totals"]["heavy_atom_pair_hit_count_total"]
            self.assertEqual(i_tot, 60 * u_tot)
            self.assertGreaterEqual(h_tot, i_tot)
            self.assertGreaterEqual(u_tot, prev_unique)
            self.assertGreaterEqual(i_tot, prev_inst)
            self.assertGreaterEqual(h_tot, prev_hits)
            prev_unique, prev_inst, prev_hits = u_tot, i_tot, h_tot

            # B-B, D-D, and A-C must be strictly 0 at ALL cutoffs
            for zero_pair in ("B-B", "D-D", "A-C"):
                self.assertEqual(entry["unique_asym_residue_pair_count"][zero_pair], 0)
                self.assertEqual(entry["instance_residue_pair_observation_count"][zero_pair], 0)
                self.assertEqual(entry["heavy_atom_pair_hit_count"][zero_pair], 0)

        # Exact 5.0 A verification
        c50 = l2a["census_by_cutoff"]["5.0"]
        self.assertEqual(c50["totals"]["unique_asym_residue_pair_count_total"], 204)
        self.assertEqual(c50["totals"]["instance_residue_pair_observation_count_total"], 12240)
        self.assertEqual(c50["totals"]["heavy_atom_pair_hit_count_total"], 66390)
        expected_50_unique = {
            "A-A": 41,
            "A-B": 1,
            "A-C": 0,
            "A-D": 32,
            "B-B": 0,
            "B-C": 54,
            "B-D": 35,
            "C-C": 7,
            "C-D": 34,
            "D-D": 0,
        }
        for cp, exp_u in expected_50_unique.items():
            self.assertEqual(c50["unique_asym_residue_pair_count"][cp], exp_u)
            self.assertEqual(c50["instance_residue_pair_observation_count"][cp], 60 * exp_u)

        self.assertEqual(len(l2a["complete_residue_pair_list_at_5A"]), 204)
        self.assertEqual(len(self.receipt["edge_schema"]["edges_at_5A"]), 204)

    def test_07_l2b_sasa_and_buried_surface_area(self) -> None:
        l2b = self.receipt["work_packages"]["L2_B_sasa_buried_area"]
        self.assertEqual(len(l2b["per_residue_delta_sasa"]), 659)
        for row in l2b["per_residue_delta_sasa"]:
            self.assertGreaterEqual(row["delta_sasa_A2"], 0.0)
            self.assertAlmostEqual(
                row["isolated_subunit_sasa_A2"] - row["assembly_sasa_A2"],
                row["delta_sasa_A2"],
                places=2,
            )

        self.assertAlmostEqual(l2b["representative_subinterface_pairwise_bsa_sum_A2"], 9570.22, places=2)
        self.assertAlmostEqual(l2b["representative_subinterface_per_side_bsa_sum_A2"], 4785.10, places=2)
        self.assertAlmostEqual(l2b["whole_capsid_total_buried_surface_area_A2"], 574213.20, places=2)
        self.assertAlmostEqual(l2b["whole_capsid_one_side_interface_area_A2"], 287106.60, places=2)

        by_class = {
            row["interface_class"]: row for row in l2b["interface_buried_surface_area_by_sub_interface"]
        }
        self.assertEqual(len(by_class), 10)
        for row in by_class.values():
            self.assertEqual(row["unordered_capsid_edge_multiplicity"], 60)
            self.assertAlmostEqual(
                row["whole_capsid_sub_interface_buried_surface_area_A2"],
                60.0 * row["total_buried_surface_area_A2"],
                places=2,
            )
            self.assertAlmostEqual(
                row["whole_capsid_one_side_bsa_A2"],
                60.0 * row["buried_surface_area_per_side_A2"],
                places=2,
            )
        self.assertAlmostEqual(by_class["A-A_C5_pentamer_op36"]["buried_surface_area_A2"], 1882.00, places=2)
        self.assertAlmostEqual(by_class["C-D_C1_intra_asu_op1"]["buried_surface_area_A2"], 1487.53, places=2)
        self.assertAlmostEqual(by_class["B-C_C3_inter_asu_hexamer_op33"]["buried_surface_area_A2"], 1365.72, places=2)
        self.assertAlmostEqual(by_class["A-B_C2_peripheral_op2"]["buried_surface_area_A2"], 48.65, places=2)

    def test_08_l2c_interface_chemistry(self) -> None:
        l2c = self.receipt["work_packages"]["L2_C_interface_chemistry"]
        self.assertEqual(len(l2c["hydrophobic_contacts_at_4p5A"]), 20)
        self.assertEqual(len(l2c["aromatic_pi_stacking_candidates"]), 10)
        self.assertEqual(len(l2c["hydrogen_bond_candidates"]), 19)
        self.assertEqual(len(l2c["salt_bridge_candidates"]), 9)
        self.assertIn("candidate_qualification_caveats", l2c)

        for entry in l2c["summary_by_sub_interface"]:
            self.assertEqual(entry["hydrophobic_candidate_pairs"], entry["hydrophobic_residue_pairs_4p5A"])
            self.assertEqual(entry["pi_stacking_candidate_pairs"], entry["aromatic_pi_stacking_candidates_7p0A"])
            self.assertEqual(
                entry["hbond_heavy_atom_distance_candidate_pairs"], entry["hydrogen_bond_candidates_3p5A"]
            )
            self.assertEqual(
                entry["salt_bridge_heavy_atom_distance_candidate_pairs"], entry["salt_bridge_candidates_4p0A"]
            )

        for hb in l2c["hydrogen_bond_candidates"]:
            self.assertLessEqual(hb["donor_acceptor_distance_A"], 3.5)
        for sb in l2c["salt_bridge_candidates"]:
            self.assertLessEqual(sb["nitrogen_oxygen_distance_A"], 4.0)
        for pi in l2c["aromatic_pi_stacking_candidates"]:
            self.assertLessEqual(pi["ring_centroid_distance_A"], 7.0)
            self.assertEqual(pi["ring_plane_angle_deg"], pi["ring_normal_angle_deg"])
            self.assertIn("min_heavy_atom_distance_A", pi)

        eq = self.receipt["epistemic_qualifications"]
        self.assertEqual(
            self.receipt["verification_disposition"],
            "remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run",
        )
        self.assertEqual(
            eq["verification_disposition"],
            "remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run",
        )
        prov = self.receipt["publication_and_ci_provenance"]
        self.assertEqual(
            prov["github_pr_url"],
            "https://github.com/aurascoper/flow-cytometry-model-bridge/pull/6",
        )
        self.assertEqual(
            prov["historical_v1_baseline_commit"],
            "955e32dc673d3a874f874ffb8965497db03ca112",
        )
        self.assertEqual(
            prov["level1_v2_repair_commit"],
            "6a09714d03fccf29905c99bddc416785781063d1",
        )
        self.assertEqual(
            prov["level2_atlas_and_crosscheck_commit"],
            "fc5e1a407e375014653b2b72f73b89559c14e2ee",
        )
        self.assertEqual(
            prov["level3_initial_publication_commit"],
            "1abeaad2e44fe25e3e4cd17b7920a51f990fc55e",
        )
        self.assertEqual(len(eq["m7_falsification_design_rules"]), 4)
        hist = self.receipt["historical_level1_superseded_evidence"]
        self.assertEqual(hist["superseded_homotypic_rationale_status"], ":Superseded")
        self.assertIn(
            "B and D both lack homotypic interfaces in the 6TAP residue-contact census",
            hist["corrected_homotypic_census_statement"],
        )

    def test_09_l2d_icosahedral_operator_axis_mapping(self) -> None:
        l2d = self.receipt["work_packages"]["L2_D_icosahedral_operator_axis_mapping"]
        self.assertEqual(
            l2d["rotation_order_spectrum"],
            {"C1_identity": 1, "C2": 15, "C3": 20, "C5": 24},
        )
        self.assertTrue(l2d["group_verification"]["all_so3_verified"])
        self.assertTrue(l2d["group_verification"]["group_closure_verified"])
        self.assertLess(l2d["group_verification"]["max_orthonormality_error"], 1e-5)

        self.assertEqual(l2d["unique_axis_line_counts"]["C5_axes"], 6)
        self.assertEqual(l2d["unique_axis_line_counts"]["C3_axes"], 10)
        self.assertEqual(l2d["unique_axis_line_counts"]["C2_axes"], 15)

        for op_row in l2d["all_60_operators_axis_table"]:
            self.assertEqual(op_row["nearest_canonical_axis_alignment_deg"], 0.0)

        for ch in ("A", "B", "C", "D", "E", "F"):
            prof = l2d["subunit_radial_and_axis_alignment"][ch]
            self.assertEqual(prof["symmetry_copy_angle_std_deg"], 0.0)

        # Cross-deposit proof that op_id is deposit-specific
        comp = {
            (row["structure_id"], row["operator_id"]): row
            for row in l2d["deposit_specific_op_id_proof"]["cross_deposit_comparison"]
        }
        self.assertEqual(comp[("6TAP", "2")]["symmetry_symbol"], "C2")
        self.assertEqual(comp[("6TAP", "36")]["symmetry_symbol"], "C5")
        self.assertEqual(comp[("1OHF", "2")]["symmetry_symbol"], "C5")
        self.assertEqual(comp[("6R24", "37")]["symmetry_symbol"], "C5")

    def test_10_l2e_nc_zn_coordination_and_contacts(self) -> None:
        l2e = self.receipt["work_packages"]["L2_E_nc_zn_spatial_context"]
        zn = l2e["zn_coordination_sphere"]
        self.assertEqual(zn["coordination_geometry"], "tetrahedral_CCHC_zinc_knuckle")
        self.assertEqual(len(zn["coordinating_ligands"]), 4)
        for lig in zn["coordinating_ligands"]:
            self.assertGreaterEqual(lig["computed_distance_A"], 1.95)
            self.assertLessEqual(lig["computed_distance_A"], 2.45)
        self.assertAlmostEqual(zn["mean_bond_angle_deg"], 109.42, places=1)
        self.assertEqual(len(l2e["nc_ca_proximity_by_cutoff"]["5.0"]["contacts"]), 16)

    def test_11_fallback_kdtree_equivalence_and_zero_dashes(self) -> None:
        rng = np.random.default_rng(42)
        pts_a = rng.uniform(0.0, 20.0, size=(80, 3))
        pts_b = rng.uniform(0.0, 20.0, size=(75, 3))
        fb_a = atlas._FallbackKDTree(pts_a)
        fb_b = atlas._FallbackKDTree(pts_b)
        res_fb = fb_a.query_ball_tree(fb_b, 4.5)
        for idx_i, neighbors in enumerate(res_fb):
            dists = np.linalg.norm(pts_b - pts_a[idx_i], axis=1)
            expected = sorted(np.where(dists <= 4.5)[0].tolist())
            self.assertEqual(sorted(neighbors), expected)

        # Zero em-dashes and zero en-dashes across all Level-2 files
        for p in (self.cif_path, self.schema_path, self.spec_path, self.script_path, self.receipt_path):
            text = p.read_text(encoding="utf-8")
            self.assertEqual(text.count("\u2014"), 0, f"Em-dash found in {p}")
            self.assertEqual(text.count("\u2013"), 0, f"En-dash found in {p}")

    def test_12_parser_edge_cases_multimodel_altloc_oper_expression_so3_and_spec_sync(self) -> None:
        # 1. Loop with embedded # comments and blank lines, multi-model NMR, altloc A/B, and H/D atoms
        cif_with_comments_altloc_models = """data_EDGE
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA A ALA A 1 1 ? 1.0 0.0 0.0 0.35 20.0 10 ALA A CA 1
# Comment inside loop_ body must not terminate loop_

ATOM 2 C CA B ALA A 1 1 ? 2.0 0.0 0.0 0.65 18.0 10 ALA A CA 1
ATOM 3 H HA . ALA A 1 1 ? 1.5 0.5 0.0 1.00 20.0 10 ALA A HA 1
ATOM 4 D DA . ALA A 1 1 ? 1.5 -0.5 0.0 1.00 20.0 10 ALA A DA 1
ATOM 5 C CA . ALA A 1 1 ? 9.0 9.0 9.0 1.00 25.0 10 ALA A CA 2
"""
        parsed = atlas.parse_mmcif_text(cif_with_comments_altloc_models)
        raw_rows = atlas.extract_loop_dicts(parsed, "_atom_site")
        self.assertEqual(len(raw_rows), 5)
        filt_m1 = atlas.filter_canonical_heavy_atom_rows(raw_rows, model_num="1", heavy_atoms_only=True)
        self.assertEqual(len(filt_m1), 1)
        self.assertEqual(filt_m1[0]["id"], "2")
        self.assertEqual(filt_m1[0]["label_alt_id"], "B")
        residues_m1 = atlas.extract_residue_nodes_from_cif(parsed, model_num="1", heavy_atoms_only=True)["A"]
        self.assertEqual(len(residues_m1), 1)
        self.assertEqual(residues_m1[0]["atom_count"], 1)

        filt_m2 = atlas.filter_canonical_heavy_atom_rows(raw_rows, model_num="2", heavy_atoms_only=True)
        self.assertEqual(len(filt_m2), 1)
        self.assertEqual(filt_m2[0]["id"], "5")
        residues_m2 = atlas.extract_residue_nodes_from_cif(parsed, model_num="2", heavy_atoms_only=True)["A"]
        self.assertEqual(len(residues_m2), 1)
        self.assertEqual(residues_m2[0]["atom_count"], 1)

        with self.assertRaises(ValueError):
            atlas.extract_residue_nodes_from_cif(parsed, model_num="99")

        # 2. Fail-closed on orphan token outside loop and multi-token single tag
        with self.assertRaises(ValueError):
            atlas.parse_mmcif_text("data_BAD\n_entry.id 6TAP EXTRA_TOKEN\n")
        with self.assertRaises(ValueError):
            atlas.parse_mmcif_text("data_BAD\n_entry.id 6TAP\nORPHAN_LINE_OUTSIDE_LOOP\n")

        # 3. Compound oper_expression expansion: single list/range and Cartesian product (1-2)(3-4)
        oper_rows = [
            {
                "id": "1",
                "matrix[1][1]": "1", "matrix[1][2]": "0", "matrix[1][3]": "0", "vector[1]": "0",
                "matrix[2][1]": "0", "matrix[2][2]": "1", "matrix[2][3]": "0", "vector[2]": "0",
                "matrix[3][1]": "0", "matrix[3][2]": "0", "matrix[3][3]": "1", "vector[3]": "0",
            },
            {
                "id": "2",
                "matrix[1][1]": "-1", "matrix[1][2]": "0", "matrix[1][3]": "0", "vector[1]": "10",
                "matrix[2][1]": "0", "matrix[2][2]": "-1", "matrix[2][3]": "0", "vector[2]": "0",
                "matrix[3][1]": "0", "matrix[3][2]": "0", "matrix[3][3]": "1", "vector[3]": "0",
            },
            {
                "id": "3",
                "matrix[1][1]": "1", "matrix[1][2]": "0", "matrix[1][3]": "0", "vector[1]": "0",
                "matrix[2][1]": "0", "matrix[2][2]": "-1", "matrix[2][3]": "0", "vector[2]": "0",
                "matrix[3][1]": "0", "matrix[3][2]": "0", "matrix[3][3]": "-1", "vector[3]": "5",
            },
        ]
        ops_prod = atlas.parse_and_expand_oper_expression("(1,2)(1,3)", oper_rows)
        self.assertEqual([op["id"] for op in ops_prod], ["1x1", "1x3", "2x1", "2x3"])
        np.testing.assert_allclose(ops_prod[3]["R"], np.diag([-1.0, 1.0, -1.0]))
        np.testing.assert_allclose(ops_prod[3]["t"], np.array([10.0, 0.0, 5.0]))
        with self.assertRaises(ValueError):
            atlas.parse_and_expand_oper_expression("(1-999)", oper_rows)

        # 4. SO(3) validation in derive_operator_axis_and_angle rejects reflection and non-orthogonal matrices
        with self.assertRaises(ValueError):
            atlas.derive_operator_axis_and_angle(np.diag([1.0, 1.0, -1.0]))
        with self.assertRaises(ValueError):
            atlas.derive_operator_axis_and_angle(np.diag([2.0, 1.0, 1.0]))

        # 5. Spec-to-receipt numerical synchronization check
        spec_text = self.spec_path.read_text(encoding="utf-8")
        self.assertIn("**66,390**", spec_text)
        self.assertIn("**31,860**", spec_text)
        self.assertIn("**120,930**", spec_text)
        self.assertIn("**202,020**", spec_text)
        self.assertIn("TTRCTYCSFRGHTFDNCRKRQKDRQEE", spec_text)
        self.assertNotIn("NTRCTYCSFRGHTFDNCRKRQKDRQEE", spec_text)
        self.assertIn("15.4794 A", spec_text)
        self.assertIn("14.8186 A", spec_text)
        self.assertIn("163.6982", spec_text)
        self.assertIn("164.9754", spec_text)

    def test_13_independent_crosscheck_scope_and_voronoi_bsa_ols_calibration(self) -> None:
        crosscheck_path = REPO_ROOT / "artifacts" / "m7_arc_level2_independent_crosscheck_receipt.json"
        cc = json.loads(crosscheck_path.read_text(encoding="utf-8"))
        self.assertEqual(
            cc["independence_scope"],
            "COMPUTATIONAL_IMPLEMENTATION_INDEPENDENT_SAME_DEPOSITED_6TAP_COORDINATES",
        )
        self.assertEqual(cc["point_symmetry"], "I")
        self.assertEqual(
            cc["verification_disposition"],
            "remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run",
        )
        self.assertEqual(
            cc["publication_and_ci_provenance"]["github_pr_url"],
            "https://github.com/aurascoper/flow-cytometry-model-bridge/pull/6",
        )
        self.assertIn("same deposited PDB 6TAP", cc["independence_scope_statement"])

        cal = cc["sasa_and_voronoi_surface_crosscheck"][
            "voronoi_vs_bsa_geometric_distinction_and_linear_calibration"
        ]
        ols_2s = cal["ols_vs_two_sided_shrake_rupley_pairwise_bsa"]
        ols_1s = cal["ols_vs_one_side_shrake_rupley_bsa"]
        self.assertAlmostEqual(ols_2s["slope"], 0.444876, places=4)
        self.assertAlmostEqual(ols_2s["intercept_A2"], -2.0833, places=2)
        self.assertGreater(ols_2s["r_squared"], 0.998)
        self.assertAlmostEqual(ols_2s["rmse_A2"], 11.1209, places=2)
        self.assertAlmostEqual(ols_2s["max_abs_residual_A2"], 24.1430, places=2)
        self.assertAlmostEqual(ols_1s["slope"], 0.889755, places=4)
        self.assertGreater(ols_1s["r_squared"], 0.998)

        subs = cc["sasa_and_voronoi_surface_crosscheck"]["sub_interface_crosschecks"]
        self.assertEqual(len(subs), 10)
        for s in subs:
            self.assertIn("predicted_voronoi_A2", s)
            self.assertIn("residual_A2", s)
            self.assertIn("voronoi_to_two_sided_bsa_ratio", s)
            self.assertIn("voronoi_to_one_side_bsa_ratio", s)
            self.assertAlmostEqual(
                s["voronoi_radical_plane_contact_area_A2"] - s["predicted_voronoi_A2"],
                s["residual_A2"],
                places=1,
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)


