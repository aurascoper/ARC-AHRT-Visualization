"""Cross-check the MATLAB dashboard sidecar against the Level 3 receipt and the Python figure sidecar.

MATLAB is not needed: these tests read the committed JSON files only. They are skipped
until matlab/renders/matlab_arc_ahrt_3d_receipt.json has been regenerated with the v2
schema (run matlab/arc_ahrt_3d_visualization.m under MATLAB).
"""
import hashlib
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MATLAB_SIDECAR = ROOT / "matlab" / "renders" / "matlab_arc_ahrt_3d_receipt.json"
PY_SIDECAR = ROOT / "docs" / "6tap_interface_atlas_figure.receipt.json"
RECEIPT = ROOT / "artifacts" / "m7_arc_level3_interface_comparison_receipt.json"
SHELLS = ROOT / "geometry" / "shells.json"

matlab = json.loads(MATLAB_SIDECAR.read_text())
IS_V2 = matlab.get("schema") == "arc-ahrt-matlab-3d-visualization-receipt-v2"


@unittest.skipUnless(IS_V2, "MATLAB sidecar is still v1; re-run the MATLAB script to regenerate it")
class MatlabSidecarCrossCheck(unittest.TestCase):
    def setUp(self):
        self.receipt = json.loads(RECEIPT.read_text())
        self.py = json.loads(PY_SIDECAR.read_text())

    def test_same_level3_receipt_as_python_figure(self):
        level3 = hashlib.sha256(RECEIPT.read_bytes()).hexdigest()
        self.assertEqual(matlab["inputs_sha256"]["level3_receipt"], level3)
        self.assertIn(level3, json.dumps(self.py["inputs_sha256"]))

    def test_residue_counts_match_python_pov_classes(self):
        classes = self.py["pov_sphere_classes"]
        self.assertEqual(matlab["residues_in_assembly"],
                         {"ca": classes["capsid"]["count"], "nc": classes["nc"]["count"], "zn": classes["zn"]["count"]})

    def test_identities_match_receipt(self):
        stub = self.receipt["constrained_graph_state_space_connectivity_and_mcmc_audit"][
            "stub_difference_conservation_identity"]["empirical_6tap_5A_cutoff_values"]
        self.assertEqual(matlab["receipt_identities"], {
            "E_NN": stub["observed_E_NN"], "E_CC": stub["observed_E_CC"], "E_CN": stub["observed_E_CN"],
            "E_same": stub["observed_E_same"], "E_CC_minus_E_NN": stub["conserved_E_CC_minus_E_NN"]})

    def test_outer_radius_read_from_shells_json(self):
        shells = json.loads(SHELLS.read_text())
        self.assertAlmostEqual(matlab["outer_radius_nm_from_shells_json"],
                               shells["6TAP"]["descriptors"]["radial_profile"]["outer_radius_nm"], places=9)


if __name__ == "__main__":
    unittest.main()
