import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "render_6tap_interface_figure.py"
RECEIPT = ROOT / "artifacts" / "m7_arc_level3_interface_comparison_receipt.json"

spec = importlib.util.spec_from_file_location("render_6tap", SCRIPT)
render = importlib.util.module_from_spec(spec)
spec.loader.exec_module(render)


class ReceiptGateTests(unittest.TestCase):
    # The renderer must exit 0 on the committed receipt and nonzero on a receipt that
    # disagrees with itself, so an unattended run cannot look green while printing an
    # error. sys.exit(str) is exit status 1; sys.exit(None) would be 0.

    def setUp(self):
        self.receipt = json.loads(RECEIPT.read_text())
        self.audit = "constrained_graph_state_space_connectivity_and_mcmc_audit"

    def assertGateAborts(self, receipt):
        with self.assertRaises(SystemExit) as cm:
            render.receipt_gates(render.interface_rows(receipt), receipt)
        self.assertNotIn(cm.exception.code, (None, 0))
        self.assertIn("receipt disagrees with itself", str(cm.exception.code))

    def test_real_receipt_passes(self):
        render.receipt_gates(render.interface_rows(self.receipt), self.receipt)

    def test_tampered_stub_value_aborts_nonzero(self):
        tampered = copy.deepcopy(self.receipt)
        stub = tampered[self.audit]["stub_difference_conservation_identity"]["empirical_6tap_5A_cutoff_values"]
        stub["conserved_E_CC_minus_E_NN"] = 35  # the retired value the gate exists to catch
        self.assertGateAborts(tampered)

    def test_scheme_observed_count_mismatch_aborts(self):
        tampered = copy.deepcopy(self.receipt)
        key = "within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc"
        tampered[self.audit][key]["observed_same_domain_contacts"] = 130
        self.assertGateAborts(tampered)

    def test_missing_fingerprint_entry_aborts_with_message(self):
        tampered = copy.deepcopy(self.receipt)
        tampered["symmetry_normalized_residue_fingerprints"]["subinterface_fingerprints"].pop()
        self.assertGateAborts(tampered)


class PovClassificationTests(unittest.TestCase):
    # Synthetic scene, coordinates in Angstrom about centre (10, 20, 30). Capsid spheres come
    # in two radius groups (4 + 2) that match no other class, so both must fold into capsid.
    CENTRE = (10.0, 20.0, 30.0)
    SPHERES = (
        [((s * 100, 0, 0), 3.06) for s in (1, -1)] + [((0, s * 100, 0), 3.06) for s in (1, -1)]
        + [((0, 0, s * 100), 2.8) for s in (1, -1)]
        + [((50, 0, 0), 1.9), ((-25, 43.30127, 0), 1.9), ((-25, -43.30127, 0), 1.9)]
        + [((0, 0, 0), 1.39)]
    )

    def write_pov(self, tmp):
        path = Path(tmp) / "scene.pov"
        lines = ["// synthetic"] + [
            f"sphere{{<{x + self.CENTRE[0]},{y + self.CENTRE[1]},{z + self.CENTRE[2]}>, {r}}}"
            for (x, y, z), r in self.SPHERES]
        path.write_text("\n".join(lines) + "\n")
        return path

    def test_classes_counts_and_radii(self):
        with tempfile.TemporaryDirectory() as tmp:
            classes = render.read_pov(self.write_pov(tmp), {"capsid": 6, "nc": 3, "zn": 1})
        self.assertEqual({k: len(v) for k, v in classes.items()}, {"capsid": 6, "nc": 3, "zn": 1})
        for name, radius_nm in (("capsid", 10.0), ("nc", 5.0), ("zn", 0.0)):
            for sphere in classes[name]:
                self.assertAlmostEqual(sphere[3], radius_nm, places=4)

    def test_inventory_mismatch_aborts_nonzero(self):
        # nc=5 matches no radius group, so all nine non-Zn spheres land in capsid against 4.
        # (nc=4 would not abort: the 4-sphere capsid group would be taken as NC and the totals
        # would still balance, the coincidence read_pov's docstring warns about.)
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit) as cm:
                render.read_pov(self.write_pov(tmp), {"capsid": 4, "nc": 5, "zn": 1})
        self.assertNotIn(cm.exception.code, (None, 0))
        self.assertIn("do not match receipt inventory", str(cm.exception.code))


if __name__ == "__main__":
    unittest.main()
