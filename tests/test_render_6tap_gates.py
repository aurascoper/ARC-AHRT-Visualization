import copy
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "render_6tap_interface_figure.py"
RECEIPT = ROOT / "artifacts" / "m7_arc_level3_interface_comparison_receipt.json"

spec = importlib.util.spec_from_file_location("render_6tap", SCRIPT)
render = importlib.util.module_from_spec(spec)
spec.loader.exec_module(render)


class ReceiptGateTests(unittest.TestCase):
    # The renderer must exit 0 on the committed receipt and nonzero on a receipt whose
    # stub-identity section disagrees with its rows, so an unattended run cannot look green
    # while printing an error. sys.exit(str) is exit status 1; sys.exit(None) would be 0.

    def setUp(self):
        self.receipt = json.loads(RECEIPT.read_text())

    def test_real_receipt_passes(self):
        render.receipt_gates(render.interface_rows(self.receipt), self.receipt)

    def test_tampered_stub_value_aborts_nonzero(self):
        tampered = copy.deepcopy(self.receipt)
        stub = tampered["constrained_graph_state_space_connectivity_and_mcmc_audit"][
            "stub_difference_conservation_identity"]["empirical_6tap_5A_cutoff_values"]
        stub["conserved_E_CC_minus_E_NN"] = 35  # the retired value the gate exists to catch
        with self.assertRaises(SystemExit) as cm:
            render.receipt_gates(render.interface_rows(tampered), tampered)
        self.assertNotIn(cm.exception.code, (None, 0))
        self.assertIn("receipt disagrees with itself", str(cm.exception.code))


if __name__ == "__main__":
    unittest.main()
