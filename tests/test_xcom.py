import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from ahrt_stage_a.xcom import XCOMDataError, XCOMTable, coefficient_closure_error, load_and_verify_manifest
from ahrt_stage_a.validation import validate_xcom

FIXTURE = Path(__file__).parent / "fixtures" / "xcom_fixture_not_physical.csv"


class XCOMTests(unittest.TestCase):
    def test_fixture_load_and_interpolate(self):
        table = XCOMTable.from_csv(FIXTURE)
        coeff = table.at(0.0031622776601683794)
        self.assertGreater(coeff.total_without_coherent, 0)
        e0, e1 = coefficient_closure_error(coeff)
        self.assertLess(e0, 1e-12)
        self.assertLess(e1, 1e-12)


    def test_xcom_attenuation_gate_on_synthetic_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            table_copy = tmp / "fixture.csv"
            table_copy.write_bytes(FIXTURE.read_bytes())
            digest = hashlib.sha256(table_copy.read_bytes()).hexdigest()
            manifest = {
                "source": "https://physics.nist.gov/PhysRefData/Xcom/",
                "material": "TEST FIXTURE ONLY, NOT WATER",
                "units": "cm2/g; energy_MeV",
                "retrieved_utc": "2026-11-03T00:00:00Z",
                "sha256": digest,
            }
            manifest_path = tmp / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            result = validate_xcom(table_copy, manifest_path, 0.005, histories=500_000, seed=9)
            self.assertEqual(result["status"], "PASS")

    def test_manifest_hash_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            table_copy = tmp / "fixture.csv"
            table_copy.write_bytes(FIXTURE.read_bytes())
            digest = hashlib.sha256(table_copy.read_bytes()).hexdigest()
            manifest = {
                "source": "https://physics.nist.gov/PhysRefData/Xcom/",
                "material": "TEST FIXTURE ONLY, NOT WATER",
                "units": "cm2/g; energy_MeV",
                "retrieved_utc": "2026-11-03T00:00:00Z",
                "sha256": digest,
            }
            manifest_path = tmp / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            self.assertEqual(load_and_verify_manifest(table_copy, manifest_path)["sha256"], digest)
            table_copy.write_text(table_copy.read_text() + "\n", encoding="utf-8")
            with self.assertRaises(XCOMDataError):
                load_and_verify_manifest(table_copy, manifest_path)


if __name__ == "__main__":
    unittest.main()
