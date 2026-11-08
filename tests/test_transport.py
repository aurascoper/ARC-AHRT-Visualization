import tempfile
import unittest
from pathlib import Path

from ahrt_stage_a.transport import PhotonTransport, TransportConfig
from ahrt_stage_a.xcom import XCOMTable

FIXTURE = Path(__file__).parent / "fixtures" / "xcom_fixture_not_physical.csv"


class TransportTests(unittest.TestCase):
    def test_smoke_transport_fixture_only(self):
        # The fixture is intentionally synthetic and is used only to test code paths.
        table = XCOMTable.from_csv(FIXTURE)
        config = TransportConfig(
            photon_energy_MeV=0.005,
            macro_radius_cm=0.1,
            macro_length_cm=0.2,
            cell_radius_cm=0.001,
            cell_length_cm=0.002,
            histories=100,
            seed=41,
        )
        with tempfile.TemporaryDirectory() as tmp:
            summary = PhotonTransport(table, config).run(tmp)
            self.assertEqual(summary["tallies"]["histories"], 100)
            self.assertTrue((Path(tmp) / "compton_electron_sources.csv").is_file())
            self.assertTrue((Path(tmp) / "run_summary.json").is_file())
            energy = summary["energy_accounting"]
            accounted = (
                energy["escaped_photon_energy_MeV"]
                + energy["compton_electron_source_energy_MeV"]
                + energy["photoelectric_local_energy_MeV"]
                + energy["below_cutoff_unresolved_energy_MeV"]
                + energy["interaction_cap_unresolved_energy_MeV"]
            )
            self.assertAlmostEqual(accounted, energy["incident_energy_MeV"], places=9)


if __name__ == "__main__":
    unittest.main()
