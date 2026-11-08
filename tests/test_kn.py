import random
import unittest

from ahrt_stage_a.klein_nishina import (
    expected_kn_moments,
    sample_compton,
    scattered_energy_MeV,
)


class KleinNishinaTests(unittest.TestCase):
    def test_energy_angle_identity(self):
        rng = random.Random(1234)
        energy = 0.661657
        for _ in range(20_000):
            sample = sample_compton(energy, rng)
            expected = scattered_energy_MeV(energy, sample.cos_theta)
            self.assertAlmostEqual(sample.photon_energy_out_MeV, expected, places=14)
            self.assertAlmostEqual(
                sample.electron_kinetic_energy_MeV + sample.photon_energy_out_MeV,
                energy,
                places=14,
            )

    def test_cs137_compton_edge_and_backscatter_energy(self):
        energy = 0.661657
        backscatter = scattered_energy_MeV(energy, -1.0)
        edge = energy - backscatter
        self.assertAlmostEqual(backscatter * 1000.0, 184.3, delta=0.2)
        self.assertAlmostEqual(edge * 1000.0, 477.3, delta=0.2)

    def test_sampler_moments_against_quadrature(self):
        energy = 0.661657
        expected = expected_kn_moments(energy, n_intervals=20_000)
        rng = random.Random(81723)
        n = 80_000
        mus = []
        energies = []
        for _ in range(n):
            sample = sample_compton(energy, rng)
            mus.append(sample.cos_theta)
            energies.append(sample.photon_energy_out_MeV)
        mean_mu = sum(mus) / n
        mean_e = sum(energies) / n
        var_mu = sum((x - mean_mu) ** 2 for x in mus) / (n - 1)
        var_e = sum((x - mean_e) ** 2 for x in energies) / (n - 1)
        self.assertLess(abs(mean_mu - expected["mean_cos_theta"]), 5.0 * (var_mu / n) ** 0.5)
        self.assertLess(abs(mean_e - expected["mean_photon_energy_out_MeV"]), 5.0 * (var_e / n) ** 0.5)


if __name__ == "__main__":
    unittest.main()
