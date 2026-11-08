import random
import unittest

from ahrt_stage_a.geometry import Cylinder


class GeometryTests(unittest.TestCase):
    def test_axis_aligned_chord(self):
        c = Cylinder(radius_cm=1.0, length_cm=2.0)
        distance = c.distance_to_exit((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0))
        self.assertAlmostEqual(distance, 2.0, places=12)

    def test_parallel_beam_entry_is_on_surface(self):
        c = Cylinder(radius_cm=2.0, length_cm=4.0)
        p = c.sample_parallel_beam_entry(random.Random(3))
        self.assertAlmostEqual(p[0] ** 2 + p[1] ** 2, 4.0, places=12)
        self.assertLessEqual(abs(p[2]), 2.0)


if __name__ == "__main__":
    unittest.main()
