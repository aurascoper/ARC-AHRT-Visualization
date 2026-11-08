import math, unittest
from ahrt_stage_a.klein_nishina import rotate_direction

class RotateDirectionRegression(unittest.TestCase):
    def test_rotate_direction_preserves_angle_and_norm(self):
        d = (1.0, 0.0, 0.0)
        out = rotate_direction(d, 0.5, 1.1)
        self.assertAlmostEqual(sum(v*v for v in out), 1.0, places=12)
        dot = sum(out[i]*d[i] for i in range(3))
        self.assertAlmostEqual(dot, 0.5, places=12)
    def test_rotate_direction_off_axis(self):
        d = (0.267261241662424, 0.534522483824849, 0.801783725737273)
        out = rotate_direction(d, -0.3, 4.0)
        self.assertAlmostEqual(sum(v*v for v in out), 1.0, places=12)
        self.assertAlmostEqual(sum(out[i]*d[i] for i in range(3)), -0.3, places=9)
