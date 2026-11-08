"""Nested finite-cylinder geometry for the Stage A photon domain."""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Cylinder:
    radius_cm: float
    length_cm: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.radius_cm) or self.radius_cm <= 0:
            raise ValueError("Cylinder radius must be finite and positive.")
        if not math.isfinite(self.length_cm) or self.length_cm <= 0:
            raise ValueError("Cylinder length must be finite and positive.")

    def contains(self, p: tuple[float, float, float], atol: float = 1e-12) -> bool:
        x, y, z = p
        return x * x + y * y <= (self.radius_cm + atol) ** 2 and abs(z) <= self.length_cm / 2 + atol

    def distance_to_exit(
        self,
        p: tuple[float, float, float],
        d: tuple[float, float, float],
    ) -> float:
        """Distance to first positive intersection with radial wall or endcap."""
        x, y, z = p
        dx, dy, dz = d
        candidates: list[float] = []
        a = dx * dx + dy * dy
        b = 2.0 * (x * dx + y * dy)
        c = x * x + y * y - self.radius_cm * self.radius_cm
        if a > 1e-30:
            disc = b * b - 4.0 * a * c
            if disc >= 0.0:
                root = math.sqrt(max(0.0, disc))
                for t in ((-b - root) / (2.0 * a), (-b + root) / (2.0 * a)):
                    if t > 1e-10:
                        zt = z + t * dz
                        if abs(zt) <= self.length_cm / 2 + 1e-10:
                            candidates.append(t)
        if abs(dz) > 1e-30:
            for zcap in (-self.length_cm / 2, self.length_cm / 2):
                t = (zcap - z) / dz
                if t > 1e-10:
                    xt, yt = x + t * dx, y + t * dy
                    if xt * xt + yt * yt <= self.radius_cm * self.radius_cm + 1e-10:
                        candidates.append(t)
        if not candidates:
            raise RuntimeError("No positive exit intersection found; point may be outside domain.")
        return min(candidates)

    def sample_parallel_beam_entry(self, rng) -> tuple[float, float, float]:
        """Uniform fluence over the cylinder's projected y-z rectangle."""
        y = rng.uniform(-self.radius_cm, self.radius_cm)
        z = rng.uniform(-self.length_cm / 2, self.length_cm / 2)
        x = -math.sqrt(max(0.0, self.radius_cm**2 - y**2))
        return (x, y, z)
