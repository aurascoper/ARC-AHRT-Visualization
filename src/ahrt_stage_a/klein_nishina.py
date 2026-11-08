"""Unpolarized free-electron Klein-Nishina angular sampling."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

ELECTRON_REST_ENERGY_MEV = 0.51099895000


@dataclass(frozen=True)
class ComptonSample:
    cos_theta: float
    phi_rad: float
    photon_energy_out_MeV: float
    electron_kinetic_energy_MeV: float


def scattered_energy_MeV(energy_in_MeV: float, cos_theta: float) -> float:
    if energy_in_MeV <= 0 or not math.isfinite(energy_in_MeV):
        raise ValueError("Incident photon energy must be finite and positive.")
    if not -1.0 <= cos_theta <= 1.0:
        raise ValueError("cos(theta) must lie in [-1, 1].")
    alpha = energy_in_MeV / ELECTRON_REST_ENERGY_MEV
    return energy_in_MeV / (1.0 + alpha * (1.0 - cos_theta))


def kn_pdf_unnormalized_mu(energy_in_MeV: float, cos_theta: float) -> float:
    """Differential KN kernel dσ/dΩ divided by r_e²; azimuth-independent."""
    e_ratio = scattered_energy_MeV(energy_in_MeV, cos_theta) / energy_in_MeV
    sin2 = max(0.0, 1.0 - cos_theta * cos_theta)
    return 0.5 * e_ratio * e_ratio * (e_ratio + 1.0 / e_ratio - sin2)


def sample_compton(energy_in_MeV: float, rng: random.Random) -> ComptonSample:
    """Rejection-sample cos(theta) from the free-electron KN law.

    Uniform proposals in mu=cos(theta) and uniform azimuth sample solid angle.
    The KN kernel divided by r_e^2 is bounded by one, so it is a valid
    rejection probability for this proposal distribution.
    """
    if energy_in_MeV <= 0 or not math.isfinite(energy_in_MeV):
        raise ValueError("Incident photon energy must be finite and positive.")
    while True:
        mu = 2.0 * rng.random() - 1.0
        if rng.random() <= kn_pdf_unnormalized_mu(energy_in_MeV, mu):
            break
    phi = 2.0 * math.pi * rng.random()
    e_out = scattered_energy_MeV(energy_in_MeV, mu)
    return ComptonSample(
        cos_theta=mu,
        phi_rad=phi,
        photon_energy_out_MeV=e_out,
        electron_kinetic_energy_MeV=energy_in_MeV - e_out,
    )


def expected_kn_moments(energy_in_MeV: float, n_intervals: int = 100_000) -> dict[str, float]:
    """Numerically integrate exact KN angular law over mu for test oracles."""
    if n_intervals < 1000:
        raise ValueError("Use at least 1000 intervals for quadrature.")
    if n_intervals % 2:
        n_intervals += 1
    h = 2.0 / n_intervals
    sums = {"norm": 0.0, "mu": 0.0, "eout": 0.0, "mu2": 0.0, "eout2": 0.0}
    for i in range(n_intervals + 1):
        mu = -1.0 + i * h
        w = 1.0 if i in (0, n_intervals) else (4.0 if i % 2 else 2.0)
        kernel = kn_pdf_unnormalized_mu(energy_in_MeV, mu)
        eout = scattered_energy_MeV(energy_in_MeV, mu)
        sums["norm"] += w * kernel
        sums["mu"] += w * kernel * mu
        sums["eout"] += w * kernel * eout
        sums["mu2"] += w * kernel * mu * mu
        sums["eout2"] += w * kernel * eout * eout
    norm = sums["norm"] * h / 3.0
    mean_mu = sums["mu"] * h / 3.0 / norm
    mean_eout = sums["eout"] * h / 3.0 / norm
    mean_mu2 = sums["mu2"] * h / 3.0 / norm
    mean_eout2 = sums["eout2"] * h / 3.0 / norm
    return {
        "mean_cos_theta": mean_mu,
        "mean_photon_energy_out_MeV": mean_eout,
        "variance_cos_theta": max(0.0, mean_mu2 - mean_mu * mean_mu),
        "variance_photon_energy_out_MeV2": max(0.0, mean_eout2 - mean_eout * mean_eout),
    }


def rotate_direction(direction: tuple[float, float, float], cos_theta: float, phi: float) -> tuple[float, float, float]:
    """Rotate a direction by polar angle arccos(cos_theta) and azimuth phi."""
    dx, dy, dz = direction
    norm = math.sqrt(dx * dx + dy * dy + dz * dz)
    if norm == 0.0:
        raise ValueError("Direction vector cannot be zero.")
    dx, dy, dz = dx / norm, dy / norm, dz / norm

    # Select a stable reference not parallel to the incoming direction.
    if abs(dz) < 0.9:
        ref = (0.0, 0.0, 1.0)
    else:
        ref = (0.0, 1.0, 0.0)
    # e1 = normalize(ref x d); e2 = d x e1
    e1 = (
        ref[1] * dz - ref[2] * dy,
        ref[2] * dx - ref[0] * dz,
        ref[0] * dy - ref[1] * dx,
    )
    e1n = math.sqrt(sum(v * v for v in e1))
    e1 = tuple(v / e1n for v in e1)
    e2 = (
        dy * e1[2] - dz * e1[1],
        dz * e1[0] - dx * e1[2],
        dx * e1[1] - dy * e1[0],
    )
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    cp, sp = math.cos(phi), math.sin(phi)
    d = (dx, dy, dz)
    out = tuple(
        cos_theta * d[i] + sin_theta * (cp * e1[i] + sp * e2[i])
        for i in range(3)
    )
    out_norm = math.sqrt(sum(v * v for v in out))
    return tuple(v / out_norm for v in out)


def electron_direction_from_compton(
    photon_direction_in: tuple[float, float, float],
    photon_direction_out: tuple[float, float, float],
    energy_in_MeV: float,
    energy_out_MeV: float,
) -> tuple[float, float, float]:
    """Electron recoil direction from momentum conservation (p*c units)."""
    vec = tuple(
        energy_in_MeV * photon_direction_in[i] - energy_out_MeV * photon_direction_out[i]
        for i in range(3)
    )
    norm = math.sqrt(sum(v * v for v in vec))
    if norm == 0.0:
        # Forward scatter at zero recoil is a measure-zero limiting case.
        return tuple(photon_direction_in[i] for i in range(3))
    return tuple(v / norm for v in vec)
