"""G0/G1 validation helpers for Stage A."""

from __future__ import annotations

import math
import random
from pathlib import Path

from .klein_nishina import expected_kn_moments, sample_compton
from .xcom import XCOMTable, coefficient_closure_error, load_and_verify_manifest


def validate_xcom(
    table_path: str | Path,
    manifest_path: str | Path,
    energy_MeV: float,
    density_g_cm3: float = 1.0,
    slab_thickness_cm: float = 1.0,
    histories: int = 3_000_000,
    seed: int = 445566,
) -> dict:
    """G0 provenance + G1 XCOM coefficient closure and attenuation sampling.

    The attenuation sampling gate compares a sampled exponential free-path
    probability with Beer-Lambert transmission calculated from the imported
    XCOM coefficient. It tests the transport implementation, not XCOM itself
    against a second independent database.
    """
    if density_g_cm3 <= 0 or slab_thickness_cm <= 0 or histories <= 0:
        raise ValueError("density, slab thickness and histories must be positive.")
    manifest = load_and_verify_manifest(table_path, manifest_path)
    table = XCOMTable.from_csv(table_path)
    coeff = table.at(energy_MeV)
    no_coh, with_coh = coefficient_closure_error(coeff)
    raw_row_closure = table.max_row_closure_errors()

    mu_linear_cm_inv = coeff.total_without_coherent * density_g_cm3
    optical_depth = mu_linear_cm_inv * slab_thickness_cm
    expected_interaction_probability = -math.expm1(-optical_depth)
    rng = random.Random(seed)
    interactions = 0
    for _ in range(histories):
        free_path = -math.log(max(rng.random(), 1e-300)) / mu_linear_cm_inv
        if free_path < slab_thickness_cm:
            interactions += 1
    observed_probability = interactions / histories
    attenuation_relative_error = abs(observed_probability - expected_interaction_probability) / max(
        expected_interaction_probability, 1e-15
    )

    closure_pass = max(
        raw_row_closure["max_relative_error_without_coherent"],
        raw_row_closure["max_relative_error_with_coherent"],
        no_coh, with_coh,
    ) <= 0.01
    attenuation_pass = attenuation_relative_error <= 0.01
    result = {
        "gate": "G0/G1-XCOM",
        "status": "PASS" if closure_pass and attenuation_pass else "FAIL",
        "material": manifest["material"],
        "source": manifest["source"],
        "data_sha256": manifest["sha256"],
        "energy_MeV": energy_MeV,
        "interpolation": "exact row" if any(math.isclose(r.energy_MeV, energy_MeV, rel_tol=1e-12) for r in table.rows) else "log-log linear between adjacent XCOM rows",
        "mu_total_without_coherent_cm2_g": coeff.total_without_coherent,
        "mu_total_with_coherent_cm2_g": coeff.total_with_coherent,
        "mu_partial_sum_without_coherent_cm2_g": coeff.partial_sum_without_coherent,
        "relative_closure_error_without_coherent": no_coh,
        "relative_closure_error_with_coherent_at_requested_energy": with_coh,
        "raw_imported_row_closure": raw_row_closure,
        "closure_acceptance_relative_error": 0.01,
        "mean_free_path_cm": 1.0 / mu_linear_cm_inv,
        "attenuation_test": {
            "density_g_cm3": density_g_cm3,
            "slab_thickness_cm": slab_thickness_cm,
            "linear_mu_cm_inv": mu_linear_cm_inv,
            "histories": histories,
            "seed": seed,
            "expected_interaction_probability_beer_lambert": expected_interaction_probability,
            "observed_interaction_probability": observed_probability,
            "relative_error": attenuation_relative_error,
            "acceptance_relative_error": 0.01,
        },
        "scope_note": "Internal closure plus sampled attenuation against Beer-Lambert using the imported XCOM coefficient. This is not an independent re-query or experimental validation of NIST XCOM.",
    }
    return result


def validate_kn(energy_MeV: float = 0.661657, samples: int = 150_000, seed: int = 81723) -> dict:
    if samples < 10_000:
        raise ValueError("Use at least 10,000 samples for the KN validation gate.")
    rng = random.Random(seed)
    expected = expected_kn_moments(energy_MeV)
    sum_mu = sum_mu2 = sum_e = sum_e2 = 0.0
    max_energy_identity_error = 0.0
    for _ in range(samples):
        s = sample_compton(energy_MeV, rng)
        sum_mu += s.cos_theta
        sum_mu2 += s.cos_theta * s.cos_theta
        sum_e += s.photon_energy_out_MeV
        sum_e2 += s.photon_energy_out_MeV ** 2
        identity = abs(
            s.photon_energy_out_MeV
            - energy_MeV / (1.0 + energy_MeV / 0.51099895000 * (1.0 - s.cos_theta))
        )
        max_energy_identity_error = max(max_energy_identity_error, identity)
    mean_mu = sum_mu / samples
    mean_e = sum_e / samples
    var_mu = max(0.0, sum_mu2 / samples - mean_mu * mean_mu)
    var_e = max(0.0, sum_e2 / samples - mean_e * mean_e)
    se_mu = math.sqrt(var_mu / samples)
    se_e = math.sqrt(var_e / samples)
    z_mu = abs(mean_mu - expected["mean_cos_theta"]) / max(se_mu, 1e-15)
    z_e = abs(mean_e - expected["mean_photon_energy_out_MeV"]) / max(se_e, 1e-15)
    passed = max(z_mu, z_e) <= 5.0 and max_energy_identity_error <= 1e-12
    return {
        "gate": "G1-Klein-Nishina",
        "status": "PASS" if passed else "FAIL",
        "energy_in_MeV": energy_MeV,
        "samples": samples,
        "seed": seed,
        "sample_mean_cos_theta": mean_mu,
        "quadrature_mean_cos_theta": expected["mean_cos_theta"],
        "cos_theta_difference_in_standard_errors": z_mu,
        "sample_mean_photon_energy_out_MeV": mean_e,
        "quadrature_mean_photon_energy_out_MeV": expected["mean_photon_energy_out_MeV"],
        "energy_difference_in_standard_errors": z_e,
        "max_energy_angle_identity_error_MeV": max_energy_identity_error,
        "acceptance": "Both sampled moments within 5 Monte Carlo standard errors of numerical KN quadrature; event energy-angle identity <= 1e-12 MeV.",
        "compton_backscatter_photon_energy_MeV": energy_MeV / (1.0 + 2.0 * energy_MeV / 0.51099895000),
        "compton_edge_electron_kinetic_energy_MeV": energy_MeV - energy_MeV / (1.0 + 2.0 * energy_MeV / 0.51099895000),
        "model_scope": "Unpolarized free-electron Klein-Nishina angular law; no bound-electron scattering function.",
    }
