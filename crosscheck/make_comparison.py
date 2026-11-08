"""Assemble AHRT-OpenMC-Crosscheck-01 deliverables from Stage A + OpenMC outputs.

Outputs (results/openmc_crosscheck/):
  interaction_channel_comparison.csv
  conditional_unconditional_probabilities.json
  energy_balance_report.json

Reads:
  Stage A: results/cs137_reference/run_summary.json (+ photon_interactions.csv)
  OpenMC:  results/openmc_crosscheck/tallies_{full,aimed}.json
"""
import csv, json, math, os

BASE = "/Users/aurascoper/Developer/AHRT"
A = os.path.join(BASE, "results", "cs137_reference", "run_summary.json")
OC = os.path.join(BASE, "results", "openmc_crosscheck")

stage_a = json.load(open(A))
om_full = json.load(open(os.path.join(OC, "tallies_full.json")))
om_aimed = json.load(open(os.path.join(OC, "tallies_aimed.json")))

# ---------- Stage A channel fractions ----------
ta = stage_a["tallies"]
a_int = ta["photon_interactions"]
a_comp = ta["compton_interactions"] / a_int
a_pe = ta["photoelectric_interactions"] / a_int

# ---------- OpenMC macro-channel fractions (full run) ----------
mac = om_full["tallies"]["macro"]
# OpenMC scores: 'scatter' = coherent+incoherent; separate available
incoh = mac["incoherent-scatter"]["mean"]; coh = mac["coherent-scatter"]["mean"]
pe = mac["photoelectric"]["mean"]; pair = mac["pair-production"]["mean"]
tot_int = incoh + coh + pe + pair
om_comp = incoh / tot_int
om_pe = pe / tot_int
channels = [
    ("compton_incoherent_frac", a_comp, om_comp),
    ("photoelectric_frac", a_pe, om_pe),
    ("coherent_scatter_frac", 0.0, coh / tot_int),
    ("pair_production_frac", 0.0, pair / tot_int if tot_int else 0.0),
]
with open(os.path.join(OC, "interaction_channel_comparison.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["channel", "stage_A_fraction", "openmc_fraction", "note"])
    for name, a, o in channels:
        note = "Rayleigh(coherent) REPORTED SEPARATELY; Stage A excludes it by design" if name == "coherent_scatter_frac" else ""
        w.writerow([name, f"{a:.6f}", f"{o:.6f}", note])

# ---------- Probability estimators ----------
# Stage A unconditional measured: 0 direct cell interactions / 100000 histories
a_u_phat, a_u_n = 0.0, ta["histories"]
# 95% Poisson upper bound for zero counts
a_u_hi = -math.log(0.05) / a_u_n
# OpenMC full (unconditional): cell tally total interactions per source photon
cellf = om_full["tallies"]["cell"]
om_u_mean = sum(cellf[s]["mean"] for s in ["scatter", "photoelectric", "pair-production"])
om_u_std = math.sqrt(sum(cellf[s]["std"] ** 2 for s in ["scatter", "photoelectric", "pair-production"]))
# OpenMC aimed (conditional): cell tally per aimed photon
cella = om_aimed["tallies"]["cell"]
om_c_mean = sum(cella[s]["mean"] for s in ["scatter", "photoelectric", "pair-production"])
om_c_std = math.sqrt(sum(cella[s]["std"] ** 2 for s in ["scatter", "photoelectric", "pair-production"]))
# Stage A analytic conditional for aimed-beam geometry: E[1-exp(-mu*chord)] with chord through cell
mu = 0.08563156600105291
Rc, halfL = 5e-4, 1e-3
# y,z uniform over cell footprint: chord = 2*sqrt(Rc^2-y^2) (z always within +-1e-3 = CELL_L/2)
import statistics
rngn = 2_000_000
mean_factor = 0.0
import random
random.seed(1)
for _ in range(rngn):
    y = random.uniform(-Rc, Rc)
    mean_factor += 1 - math.exp(-2 * mu * math.sqrt(Rc * Rc - y * y))
a_cond_analytic = mean_factor / rngn

probs = {
    "P_unconditional_StageA": {"phat": a_u_phat, "n": a_u_n, "poisson95_hi": a_u_hi,
        "note": "0 observed direct cell interactions in 100k full-beam histories"},
    "P_unconditional_OpenMC": {"mean_per_source": om_u_mean, "std": om_u_std,
        "source": "full-beam run n=%d" % om_full["histories"]},
    "P_conditional_StageA_analytic": {"value": a_cond_analytic,
        "note": "E[1-exp(-mu*chord)] over aimed beam footprint; mu from imported XCOM row",
        "frozen_baseline_number": 8.6e-5},
    "P_conditional_OpenMC": {"mean_per_source": om_c_mean, "std": om_c_std,
        "source": "aimed-beam run n=%d" % om_aimed["histories"]},
    "definitions": {
        "P_conditional": "P(primary photon interacts in ROI | photon crosses ROI footprint)",
        "P_unconditional": "P(primary photon interacts in ROI | photon launched from full source)"},
}
json.dump(probs, open(os.path.join(OC, "conditional_unconditional_probabilities.json"), "w"), indent=2)

# ---------- Energy accounting ----------
acct_a = stage_a["energy_accounting"]
mac_heat = mac["heating"]["mean"]  # eV per source photon
cell_heat = cellf["heating"]["mean"]
eb = {
  "stage_A": acct_a,
  "openmc_full": {
     "macro_heating_eV_per_source": mac_heat,
     "cell_heating_eV_per_source": cell_heat,
     "note": "OpenMC led deposits recoil-electron energy locally; NOT a transported dose"},
  "closure": {
     "stage_A_incident_MeV": acct_a["incident_energy_MeV"],
     "stage_A_residuum_MeV": acct_a["incident_energy_MeV"] - acct_a["escaped_photon_energy_MeV"]
                               - acct_a["compton_electron_source_energy_MeV"] - acct_a["photoelectric_local_energy_MeV"],
     "target": "<=0.5% discrepancies subject to scoring-definition differences"},
}
json.dump(eb, open(os.path.join(OC, "energy_balance_report.json"), "w"), indent=2)
print("deliverables written:", os.listdir(OC))
print(json.dumps({k: v for k, v in probs.items() if "conditional" in k}, indent=1)[:800])
