#!/usr/bin/env python3
"""Render the 6TAP interface figure from the Level 3 receipt and the ChimeraX POV export.

Every plotted number is read from one of the two input files; none is typed in here.
POV coordinates are taken to be in Angstrom (ChimeraX export) and shown in nm.

Usage (needs matplotlib):
  python scripts/render_6tap_interface_figure.py --pov /path/to/ARC-AHRT.pov
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECEIPT = ROOT / "artifacts" / "m7_arc_level3_interface_comparison_receipt.json"
DEFAULT_OUT = ROOT / "docs" / "6tap_interface_atlas_figure.png"

SLAB_HALF_NM = 2.0  # half-thickness of the central slab shown in panel D

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
CONTACT = [  # (receipt key, legend label, colour)
    ("NTD_NTD_pairs", "NTD-NTD (same domain)", "#2a78d6"),
    ("CTD_CTD_pairs", "CTD-CTD (same domain)", "#eb6834"),
    ("NTD_CTD_pairs", "NTD-CTD (cross domain)", "#1baf7a"),
]
NC_COLOR, ZN_COLOR = "#6250d6", "#eda100"

SPHERE = re.compile(r"sphere\{<([-\d.eE]+),([-\d.eE]+),([-\d.eE]+)>,\s*([\d.]+)")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interface_rows(receipt: dict) -> list[dict]:
    profiles = receipt["normalized_interface_chemistry_and_area_profiles"]["subinterface_profiles"]
    rows = []
    for p in profiles:
        dap = p["domain_architecture_profile"]
        descriptor = " ".join(p["interface_class"].split("_")[2:-1])
        rows.append({
            "name": p["interface_class"],
            "family": p["interface_fingerprint_family"].split("_")[1],
            "label": f"{p['chain_pair']} {p['symmetry_symbol'].split('_')[0]} {descriptor}",
            **{key: dap[key] for key, _, _ in CONTACT},
        })
    return rows


def receipt_gates(rows: list[dict], receipt: dict) -> None:
    """Abort if independent sections of the receipt disagree about the numbers plotted.

    Row gate: each row's contact total (NN + CC + CN, from the chemistry profiles)
    equals the residue-pair count recorded separately under
    symmetry_normalized_residue_fingerprints.
    Sum gate: the rows summed over all sub-interfaces equal the E_NN, E_CC and E_CN
    recorded in the stub-identity section, and that section's derived values
    (E_same, conserved E_CC - E_NN) follow from them.
    """
    fingerprints = receipt["symmetry_normalized_residue_fingerprints"]["subinterface_fingerprints"]
    pairs = {f["interface_class"]: f["unique_residue_pair_count_at_5A"] for f in fingerprints}
    for r in rows:
        total = sum(r[key] for key, _, _ in CONTACT)
        if total != pairs[r["name"]]:
            sys.exit(f"receipt disagrees with itself: {r['name']} profile total {total} "
                     f"!= fingerprint pair count {pairs[r['name']]}")
    stub = receipt["constrained_graph_state_space_connectivity_and_mcmc_audit"][
        "stub_difference_conservation_identity"]["empirical_6tap_5A_cutoff_values"]
    nn, cc, cn = (sum(r[key] for r in rows) for key, _, _ in CONTACT)
    derived = {"observed_E_NN": nn, "observed_E_CC": cc, "observed_E_CN": cn,
               "observed_E_same": nn + cc, "conserved_E_CC_minus_E_NN": cc - nn}
    recorded = {key: stub[key] for key in derived}
    if derived != recorded:
        sys.exit(f"receipt disagrees with itself: summed rows give {derived}, "
                 f"stub-identity section records {recorded}")


def expected_sphere_counts(receipt: dict) -> dict[str, int]:
    """Residue spheres the full assembly should contain, per class, from the mmCIF inventory."""
    prov = receipt["mmcif_provenance"]
    operators = prov["inventory_reconciliation"]["operators_count"]
    role = {"CA_CAPSID_SHELL": "capsid", "NC_ZINC_KNUCKLE_FRAGMENT": "nc", "ZN_COORDINATION_LIGAND": "zn"}
    counts: dict[str, int] = defaultdict(int)
    for asym in prov["struct_asym_inventory"]:
        counts[role[asym["functional_component_id"]]] += asym["residue_count_per_asu"] * operators
    return dict(counts)


def read_pov(path: Path, expected: dict[str, int]) -> dict[str, list[tuple]]:
    """Spheres per class as (x, y, z, radius, line) in nm, centred on the assembly centroid.

    The export carries no labels, so classes are assigned by matching each
    sphere-radius group's count to the receipt's inventory; a mismatch aborts.
    """
    by_radius: dict[float, list[tuple]] = defaultdict(list)
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        m = SPHERE.search(line)
        if m:
            x, y, z, r = map(float, m.groups())
            by_radius[round(r, 2)].append((x, y, z, lineno))
    everything = [s for group in by_radius.values() for s in group]
    centre = [sum(s[i] for s in everything) / len(everything) for i in range(3)]
    classes: dict[str, list[tuple]] = defaultdict(list)
    for group in by_radius.values():
        name = next((k for k in ("zn", "nc") if len(group) == expected[k]), "capsid")
        for x, y, z, lineno in group:
            dx, dy, dz = ((v - c) * 0.1 for v, c in zip((x, y, z), centre))
            classes[name].append((dx, dy, dz, math.sqrt(dx * dx + dy * dy + dz * dz), lineno))
    found = {k: len(v) for k, v in classes.items()}
    if found != expected:
        sys.exit(f"POV sphere classes {found} do not match receipt inventory {expected}")
    return dict(classes)


def style(ax, title: str, subtitle: str = "") -> None:
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold", color=INK, pad=22 if subtitle else 8)
    if subtitle:
        ax.text(0, 1.02, subtitle, transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c3c2b7")
    ax.tick_params(colors=MUTED, labelsize=8.5, length=0)


def panel_interfaces(ax, rows: list[dict]) -> None:
    rows = sorted(rows, key=lambda r: (r["family"], -sum(r[k] for k, _, _ in CONTACT)))[::-1]
    left = [0] * len(rows)
    for key, label, colour in CONTACT:
        vals = [r[key] for r in rows]
        ax.barh(range(len(rows)), vals, left=left, color=colour, edgecolor=SURFACE, linewidth=1.5,
                height=0.62, label=label)
        left = [a + b for a, b in zip(left, vals)]
    for i, r in enumerate(rows):
        same = r["NTD_NTD_pairs"] + r["CTD_CTD_pairs"]
        ax.text(left[i] + 0.8, i, f"{same}/{left[i]}", va="center", fontsize=8.5, color=INK)
    ax.set_yticks(range(len(rows)), [f"{r['label']}  ({r['family']})" for r in rows])
    ax.tick_params(axis="y", labelcolor=INK2)
    ax.set_xlim(0, max(left) * 1.14)
    ax.set_xlabel("Residue-pair contacts at d ≤ 5.0 Å", fontsize=9, color=INK2)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    style(ax, "A. Contacts per sub-interface", "End label: same-domain / total. Family in parentheses.")


def panel_families(ax, rows: list[dict], stub: dict) -> None:
    groups = sorted({r["family"] for r in rows}) + ["All"]
    bottom = [0] * len(groups)
    for key, _, colour in CONTACT:
        vals = [sum(r[key] for r in rows if g in (r["family"], "All")) for g in groups]
        ax.bar(range(len(groups)), vals, bottom=bottom, color=colour, edgecolor=SURFACE,
               linewidth=1.5, width=0.5)
        for i, v in enumerate(vals):
            if v >= 12:
                ax.text(i, bottom[i] + v / 2, str(v), ha="center", va="center", fontsize=8.5,
                        color="#ffffff")
        bottom = [a + b for a, b in zip(bottom, vals)]
    for i, g in enumerate(groups):
        same = sum(r["NTD_NTD_pairs"] + r["CTD_CTD_pairs"] for r in rows if g in (r["family"], "All"))
        ax.text(i, bottom[i] + 3, f"{same}/{bottom[i]}", ha="center", fontsize=8.5, color=INK)
    ax.set_xticks(range(len(groups)), [g if g == "All" else f"Family {g}" for g in groups])
    ax.tick_params(axis="x", labelcolor=INK2)
    ax.set_ylim(0, max(bottom) * 1.1)
    ax.set_ylabel("Residue-pair contacts", fontsize=9, color=INK2)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    style(ax, "B. Domain split by interface family",
          f"Conserved under both nulls: E_CC − E_NN = {stub['observed_E_CC']} − {stub['observed_E_NN']}"
          f" = {stub['conserved_E_CC_minus_E_NN']}. Colours as in A.")


def panel_nulls(ax, audit: dict) -> None:
    schemes = [
        ("Scheme I\nglobal degree-preserving", audit["global_degree_preserving_simple_graph_double_edge_swap_mcmc"]),
        ("Scheme II\nwithin sub-interface", audit["within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc"]),
    ]
    observed = schemes[0][1]["observed_same_domain_contacts"]
    for i, (_, s) in enumerate(schemes):
        mean, sd = s["mcmc_mean_same_domain_contacts"], s["mcmc_std_same_domain_contacts"]
        ax.plot([s["mcmc_min_same_domain_contacts"], s["mcmc_max_same_domain_contacts"]], [i, i],
                color=MUTED, linewidth=1.2, solid_capstyle="butt")
        ax.plot([mean - sd, mean + sd], [i, i], color=INK2, linewidth=5, solid_capstyle="butt")
        ax.plot(mean, i, "o", color=INK, markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.5)
        note = (f"mean {mean:.2f} ± {sd:.2f},  z = {s['same_domain_z_score']:.2f},  "
                f"p = {s['empirical_p_value_phipson_smyth']:.3g}")
        if "inferential_role" in s:
            note += "\nreceipt role: " + s["inferential_role"].replace("_", " ").lower()
        ax.text(0.02, i + 0.12, note, transform=ax.get_yaxis_transform(), fontsize=8.5, color=INK,
                va="bottom")
    ax.axvline(observed, color="#d03b3b", linewidth=1.5)
    ax.text(observed - 0.6, -0.42, f"observed {observed}", color=INK, fontsize=8.5, va="top", ha="right")
    ax.set_yticks(range(len(schemes)), [name for name, _ in schemes])
    ax.tick_params(axis="y", labelcolor=INK2)
    ax.set_ylim(-0.6, len(schemes) - 0.25)
    ax.set_xlabel("Same-domain contacts E_same", fontsize=9, color=INK2)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    style(ax, "C. Observed same-domain count against two nulls",
          "Dot: MCMC mean. Thick bar: ± 1 s.d. Thin line: sampled min to max.")


def panel_pov(ax, classes: dict[str, list[tuple]], pov_name: str) -> dict:
    span = {k: (min(s[3] for s in v), max(s[3] for s in v)) for k, v in classes.items()}
    slab = {k: [s for s in v if abs(s[2]) < SLAB_HALF_NM] for k, v in classes.items()}
    ax.scatter([s[0] for s in slab["capsid"]], [s[1] for s in slab["capsid"]], s=5, color=MUTED,
               linewidths=0, label=f"Capsid residues ({len(classes['capsid']):,})")
    ax.scatter([s[0] for s in slab["nc"]], [s[1] for s in slab["nc"]], s=16, color=NC_COLOR,
               linewidths=0, label=f"NC residues ({len(classes['nc']):,}), centres "
                                    f"{span['nc'][0]:.1f} to {span['nc'][1]:.1f} nm")
    ax.scatter([s[0] for s in slab["zn"]], [s[1] for s in slab["zn"]], s=60, color=ZN_COLOR,
               edgecolors=SURFACE, linewidths=1.5, label=f"Zn ({len(classes['zn'])}), {span['zn'][0]:.1f} nm")
    for radius in span["capsid"]:
        ax.add_patch(plt.Circle((0, 0), radius, fill=False, color=INK2, linewidth=0.9, linestyle=(0, (4, 3))))
    inner, outer = span["capsid"]
    ax.text(0, inner - 1.0, f"{inner:.1f} nm", ha="center", va="top", fontsize=8.5, color=INK)
    ax.text(0, outer + 0.6, f"{outer:.1f} nm", ha="center", va="bottom", fontsize=8.5, color=INK)
    ax.margins(0.1)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x (nm)", fontsize=9, color=INK2)
    ax.set_ylabel("y (nm)", fontsize=9, color=INK2)
    legend = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=8, frameon=False)
    for handle in legend.legend_handles:
        handle.set_sizes([36])
    style(ax, f"D. Central slab of {pov_name}",
          f"|z| < {SLAB_HALF_NM:g} nm, one sphere per residue. Dashed: innermost and outermost capsid centres.")
    return {k: {"count": len(v), "centre_radius_nm": [round(x, 3) for x in span[k]],
                "lines": [min(s[4] for s in v), max(s[4] for s in v)], "in_slab": len(slab[k])}
            for k, v in classes.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pov", type=Path, required=True, help="ChimeraX POV-Ray export of the 6TAP assembly")
    ap.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    receipt = json.loads(args.receipt.read_text())
    rows = interface_rows(receipt)
    receipt_gates(rows, receipt)
    audit = receipt["constrained_graph_state_space_connectivity_and_mcmc_audit"]
    stub = audit["stub_difference_conservation_identity"]["empirical_6tap_5A_cutoff_values"]
    classes = read_pov(args.pov, expected_sphere_counts(receipt))
    hashes = {args.receipt.name: sha256(args.receipt), args.pov.name: sha256(args.pov)}

    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5), dpi=160, facecolor=SURFACE,
                             gridspec_kw={"width_ratios": [1.25, 1]}, layout="constrained")
    panel_interfaces(axes[0][0], rows)
    panel_families(axes[0][1], rows, stub)
    panel_nulls(axes[1][0], audit)
    pov_summary = panel_pov(axes[1][1], classes, args.pov.name)
    fig.suptitle(f"PDB {receipt['mmcif_provenance']['structure_id']} capsid: interface contacts, "
                 "null-model comparison and assembly cross-section",
                 x=0.01, ha="left", fontsize=13, fontweight="bold", color=INK)
    fig.supxlabel("Sources: " + ";  ".join(f"{n} sha256 {h[:12]}" for n, h in hashes.items()),
                  x=0.01, ha="left", fontsize=8, color=MUTED)
    fig.savefig(args.out)
    args.out.with_suffix(".receipt.json").write_text(json.dumps({
        "schema": "arc-ahrt-6tap-interface-figure-receipt-v1",
        "figure": {"file": args.out.name, "sha256": sha256(args.out)},
        "inputs_sha256": hashes,
        "reproducibility_note": f"{args.pov.name} is a ChimeraX export that is not tracked in this repository; "
                                "the figure can be regenerated only with the file matching its sha256 above.",
        "slab_half_thickness_nm": SLAB_HALF_NM,
        "pov_sphere_classes": pov_summary,
    }, indent=2) + "\n")
    print(f"wrote {args.out} and {args.out.with_suffix('.receipt.json').name}")


if __name__ == "__main__":
    main()
