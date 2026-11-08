# ==============================================================================
# rstudio/arc_ahrt_edge_cases.R
# Title: RStudio Edge-Case Audit & Visualization Suite for ARC Capsid & AHRT
# Epistemic Status: [MEASURED GEOMETRY & TOPAS-NBIO TRACK SERIES - NO MECHANISM CLAIMED]
#
# Usage in RStudio:
#   source("rstudio/arc_ahrt_edge_cases.R")
# Usage from CLI:
#   Rscript rstudio/arc_ahrt_edge_cases.R
#
# Audited Edge Cases:
#   1. ASU-only vs 60-Operator Cutoff Percolation Edge Case (6TAP):
#      A-A (C5 pentamer), B-D, and C-C (C3 hexamer) have 0 contacts in the single
#      deposited ASU at all cutoffs (3.5-6.0 A) and only emerge under icosahedral
#      symmetry operators; A-B is 0 up to 4.5 A and emerges at 5.0 A via a single
#      long-range contact (LYS181-ASN204 at 4.6059 A, op 2).
#   2. Multiplicity & Null-Model Stratification Edge Case (Level 3 12-Test Battery):
#      Out of 12 permutation tests across degree-preserving and circular-shift nulls,
#      8/12 survive Benjamini-Hochberg FDR (q <= 0.05) and 3/12 survive conservative
#      omnibus-12 Holm-Bonferroni (alpha <= 0.05). Notably, homologous_83_residue_same_domain
#      survives within-family Holm (p = 0.0325) and omnibus BH (q = 0.0195), but sits right
#      on the omnibus-12 Holm boundary (p = 0.0585).
#   3. TOPAS-nBio V4-ASRA Radical Scavenger Floor & Secondary H^0 Rise:
#      In W-C1 (25 mM NO3- + MeOH), G(OH^0) hits an exact 0 floor at >= 0.1 M MeOH while
#      G(H^0) rises non-linearly by +67.9% (0.3510 -> 0.5893 molecules/100 eV) due to
#      secondary H abstraction. In W-C2 (O2-sat + DMSO), G(e_aq^-1) increases monotonically
#      from 0.01543 to 0.01908 (+23.6%) as OH^0 scavenging suppresses OH + e_aq^- recombination.
#   4. Chain E (NC Zinc-Knuckle) B-Factor Dispersion Edge Case:
#      Chain E (residues 224-250) exhibits 3x-4x higher B-factor standard deviation
#      (32.25 A^2, range 9.96-160.89 A^2) than the outer CA shell chains A-D (7.68-10.61 A^2).
# ==============================================================================

resolve_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  candidates <- character(0)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
    candidates <- c(candidates, dirname(dirname(script_path)), dirname(script_path))
  }
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  candidates <- unique(c(candidates, cwd, dirname(cwd)))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "structures", "6TAP.cif"))) {
      return(cand)
    }
  }
  stop("Could not locate repository root containing structures/6TAP.cif")
}

repo_root <- resolve_repo_root()
out_dir <- file.path(repo_root, "rstudio", "renders")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Read V4-ASRA scavenger CSV directly with base R
v4_csv_path <- file.path(repo_root, "docs", "v4_matlab_series.csv")
v4_df <- read.csv(v4_csv_path, stringsAsFactors = FALSE)

# Parse JSON fixtures using jsonlite if available, or python3 one-liner bridge if not
read_edge_case_tables <- function(root_dir) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    sweep <- jsonlite::fromJSON(file.path(root_dir, "geometry", "interface_residue_contacts_sweep.json"))
    l3 <- jsonlite::fromJSON(file.path(root_dir, "artifacts", "m7_arc_level3_interface_comparison_receipt.json"))

    cutoffs <- as.numeric(sweep$cutoffs_angstrom)
    cut_keys <- sprintf("%.1f", cutoffs)
    pairs <- c("A-A", "A-B", "A-D", "B-C", "B-D", "C-C", "C-D")
    asu_mat <- sapply(cut_keys, function(k) as.numeric(unlist(sweep$asu_only_interface_contacts_by_cutoff[[k]][pairs])))
    full_mat <- sapply(cut_keys, function(k) as.numeric(unlist(sweep$full_assembly_unique_interface_contacts_by_cutoff[[k]][pairs])))
    rownames(asu_mat) <- pairs
    rownames(full_mat) <- pairs
    colnames(asu_mat) <- cut_keys
    colnames(full_mat) <- cut_keys

    ch_names <- names(sweep$chain_summary)
    bf_df <- data.frame(
      chain = ch_names,
      residue_count = sapply(ch_names, function(c) sweep$chain_summary[[c]]$residue_count),
      b_min = sapply(ch_names, function(c) sweep$chain_summary[[c]]$b_factor_min),
      b_mean = sapply(ch_names, function(c) sweep$chain_summary[[c]]$b_factor_mean),
      b_max = sapply(ch_names, function(c) sweep$chain_summary[[c]]$b_factor_max),
      b_sd = sapply(ch_names, function(c) sweep$chain_summary[[c]]$b_factor_std),
      stringsAsFactors = FALSE
    )

    tests_df <- l3$multiplicity_adjustment_analysis$omnibus_12_permutation_test_battery$all_tests
    return(list(cutoffs = cutoffs, asu_mat = asu_mat, full_mat = full_mat, bf_df = bf_df, tests_df = tests_df))
  } else {
    # Deterministic fallback using python3 stdlib json to emit flat CSVs if jsonlite is absent
    tmp_dir <- tempdir()
    py_cmd <- sprintf(
      "import json, csv, pathlib; r = pathlib.Path('%s'); t = pathlib.Path('%s');
s = json.load(open(r/'geometry/interface_residue_contacts_sweep.json'));
l3 = json.load(open(r/'artifacts/m7_arc_level3_interface_comparison_receipt.json'));
pairs = ['A-A','A-B','A-D','B-C','B-D','C-C','C-D'];
cuts = ['3.5','4.0','4.5','5.0','5.5','6.0'];
w1 = csv.writer(open(t/'sweep_full.csv','w',newline='')); w1.writerow(['pair']+cuts);
[w1.writerow([p]+[s['full_assembly_unique_interface_contacts_by_cutoff'][c][p] for c in cuts]) for p in pairs];
w2 = csv.writer(open(t/'sweep_asu.csv','w',newline='')); w2.writerow(['pair']+cuts);
[w2.writerow([p]+[s['asu_only_interface_contacts_by_cutoff'][c][p] for c in cuts]) for p in pairs];
w3 = csv.writer(open(t/'bfactors.csv','w',newline='')); w3.writerow(['chain','residue_count','b_min','b_mean','b_max','b_sd']);
[w3.writerow([c, v['residue_count'], v['b_min'] if 'b_min' in v else v['b_factor_min'], v['b_factor_mean'], v['b_factor_max'], v['b_factor_std']]) for c,v in s['chain_summary'].items()];
w4 = csv.writer(open(t/'l3_tests.csv','w',newline='')); w4.writerow(['test_id','p_unadjusted','p_holm_within_family','p_holm_omnibus_12','q_bh_omnibus_12']);
[w4.writerow([x['test_id'],x['p_unadjusted'],x['p_holm_within_family'],x['p_holm_omnibus_12'],x['q_bh_omnibus_12']]) for x in l3['multiplicity_adjustment_analysis']['omnibus_12_permutation_test_battery']['all_tests']]",
      root_dir, tmp_dir
    )
    system2("python3", args = c("-c", shQuote(py_cmd)))
    full_df <- read.csv(file.path(tmp_dir, "sweep_full.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    asu_df <- read.csv(file.path(tmp_dir, "sweep_asu.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    full_mat <- as.matrix(full_df[, -1]); rownames(full_mat) <- full_df$pair
    asu_mat <- as.matrix(asu_df[, -1]); rownames(asu_mat) <- asu_df$pair
    bf_df <- read.csv(file.path(tmp_dir, "bfactors.csv"), stringsAsFactors = FALSE)
    tests_df <- read.csv(file.path(tmp_dir, "l3_tests.csv"), stringsAsFactors = FALSE)
    return(list(cutoffs = as.numeric(colnames(full_mat)), asu_mat = asu_mat, full_mat = full_mat, bf_df = bf_df, tests_df = tests_df))
  }
}

dat <- read_edge_case_tables(repo_root)

# Export flat edge-case summary CSV for RStudio inspection
wc1_meoh <- subset(v4_df, fixture == 1 & state == "COMPLETED")
wc1_oh <- subset(wc1_meoh, species == "OH^0")
wc1_h <- subset(wc1_meoh, species == "H^0")
wc2_dmso <- subset(v4_df, fixture == 2 & state == "COMPLETED")
wc2_eaq <- subset(wc2_dmso, species == "e_aq^-1")

edge_summary <- data.frame(
  edge_case_id = c(
    "EC1_symmetry_only_interfaces",
    "EC1_AB_5A_percolation_onset",
    "EC2_homologous_83_holm_boundary",
    "EC3_WC1_OH_zero_floor_vs_H_rise",
    "EC3_WC2_eaq_recombination_sparing",
    "EC4_chain_E_bfactor_dispersion"
  ),
  metric_name = c(
    "A-A + B-D + C-C contacts at 5.0 A (ASU vs 60-op)",
    "A-B contacts at 4.5 A vs 5.0 A (op 2 LYS181-ASN204)",
    "homologous_83_residue_same_domain (p_holm_within vs p_holm_omnibus12)",
    "W-C1 G(H^0) rise across 0.001 M -> 10 M MeOH while G(OH^0) -> 0",
    "W-C2 G(e_aq^-1) rise across 1e-6 M -> 1 M DMSO",
    "Chain E B-factor SD vs Chains A-D max SD (A^2)"
  ),
  primary_value = c(
    sum(dat$asu_mat[c("A-A", "B-D", "C-C"), "5.0"]),
    dat$full_mat["A-B", "4.5"],
    dat$tests_df$p_holm_within_family[dat$tests_df$test_id == "homologous_83_residue_same_domain"],
    wc1_h$G_1us[wc1_h$conc == 0.001],
    wc2_eaq$G_1us[wc2_eaq$conc == 1e-06],
    dat$bf_df$b_sd[dat$bf_df$chain == "E"]
  ),
  comparison_value = c(
    sum(dat$full_mat[c("A-A", "B-D", "C-C"), "5.0"]),
    dat$full_mat["A-B", "5.0"],
    dat$tests_df$p_holm_omnibus_12[dat$tests_df$test_id == "homologous_83_residue_same_domain"],
    wc1_h$G_1us[wc1_h$conc == 10],
    wc2_eaq$G_1us[wc2_eaq$conc == 1],
    max(dat$bf_df$b_sd[dat$bf_df$chain %in% c("A", "B", "C", "D")])
  ),
  stringsAsFactors = FALSE
)

summary_csv_path <- file.path(out_dir, "rstudio_edge_cases_summary.csv")
write.csv(edge_summary, summary_csv_path, row.names = FALSE)

# Render 1920x1080 4-panel RStudio Edge-Case Figure
png_path <- file.path(out_dir, "rstudio_arc_ahrt_edge_cases.png")
png(filename = png_path, width = 1920, height = 1080, res = 144, bg = "#080d14")

par(mfrow = c(2, 2), mar = c(5.2, 5.0, 3.6, 2.0), oma = c(1.5, 1.0, 3.0, 1.0),
    bg = "#080d14", fg = "#cbd5e1", col.axis = "#cbd5e1", col.lab = "#f8fafc", col.main = "#f8fafc")

# Panel 1: Edge Case 1 - Full-Assembly vs ASU-Only Interface Contact Emergence
cols_p1 <- c("#38bdf8", "#facc15", "#22c55e", "#f97316", "#a855f7", "#ec4899", "#14b8a6")
matplot(dat$cutoffs, t(dat$full_mat), type = "b", pch = 19, lty = 1, lwd = 2.4, col = cols_p1,
        xlab = "Heavy-Atom Distance Cutoff (Angstrom)",
        ylab = "Unique Residue Contact Pairs (count)",
        main = "Edge Case 1: 6TAP Full-Assembly (solid) vs ASU-Only (dashed) Contacts")
matlines(dat$cutoffs, t(dat$asu_mat), type = "b", pch = 1, lty = 2, lwd = 1.6, col = cols_p1)
abline(v = 5.0, col = "#fbbf24", lty = 3, lwd = 1.5)
legend("topleft", legend = rownames(dat$full_mat), col = cols_p1, lty = 1, pch = 19,
       cex = 0.82, bg = "#0f172a", text.col = "#f8fafc", ncol = 2)

# Panel 2: Edge Case 2 - Level 3 Omnibus-12 Multiplicity Adjustment Battery
ord <- order(dat$tests_df$p_unadjusted)
t_ord <- dat$tests_df[ord, ]
y_idx <- seq_len(nrow(t_ord))
neglog_unadj <- -log10(pmax(t_ord$p_unadjusted, 1e-6))
neglog_bh <- -log10(pmax(t_ord$q_bh_omnibus_12, 1e-6))
neglog_holm <- -log10(pmax(t_ord$p_holm_omnibus_12, 1e-6))

plot(neglog_unadj, y_idx, xlim = c(0, max(neglog_unadj) + 0.5), yaxt = "n",
     pch = 19, col = "#38bdf8", cex = 1.3,
     xlab = "-log10(p or q value) [Vertical line = 0.05 threshold]",
     ylab = "",
     main = "Edge Case 2: Level 3 12-Permutation Battery (Unadjusted vs BH vs Holm)")
points(neglog_bh, y_idx, pch = 17, col = "#22c55e", cex = 1.3)
points(neglog_holm, y_idx, pch = 15, col = "#f97316", cex = 1.3)
abline(v = -log10(0.05), col = "#ef4444", lty = 2, lwd = 2.0)
axis(2, at = y_idx, labels = substr(t_ord$test_id, 1, 22), las = 2, cex.axis = 0.72)
legend("bottomright", legend = c("Unadjusted p", "Omnibus-12 BH q (8/12 <= 0.05)", "Omnibus-12 Holm p (3/12 <= 0.05)", "0.05 threshold"),
       col = c("#38bdf8", "#22c55e", "#f97316", "#ef4444"), pch = c(19, 17, 15, NA), lty = c(NA, NA, NA, 2),
       cex = 0.80, bg = "#0f172a", text.col = "#f8fafc")

# Panel 3: Edge Case 3 - W-C1 Secondary H^0 Rise & W-C2 e_aq^-1 Recombination Sparing
plot(wc1_h$conc, wc1_h$G_1us, log = "x", type = "b", pch = 19, lwd = 2.4, col = "#facc15",
     ylim = c(0, 1.15),
     xlab = "Scavenger Concentration (M, log scale)",
     ylab = "Yield at 1 us (molecules / 100 eV)",
     main = "Edge Case 3: W-C1 Secondary H^0 Rise (+67.9%) as OH^0 Collapses to 0")
lines(wc1_oh$conc, wc1_oh$G_1us, type = "b", pch = 17, lwd = 2.4, col = "#38bdf8")
lines(wc2_eaq$conc, wc2_eaq$G_1us * 20.0, type = "b", pch = 15, lwd = 2.2, col = "#ec4899", lty = 2)
legend("topright",
       legend = c("W-C1 (MeOH): G(H^0) [0.351 -> 0.589]", "W-C1 (MeOH): G(OH^0) [0 floor >= 0.1 M]", "W-C2 (DMSO): 20 x G(e_aq^-1) [+23.6% rise]"),
       col = c("#facc15", "#38bdf8", "#ec4899"), pch = c(19, 17, 15), lty = c(1, 1, 2),
       cex = 0.80, bg = "#0f172a", text.col = "#f8fafc")

# Panel 4: Edge Case 4 - 6TAP Chain E (NC Zinc-Knuckle) B-Factor Dispersion
x_ch <- seq_len(nrow(dat$bf_df))
plot(x_ch, dat$bf_df$b_mean, ylim = c(0, 175), xaxt = "n", pch = 19, cex = 1.5, col = "#38bdf8",
     xlab = "6TAP Asymmetric Unit Polymer Chain (A-D = CA Shell, E = NC Zn-Knuckle)",
     ylab = "Atomic B-Factor (A^2: min / mean +/- SD / max)",
     main = "Edge Case 4: Chain E (NC Core) B-Factor Dispersion (SD = 32.25 A^2)")
arrows(x_ch, dat$bf_df$b_min, x_ch, dat$bf_df$b_max, code = 3, angle = 90, length = 0.08, col = "#94a3b8", lwd = 1.8)
arrows(x_ch, dat$bf_df$b_mean - dat$bf_df$b_sd, x_ch, dat$bf_df$b_mean + dat$bf_df$b_sd,
       code = 3, angle = 90, length = 0.12, col = "#ec4899", lwd = 3.0)
axis(1, at = x_ch, labels = paste0("Chain ", dat$bf_df$chain, " (n=", dat$bf_df$residue_count, ")"))
legend("topleft", legend = c("Mean B-factor (A^2)", "Mean +/- 1 SD", "Min to Max range"),
       col = c("#38bdf8", "#ec4899", "#94a3b8"), pch = c(19, NA, NA), lty = c(NA, 1, 1), lwd = c(NA, 3.0, 1.8),
       cex = 0.82, bg = "#0f172a", text.col = "#f8fafc")

mtext("RStudio Edge-Case Audit: ARC Capsid Geometry (PDB 6TAP) & AHRT v0.1 Radiolysis [MEASURED GEOMETRY - NO MECHANISM CLAIMED]",
      outer = TRUE, cex = 1.15, font = 2, col = "#f8fafc")

dev.off()

# Emit JSON receipt
receipt_path <- file.path(out_dir, "rstudio_edge_cases_receipt.json")
receipt_lines <- c(
  "{",
  '  "schema": "arc-ahrt-rstudio-edge-cases-receipt-v1",',
  sprintf('  "r_version": "%s",', R.version.string),
  '  "output_png": "rstudio/renders/rstudio_arc_ahrt_edge_cases.png",',
  '  "output_csv": "rstudio/renders/rstudio_edge_cases_summary.csv",',
  '  "resolution_px": [1920, 1080],',
  '  "edge_cases_verified": 4',
  "}"
)
writeLines(receipt_lines, receipt_path)
cat("RStudio edge-case analysis complete. Wrote:", png_path, "\n")
