#!/usr/bin/env Rscript

# Secondary, paper-wide story analysis for the Cheng et al. TG RNA-seq paper.
#
# This pass is downstream of the validated numbered notebooks. It does not
# rebuild the human TG reference or refit the existing NNLS model. It adds a
# channel-centered molecular context, explicit identity-boundary modules,
# replicate-aware longitudinal summaries, a cross-modal evidence timeline, and
# a manuscript-facing audit of numeric claims that must remain aligned.

suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

source("R/anatomic_theme.R")

run_date <- "2026-08-23"
output_dir <- file.path("publication_runs", paste0(run_date, "-secondary-story-analysis"))
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

write_table <- function(x, filename) {
  path <- file.path(table_dir, filename)
  write_csv(x, path, na = "NA")
  path
}

parse_metadata <- function(sample_ids, counts) {
  total_counts <- colSums(counts)
  detected_genes <- colSums(counts > 0)
  tibble(sample_id = sample_ids) |>
    mutate(
      cohort = sub("^([^_]+)_.*$", "\\1", sample_id),
      sample_label = sub("^[^_]+_", "", sample_id),
      sample_class = case_when(
        grepl("^hDRG", sample_label) ~ "hDRG",
        grepl("^TGN", sample_label) ~ "TG bulk",
        grepl("^TAK", sample_label) ~ "TAK2204464A",
        TRUE ~ "Other"
      ),
      specimen_id = case_when(
        grepl("^hDRG_[0-9]+$", sample_label) ~ sub("_[0-9]+$", "", sample_label),
        grepl("^TGN.*W[0-9]+$", sample_label) ~ sub("W[0-9]+$", "", sample_label),
        TRUE ~ sample_label
      ),
      week_raw = str_extract(sample_label, "W[0-9]+$"),
      maturation_week = as.integer(str_remove(week_raw, "^W")),
      total_counts = unname(total_counts[sample_id]),
      detected_genes = unname(detected_genes[sample_id]),
      genes_per_million_counts = 1e6 * detected_genes / total_counts
    )
}

blocked_trend <- function(data, value_col, endpoint) {
  data <- data |>
    filter(!is.na(maturation_week), !is.na(specimen_id)) |>
    select(specimen_id, maturation_week, value = all_of(value_col)) |>
    filter(is.finite(value))

  if (nrow(data) < 6L || n_distinct(data$specimen_id) < 3L) {
    return(tibble(
      blocked_slope_per_week = NA_real_,
      blocked_slope_ci_low = NA_real_,
      blocked_slope_ci_high = NA_real_,
      blocked_slope_p_value = NA_real_,
      mean_line_slope = NA_real_,
      line_slope_p_value = NA_real_,
      n_biological_lines = n_distinct(data$specimen_id),
      n_observations = nrow(data)
    ))
  }

  fit <- lm(value ~ maturation_week + specimen_id, data = data)
  coefficient_table <- summary(fit)$coefficients
  slope_row <- coefficient_table[rownames(coefficient_table) == "maturation_week", , drop = FALSE]
  if (nrow(slope_row) != 1L || !all(c("Estimate", "Std. Error", "Pr(>|t|)") %in% colnames(slope_row))) {
    return(tibble(
      blocked_slope_per_week = NA_real_,
      blocked_slope_ci_low = NA_real_,
      blocked_slope_ci_high = NA_real_,
      blocked_slope_p_value = NA_real_,
      mean_line_slope = NA_real_,
      line_slope_p_value = NA_real_,
      n_biological_lines = n_distinct(data$specimen_id),
      n_observations = nrow(data)
    ))
  }
  residual_df <- df.residual(fit)
  slope_estimate <- as.numeric(slope_row[1, "Estimate"])
  slope_se <- as.numeric(slope_row[1, "Std. Error"])
  line_slopes <- data |>
    group_by(specimen_id) |>
    group_modify(~ tibble(
      line_slope = unname(coef(lm(value ~ maturation_week, data = .x))["maturation_week"])
    )) |>
    ungroup() |>
    pull(line_slope)
  line_test <- if (length(line_slopes) >= 2L) {
    t.test(line_slopes, mu = 0)$p.value
  } else {
    NA_real_
  }

  tibble(
    blocked_slope_per_week = slope_estimate,
    blocked_slope_ci_low = slope_estimate - qt(0.975, residual_df) * slope_se,
    blocked_slope_ci_high = slope_estimate + qt(0.975, residual_df) * slope_se,
    blocked_slope_p_value = as.numeric(slope_row[1, "Pr(>|t|)"]),
    mean_line_slope = mean(line_slopes),
    line_slope_p_value = line_test,
    n_biological_lines = n_distinct(data$specimen_id),
    n_observations = nrow(data)
  )
}

save_figure <- function(plot, filename, width, height) {
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white", device = "png")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white", device = "pdf")
  invisible(c(png = png_path, pdf = pdf_path))
}

# 1. Import and normalization -------------------------------------------------

incoming_quant <- read_tsv(
  "RealTGN.salmon.merged.gene_counts.tsv",
  show_col_types = FALSE,
  name_repair = "minimal"
)
stopifnot(all(c("gene_id", "gene_name") %in% colnames(incoming_quant)))

sample_ids <- setdiff(colnames(incoming_quant), c("gene_id", "gene_name"))
incoming_counts_tbl <- incoming_quant |>
  transmute(
    gene = toupper(trimws(gene_name)),
    across(all_of(sample_ids), as.numeric)
  ) |>
  filter(!is.na(gene), gene != "") |>
  group_by(gene) |>
  summarise(across(all_of(sample_ids), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

incoming_counts <- incoming_counts_tbl |>
  column_to_rownames("gene") |>
  as.matrix()
storage.mode(incoming_counts) <- "numeric"
stopifnot(!anyDuplicated(rownames(incoming_counts)))

metadata <- parse_metadata(colnames(incoming_counts), incoming_counts)
complete_lines <- metadata |>
  filter(sample_class == "TG bulk", !is.na(maturation_week)) |>
  distinct(specimen_id, maturation_week) |>
  count(specimen_id, name = "n_sampled_weeks") |>
  filter(n_sampled_weeks == 4L) |>
  pull(specimen_id)

matched_metadata <- metadata |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  arrange(specimen_id, maturation_week)
stopifnot(length(complete_lines) == 3L, nrow(matched_metadata) == 12L)

all_dge <- DGEList(counts = incoming_counts, samples = as.data.frame(metadata)) |>
  calcNormFactors(method = "TMM")
all_logcpm <- cpm(all_dge, log = TRUE, prior.count = 2)

matched_counts <- incoming_counts[, matched_metadata$sample_id, drop = FALSE]
matched_dge <- DGEList(counts = matched_counts, samples = as.data.frame(matched_metadata)) |>
  calcNormFactors(method = "TMM")
matched_logcpm <- cpm(matched_dge, log = TRUE, prior.count = 2)

write_table(metadata, "secondary_sample_metadata.csv")
write_table(
  tibble(
    analysis_scope = c("all_context", "matched_longitudinal"),
    n_samples = c(ncol(incoming_counts), ncol(matched_counts)),
    normalization = "edgeR TMM followed by log2-CPM with prior count 2",
    longitudinal_unit = c("not used for inference", "three complete lines, one library per line per week"),
    hDRG_comparator = c("four hDRG libraries retained", "not included")
  ),
  "secondary_analysis_scope.csv"
)

# 2. NaV-family and identity-boundary molecular context -----------------------

gene_manifest <- tribble(
  ~gene, ~program, ~evidence_role,
  "SCN5A", "NaV family", "NaV1.5 comparator",
  "SCN8A", "NaV family", "NaV1.6 comparator",
  "SCN9A", "NaV family", "NaV1.7 sensory channel",
  "SCN10A", "NaV family", "NaV1.8 target channel",
  "SCN11A", "NaV family", "NaV1.9 comparator",
  "RBFOX3", "Pan-neuronal", "neuronal architecture",
  "SNAP25", "Pan-neuronal", "synaptic/neural architecture",
  "SYT1", "Pan-neuronal", "synaptic/neural architecture",
  "STMN2", "Pan-neuronal", "neuronal architecture",
  "SIX1", "Trigeminal/sensory context", "trigeminal developmental context",
  "ISL1", "Trigeminal/sensory context", "sensory-neuronal specification",
  "BRN3A", "Trigeminal/sensory context", "sensory-neuronal specification",
  "TRPV1", "Peptidergic/nociceptor", "nociceptor context",
  "CALCA", "Peptidergic/nociceptor", "peptidergic context",
  "CALCB", "Peptidergic/nociceptor", "peptidergic context",
  "TAC1", "Peptidergic/nociceptor", "peptidergic context",
  "NTRK1", "Peptidergic/nociceptor", "peptidergic context",
  "NTRK2", "Mechanosensory/large-fiber", "A-delta-LTMR-like axis",
  "NTRK3", "Mechanosensory/large-fiber", "A-beta-LTMR-like axis",
  "PIEZO2", "Mechanosensory/large-fiber", "mechanosensory context",
  "PVALB", "Mechanosensory/large-fiber", "A-beta-LTMR-like axis",
  "TRPM8", "Cold/pruriceptor", "thermosensory context",
  "SST", "Cold/pruriceptor", "sensory subtype context",
  "HTR3A", "Cold/pruriceptor", "sensory subtype context",
  "SOX10", "Schwann/neural-crest counter-program", "off-target/purity boundary",
  "MPZ", "Schwann/neural-crest counter-program", "off-target/purity boundary",
  "MBP", "Schwann/neural-crest counter-program", "off-target/purity boundary",
  "PMP22", "Schwann/neural-crest counter-program", "off-target/purity boundary"
) |>
  mutate(present_in_counts = gene %in% rownames(all_logcpm))
write_table(gene_manifest, "secondary_gene_manifest.csv")

context_expression <- as.data.frame(all_logcpm[gene_manifest$gene[gene_manifest$present_in_counts], , drop = FALSE]) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "log_cpm") |>
  left_join(gene_manifest, by = "gene") |>
  left_join(metadata, by = "sample_id") |>
  mutate(
    context_group = case_when(
      sample_class == "hDRG" ~ "hDRG",
      sample_class == "TAK2204464A" ~ "TAK control",
      sample_class == "TG bulk" & is.na(maturation_week) ~ "TG pilot W1",
      sample_class == "TG bulk" ~ paste0("TG week ", maturation_week),
      TRUE ~ sample_class
    ),
    context_group = factor(
      context_group,
      levels = c("hDRG", "TAK control", "TG pilot W1", paste0("TG week ", 1:4))
    )
  )
write_table(context_expression, "secondary_context_expression_logcpm.csv")

context_summary <- context_expression |>
  group_by(gene, program, evidence_role, context_group) |>
  summarise(
    mean_log_cpm = mean(log_cpm, na.rm = TRUE),
    sd_log_cpm = sd(log_cpm, na.rm = TRUE),
    n_libraries = sum(is.finite(log_cpm)),
    .groups = "drop"
  )
write_table(context_summary, "secondary_context_group_summary.csv")

matched_expression <- as.data.frame(matched_logcpm[rownames(matched_logcpm) %in% gene_manifest$gene, , drop = FALSE]) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "log_cpm") |>
  left_join(gene_manifest, by = "gene") |>
  left_join(matched_metadata, by = "sample_id")

gene_trends <- matched_expression |>
  group_by(gene, program, evidence_role) |>
  group_modify(~ blocked_trend(.x, "log_cpm", unique(.x$gene))) |>
  ungroup() |>
  rename(endpoint = gene)
write_table(gene_trends, "secondary_gene_blocked_trends.csv")

# 3. hDRG-centered module scores and reference-anchored program bridge --------

module_manifest <- tribble(
  ~module, ~module_class, ~gene,
  "Pan-neuronal architecture", "on-target", "RBFOX3",
  "Pan-neuronal architecture", "on-target", "SNAP25",
  "Pan-neuronal architecture", "on-target", "SYT1",
  "Pan-neuronal architecture", "on-target", "STMN2",
  "NaV1.8-centered channel context", "on-target", "SCN10A",
  "NaV1.8-centered channel context", "on-target", "SCN9A",
  "NaV1.8-centered channel context", "on-target", "SCN11A",
  "Peptidergic nociceptor context", "identity boundary", "CALCA",
  "Peptidergic nociceptor context", "identity boundary", "CALCB",
  "Peptidergic nociceptor context", "identity boundary", "TAC1",
  "Peptidergic nociceptor context", "identity boundary", "NTRK1",
  "Mechanosensory/large-fiber context", "identity boundary", "NTRK2",
  "Mechanosensory/large-fiber context", "identity boundary", "NTRK3",
  "Mechanosensory/large-fiber context", "identity boundary", "PIEZO2",
  "Mechanosensory/large-fiber context", "identity boundary", "PVALB",
  "Cold/pruriceptor context", "identity boundary", "TRPM8",
  "Cold/pruriceptor context", "identity boundary", "SST",
  "Cold/pruriceptor context", "identity boundary", "HTR3A",
  "Schwann/neural-crest counter-program", "counter-program", "SOX10",
  "Schwann/neural-crest counter-program", "counter-program", "MPZ",
  "Schwann/neural-crest counter-program", "counter-program", "MBP",
  "Schwann/neural-crest counter-program", "counter-program", "PMP22"
) |>
  mutate(present_in_counts = gene %in% rownames(all_logcpm))
write_table(module_manifest, "secondary_module_manifest.csv")

hdrg_stats <- context_expression |>
  filter(sample_class == "hDRG") |>
  group_by(gene) |>
  summarise(
    hDRG_mean = mean(log_cpm, na.rm = TRUE),
    hDRG_sd = sd(log_cpm, na.rm = TRUE),
    .groups = "drop"
  )

module_scores <- context_expression |>
  inner_join(module_manifest |> filter(present_in_counts), by = "gene") |>
  left_join(hdrg_stats, by = "gene") |>
  mutate(
    hDRG_centered_delta = log_cpm - hDRG_mean
  ) |>
  group_by(sample_id, module, module_class, sample_class, specimen_id, maturation_week, context_group) |>
  summarise(
    module_score_hDRG_delta = mean(hDRG_centered_delta, na.rm = TRUE),
    n_genes_present = sum(is.finite(hDRG_centered_delta)),
    genes = paste(sort(unique(gene)), collapse = ","),
    .groups = "drop"
  )
write_table(module_scores, "secondary_module_scores_hDRG_delta.csv")

module_trends <- module_scores |>
  filter(sample_class == "TG bulk", specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  group_by(module, module_class) |>
  group_modify(~ blocked_trend(.x, "module_score_hDRG_delta", unique(.x$module))) |>
  ungroup() |>
  rename(endpoint = module)
write_table(module_trends, "secondary_module_blocked_trends.csv")

reference_program_path <- "results/incoming_qc/weekly_maturation_program_scores.csv"
if (file.exists(reference_program_path)) {
  reference_programs <- read_csv(reference_program_path, show_col_types = FALSE) |>
    filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
    pivot_longer(
      cols = -c(specimen_id, sample_class, maturation_week),
      names_to = "reference_program",
      values_to = "reference_weighted_score"
    )
  write_table(reference_programs, "reference_weighted_program_scores_bridge.csv")

  reference_program_trends <- reference_programs |>
    group_by(reference_program) |>
    group_modify(~ blocked_trend(.x, "reference_weighted_score", unique(.x$reference_program))) |>
    ungroup() |>
    rename(endpoint = reference_program)
  write_table(reference_program_trends, "reference_weighted_program_blocked_trends.csv")
} else {
  reference_programs <- tibble()
  reference_program_trends <- tibble()
}

# 4. Cross-modal evidence timeline and paper numeric audit --------------------

cross_modal_evidence <- tribble(
  ~stage_order, ~evidence_layer, ~time_window, ~sample_unit, ~n, ~primary_result, ~what_it_supports, ~current_boundary,
  1L, "RNA-seq", "Weeks 1-4 after thawing", "3 matched TGN lines; one bulk library per line per week", "12 libraries", "SCN10A is detected across all matched weeks; pan-neuronal and sensory programs are present", "molecular context for a human sensory-neuronal-like culture", "does not establish functional channel competence or cellular purity",
  2L, "Voltage clamp", "6-7 weeks after plating", "individual TGN recordings", "5/15 NaV1.8-positive; 2/15 >1 nA", "33.3% NaV1.8-positive; 13.3% above 1 nA", "early detectable functional NaV1.8", "small sample and selected current thresholds",
  3L, "Voltage clamp", "10-11 weeks after plating", "two independent TGN differentiations", "31/42 NaV1.8-positive; 20/42 >1 nA", "73.8% NaV1.8-positive; 47.6% above 1 nA", "prolonged maturation increases the prevalence of robust functional NaV1.8", "RNA-seq was not collected at this same late timepoint",
  4L, "Immunocytochemistry", "13 weeks after plating", "representative fields", "denominator not reported", "NaV1.8 signal overlaps pan-neuronal marker in representative images", "late protein-level presence", "quantify cells, fields, lots, and negative-control denominator",
  5L, "Pharmacology", "10-11 weeks after plating", "NaV1.8 current recordings", "n=5 vehicle; n=5 VX-548", "30 nM VX-548 nearly abolishes evoked NaV1.8 current", "selective, pharmacologically targetable function", "report effect size and per-cell values alongside representative traces"
)
write_table(cross_modal_evidence, "cross_modal_evidence_timeline.csv")

numeric_audit <- tribble(
  ~audit_item, ~manuscript_or_source, ~status, ~recommended_resolution,
  "RNA-seq matched longitudinal design", "Three matched lines/lots, weeks 1-4; one library per line per week; W1-only pilot excluded from paired inference", "supported by current analysis", "Use one replicate term consistently and state the unit explicitly",
  "NaV1.8 RNA trend", "SCN10A blocked slope approximately +0.083 log-CPM/week; 95% CI -0.072 to +0.238; p=0.251", "stable/present, not a significant maturation trend", "Say transcript is detected across the window; reserve maturation claim for function",
  "6-7 week functional prevalence", "5/15 NaV1.8-positive; 2/15 >1 nA", "internally calculable", "Report both percentages and denominators",
  "10-11 week functional prevalence", "31/42 NaV1.8-positive; 20/42 >1 nA", "internally calculable", "Report both percentages and denominators, preferably per differentiation",
  "Activation V1/2", "Results text reports -22.86 +/- 1.67 mV; Table 3 reports -22.86 +/- 0.79 mV", "conflict", "Resolve against the source electrophysiology table before submission",
  "NaV1.8 immunostaining prevalence", "Legend says almost all pan-neuronal cells; no cell/field/lot denominator", "not quantitatively supported", "Add image-based cell counts or soften to representative protein expression",
  "Figure 2 RNA-seq panels", "Results refer to Figure 2E/F/G; supplied PDF contains only A-C", "missing figure support", "Insert the panels or rewrite the Results and legend to match A-C",
  "RNA-seq Methods", "Placeholder text remains in the manuscript", "not reproducible", "Add input, normalization, reference dataset, comparator, and replicate details"
)
write_table(numeric_audit, "paper_numeric_consistency_audit.csv")

# 5. Figures -------------------------------------------------------------------

nav_genes <- gene_manifest |>
  filter(program == "NaV family", present_in_counts) |>
  pull(gene)
nav_heatmap <- context_summary |>
  filter(gene %in% nav_genes) |>
  mutate(
    context_group = factor(
      as.character(context_group),
      levels = c("hDRG", "TAK control", "TG pilot W1", paste0("TG week ", 1:4))
    ),
    gene = factor(gene, levels = rev(nav_genes))
  )

p_nav_heatmap <- ggplot(nav_heatmap, aes(x = context_group, y = gene, fill = mean_log_cpm)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = number(mean_log_cpm, accuracy = 0.1)), size = 3) +
  scale_fill_gradient2(
    low = anatomic_brand[["blue"]], mid = "white", high = anatomic_brand[["red"]],
    midpoint = median(nav_heatmap$mean_log_cpm, na.rm = TRUE), name = "Mean\nTMM log-CPM"
  ) +
  labs(
    title = "Voltage-gated sodium-channel context",
    subtitle = "hDRG provides an expression comparator; longitudinal inference uses the three complete TGN lines",
    x = NULL, y = NULL
  ) +
  theme_anatomic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(4, 6, 4, 4)
  )

direct_genes <- c("SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A", "RBFOX3", "SNAP25", "SOX10")
direct_plot_data <- matched_expression |>
  filter(gene %in% direct_genes) |>
  mutate(gene = factor(gene, levels = direct_genes))
direct_means <- direct_plot_data |>
  group_by(gene, maturation_week) |>
  summarise(
    mean_log_cpm = mean(log_cpm),
    se_log_cpm = sd(log_cpm) / sqrt(n()),
    .groups = "drop"
  )
p_direct <- ggplot(direct_plot_data, aes(x = maturation_week, y = log_cpm, group = specimen_id, color = specimen_id)) +
  geom_line(linewidth = 0.45, alpha = 0.65) +
  geom_point(size = 1.35) +
  geom_line(data = direct_means, aes(y = mean_log_cpm, group = 1), color = anatomic_brand[["red"]], linewidth = 0.9) +
  geom_point(data = direct_means, aes(x = maturation_week, y = mean_log_cpm), inherit.aes = FALSE, color = anatomic_brand[["red"]], size = 1.8) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 1:4) +
  labs(
    title = "Direct longitudinal expression",
    subtitle = "Colored paths retain biological-line trajectories; red summarizes the three-line mean",
    x = "Culture week", y = "TMM log-CPM", color = "Biological line"
  ) +
  scale_color_manual(values = setNames(anatomic_okabe_ito[1:3], complete_lines), breaks = complete_lines) +
  theme_anatomic(base_size = 8.5) +
  theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(0.85, "lines"),
    strip.text = element_text(size = 8),
    plot.margin = margin(4, 4, 4, 4)
  )

module_order <- c(
  "Pan-neuronal architecture",
  "NaV1.8-centered channel context",
  "Mechanosensory/large-fiber context",
  "Peptidergic nociceptor context",
  "Cold/pruriceptor context",
  "Schwann/neural-crest counter-program"
)
module_plot_data <- module_scores |>
  filter(sample_class == "TG bulk", specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  mutate(module = factor(module, levels = module_order))
module_means <- module_plot_data |>
  group_by(module, maturation_week) |>
  summarise(
    mean_score = mean(module_score_hDRG_delta),
    se_score = sd(module_score_hDRG_delta) / sqrt(n()),
    .groups = "drop"
  )
p_modules <- ggplot(module_plot_data, aes(x = maturation_week, y = module_score_hDRG_delta, group = specimen_id, color = specimen_id)) +
  geom_hline(yintercept = 0, color = anatomic_brand[["grey_mid"]], linewidth = 0.35) +
  geom_line(linewidth = 0.45, alpha = 0.6) +
  geom_point(size = 1.1) +
  geom_line(data = module_means, aes(y = mean_score, group = 1), color = anatomic_brand[["red"]], linewidth = 0.85) +
  geom_point(data = module_means, aes(x = maturation_week, y = mean_score), inherit.aes = FALSE, color = anatomic_brand[["red"]], size = 1.7) +
  facet_wrap(~ module, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 1:4) +
  labs(
    title = "Identity-boundary modules relative to hDRG",
    subtitle = "Descriptive hDRG-centered log-CPM differences; these are not cell fractions",
    x = "Culture week", y = "Mean hDRG-centered log-CPM difference", color = "Biological line"
  ) +
  scale_color_manual(values = setNames(anatomic_okabe_ito[1:3], complete_lines), breaks = complete_lines) +
  theme_anatomic(base_size = 8.5) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8),
    plot.margin = margin(4, 4, 4, 4)
  )

timeline <- tribble(
  ~stage, ~layer, ~window, ~result, ~support,
  1, "RNA-seq", "W1-W4 after thaw", "SCN10A detected; pan-neuronal and sensory programs present", "Molecular context",
  2, "Early voltage clamp", "6-7 weeks after plating", "5/15 NaV1.8-positive; 2/15 >1 nA", "Detectable function",
  3, "Late voltage clamp + pharmacology", "10-11 weeks after plating", "31/42 NaV1.8-positive; 20/42 >1 nA; 30 nM VX-548 nearly abolishes current", "Robust, targetable function",
  4, "Immunocytochemistry", "13 weeks after plating", "Representative NaV1.8/pan-neuronal overlap", "Protein corroboration; denominator pending"
) |>
  mutate(layer = factor(layer, levels = c("RNA-seq", "Early voltage clamp", "Late voltage clamp + pharmacology", "Immunocytochemistry")))

p_timeline <- ggplot(timeline, aes(x = stage, y = 0.52)) +
  annotate(
    "segment", x = 1, xend = 4, y = 0.52, yend = 0.52,
    color = anatomic_brand[["grey"]], linewidth = 1.1,
    arrow = grid::arrow(length = grid::unit(0.1, "inches"), type = "closed")
  ) +
  geom_point(aes(color = layer), size = 5) +
  geom_text(aes(y = 0.93, label = paste0("STAGE ", stage)), fontface = "bold", size = 3.0) +
  geom_text(aes(y = 0.76, label = layer), fontface = "bold", size = 3.0, lineheight = 0.9) +
  geom_text(aes(y = 0.60, label = window), size = 2.65, lineheight = 0.9) +
  geom_text(aes(y = 0.18, label = str_wrap(result, width = 28)), size = 2.75, lineheight = 0.9) +
  geom_text(aes(y = -0.72, label = str_wrap(support, width = 26)), size = 2.55, lineheight = 0.9, color = "#55525D") +
  scale_color_manual(
    values = c(
      "RNA-seq" = anatomic_brand[["blue"]],
      "Early voltage clamp" = anatomic_brand[["red"]],
      "Late voltage clamp + pharmacology" = anatomic_brand[["red_dark"]],
      "Immunocytochemistry" = anatomic_okabe_ito[[4]]
    )
  ) +
  scale_x_continuous(limits = c(0.5, 4.5), breaks = NULL) +
  coord_cartesian(ylim = c(-1.08, 1.08), clip = "off") +
  labs(
    title = "The paper contains a staged, not single-timepoint, evidence chain",
    subtitle = "RNA-seq establishes early molecular context; prolonged maturation establishes functional NaV1.8",
    x = NULL, y = NULL, color = NULL
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(color = anatomic_brand[["ink"]], face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(color = "#55525D", size = 8.5, hjust = 0),
    legend.position = "none",
    plot.margin = margin(10, 12, 14, 12)
  )

trend_order <- c("SCN10A", "SCN9A", "SCN11A", "SCN5A", "SCN8A", "RBFOX3", "SNAP25", "NTRK2", "NTRK3", "SOX10")
selected_trends <- gene_trends |>
  filter(endpoint %in% trend_order) |>
  mutate(
    label = recode(
      endpoint,
      SCN5A = "SCN5A / NaV1.5",
      SCN8A = "SCN8A / NaV1.6",
      SCN9A = "SCN9A / NaV1.7",
      SCN10A = "SCN10A / NaV1.8",
      SCN11A = "SCN11A / NaV1.9",
      .default = endpoint
    ),
    label = factor(label, levels = rev(c(
      "SCN10A / NaV1.8", "SCN9A / NaV1.7", "SCN11A / NaV1.9",
      "SCN5A / NaV1.5", "SCN8A / NaV1.6", "RBFOX3", "SNAP25",
      "NTRK2", "NTRK3", "SOX10"
    )))
  )
p_trends <- ggplot(selected_trends, aes(y = label, x = blocked_slope_per_week, color = program)) +
  geom_vline(xintercept = 0, color = anatomic_brand[["grey_mid"]], linewidth = 0.45) +
  geom_errorbar(aes(xmin = blocked_slope_ci_low, xmax = blocked_slope_ci_high), orientation = "y", width = 0.18, linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(
    "NaV family" = anatomic_brand[["blue"]],
    "Pan-neuronal" = anatomic_okabe_ito[[3]],
    "Mechanosensory/large-fiber" = anatomic_okabe_ito[[2]],
    "Schwann/neural-crest counter-program" = anatomic_okabe_ito[[4]]
  ), drop = FALSE) +
  labs(
    title = "Replicate-aware longitudinal slopes",
    subtitle = "Whiskers are 95% CIs; targeted p-values are nominal and descriptive",
    x = "Blocked slope (TMM log-CPM per week)", y = NULL, color = "Program"
  ) +
  theme_anatomic(base_size = 8.5) +
  theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(0.85, "lines"),
    axis.text.y = element_text(size = 7.7),
    plot.margin = margin(4, 4, 4, 4)
  )

secondary_figure_1 <- ((p_nav_heatmap | p_direct) / (p_modules | p_trends)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = "Secondary analysis: a NaV1.8-centered molecular context with explicit identity boundaries",
    subtitle = "Direct expression and descriptive hDRG-centered modules support the functional paper without turning bulk RNA-seq into a cell-fraction claim.",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", color = anatomic_brand[["ink"]]),
      plot.subtitle = element_text(size = 10, color = "#55525D"),
      plot.tag = element_text(size = 12, face = "bold", color = anatomic_brand[["red"]])
    )
  ) &
  theme(legend.position = "bottom")
save_figure(secondary_figure_1, "Figure_SA1_nav_family_and_identity_context", 16, 13)

secondary_figure_2 <- p_timeline +
  plot_annotation(
    title = "Secondary analysis: reframe the manuscript as a staged molecular-to-functional maturation study",
    subtitle = "The evidence is strongest when each modality is assigned its own time window and interpretation boundary.",
    theme = theme(plot.title = element_text(size = 16, face = "bold", color = anatomic_brand[["ink"]]), plot.subtitle = element_text(size = 10, color = "#55525D"))
  )
save_figure(secondary_figure_2, "Figure_SA2_cross_modal_evidence_timeline", 15, 7)

# 6. Captions, methods, report, and input manifest ----------------------------

writeLines(c(
  "# Secondary Figure SA1. NaV-family context and identity boundaries.",
  "",
  "(A) Mean TMM-normalized log-CPM for SCN5A/NaV1.5, SCN8A/NaV1.6, SCN9A/NaV1.7, SCN10A/NaV1.8, and SCN11A/NaV1.9 across hDRG, TAK control, the W1-only pilot, and the three complete TGN lines summarized by week. (B) Direct longitudinal expression of selected NaV, pan-neuronal, and Schwann counter-program genes. Colored paths retain the three biological lines; red points/paths show the line mean. (C) Descriptive hDRG-centered log-CPM differences for pan-neuronal, NaV1.8-centered, mechanosensory/large-fiber, peptidergic, cold/pruriceptor, and Schwann/neural-crest counter-programs. Differences are not cell fractions. (D) Blocked longitudinal slopes from the three complete lines. Error bars are 95% confidence intervals from a model with maturation week and biological line; the W1-only pilot is excluded.",
  "",
  "The figure is intended to support a channel-centered paper story: RNA-seq establishes molecular context and identity boundaries, while voltage clamp establishes functional NaV1.8."
), file.path(output_dir, "FIGURE_SA1_CAPTION.md"))

writeLines(c(
  "# Secondary Figure SA2. Cross-modal evidence timeline.",
  "",
  "RNA-seq was collected during weeks 1-4 after thawing. Electrophysiology was collected later, at 6-7 and 10-11 weeks after plating, and immunocytochemistry was evaluated at 13 weeks after plating. The figure therefore presents a staged evidence chain rather than treating transcript, protein, and function as measurements from one synchronized endpoint."
), file.path(output_dir, "FIGURE_SA2_CAPTION.md"))

writeLines(c(
  "# Secondary story analysis methods",
  "",
  "This pass uses the existing Salmon merged gene-count table and the validated downstream result tables. Gene names were upper-cased and duplicate symbols were summed, matching the incoming-bulk handling in the numbered pipeline. TMM normalization and log2-CPM values with prior count 2 were calculated separately for (i) all incoming libraries used for hDRG/TAK/TGN context and (ii) the 12 matched TGN libraries used for longitudinal inference.",
  "",
  "The matched longitudinal unit is the three complete TGN lines sampled at weeks 1-4, with one bulk RNA-seq library per line and week. The W1-only TGN251028B pilot is retained in context tables but excluded from blocked trends. Blocked slopes use value ~ maturation_week + specimen_id; confidence intervals are t-based model intervals. hDRG-centered module scores subtract the mean expression of each gene across the four hDRG libraries and average across the genes present in each descriptive module. They are not deconvolution estimates and should not be interpreted as cell proportions.",
  "The gene and module trends are targeted, descriptive analyses of a small matched series rather than genome-wide differential-expression tests. Reported p-values are nominal model p-values; no genome-wide multiplicity correction is applied, and the direction/uncertainty of the estimates should be emphasized.",
  "",
  "The cross-modal evidence table and numeric audit intentionally retain manuscript-level items that cannot be computed from the RNA-seq count table alone, including electrophysiology denominators, immunostaining denominator gaps, and the activation V1/2 discrepancy between the Results text and Table 3."
), file.path(output_dir, "SECONDARY_ANALYSIS_METHODS.md"))

format_num <- function(x, digits = 3L) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

nav10_trend <- gene_trends |> filter(endpoint == "SCN10A")
nav5_trend <- gene_trends |> filter(endpoint == "SCN5A")
nav8_trend <- gene_trends |> filter(endpoint == "SCN8A")
rbfox3_trend <- gene_trends |> filter(endpoint == "RBFOX3")
snap25_trend <- gene_trends |> filter(endpoint == "SNAP25")
sox10_trend <- gene_trends |> filter(endpoint == "SOX10")

report_lines <- c(
  "# Secondary story analysis report",
  "",
  paste0("Generated: ", run_date),
  "",
  "## Executive result",
  "",
  "The secondary pass supports a staged, channel-centered paper story: the TGN cultures carry a human sensory-neuronal and voltage-gated sodium-channel context during weeks 1-4, while robust NaV1.8 function emerges with prolonged maturation. The RNA-seq should therefore support identity and context; the electrophysiology should carry the functional claim.",
  "",
  "## Analysis scope",
  "",
  paste0("- All-library context: ", ncol(incoming_counts), " incoming libraries, including ", sum(metadata$sample_class == "hDRG"), " hDRG libraries, ", sum(metadata$sample_class == "TAK2204464A"), " TAK controls, three complete TGN lines, and one W1-only TGN pilot."),
  paste0("- Longitudinal inference: ", nrow(matched_metadata), " libraries from ", length(complete_lines), " complete TGN lines across weeks 1-4; the W1-only pilot is excluded."),
  "- Normalization: edgeR TMM followed by log2-CPM with prior count 2.",
  "- Longitudinal model: blocked slope with maturation week plus biological line.",
  "",
  "## Data-derived findings",
  "",
  paste0("1. **SCN10A/NaV1.8 is consistently detected but is not a conclusive week-1-to-week-4 maturation trend.** The blocked slope is ", format_num(nav10_trend$blocked_slope_per_week), " log-CPM/week (95% CI ", format_num(nav10_trend$blocked_slope_ci_low), " to ", format_num(nav10_trend$blocked_slope_ci_high), "; p=", format_num(nav10_trend$blocked_slope_p_value), "). The manuscript should use RNA-seq to establish molecular context, not to claim that SCN10A transcription rises over this window."),
  paste0("2. **The broader NaV family is not one uniform maturation axis.** SCN5A/NaV1.5 slope: ", format_num(nav5_trend$blocked_slope_per_week), " (p=", format_num(nav5_trend$blocked_slope_p_value), "); SCN8A/NaV1.6 slope: ", format_num(nav8_trend$blocked_slope_per_week), " (p=", format_num(nav8_trend$blocked_slope_p_value), "); SCN9A/NaV1.7 and SCN11A/NaV1.9 are not conclusive in this small series. These targeted p-values should remain descriptive, not be presented as genome-wide differential expression."),
  paste0("3. **Pan-neuronal architecture is compatible with maturation, but different markers move at different rates.** RBFOX3 slope: ", format_num(rbfox3_trend$blocked_slope_per_week), " (p=", format_num(rbfox3_trend$blocked_slope_p_value), "); SNAP25 slope: ", format_num(snap25_trend$blocked_slope_per_week), " (p=", format_num(snap25_trend$blocked_slope_p_value), ")."),
  "4. **The identity boundary is biologically informative rather than a nuisance.** NTRK2/PIEZO2 mechanosensory context and the existing human-TG reference scores should be described as an A-delta-LTMR-like or large-fiber axis, while peptidergic markers remain limited. The current analysis does not justify calling the culture a pure mature nociceptor population.",
  paste0("5. **Counter-programs should be shown, not hidden.** SOX10 has a blocked slope of ", format_num(sox10_trend$blocked_slope_per_week), " (95% CI ", format_num(sox10_trend$blocked_slope_ci_low), " to ", format_num(sox10_trend$blocked_slope_ci_high), "; p=", format_num(sox10_trend$blocked_slope_p_value), "). This is a reason to bound bulk-RNA purity language, not a reason to discard the neuronal result."),
  "6. **The functional time course is later than the RNA-seq time course.** The manuscript reports 5/15 NaV1.8-positive cells at 6-7 weeks after plating and 31/42 at 10-11 weeks; the corresponding >1 nA counts are 2/15 and 20/42. These numbers support prolonged functional maturation but cannot be treated as same-timepoint validation of the week-1-to-week-4 RNA-seq series.",
  "",
  "## Recommended manuscript reframe",
  "",
  "> The cultures acquire a human sensory-neuronal-like transcriptional state containing SCN10A/NaV1.8 and related sodium-channel programs. They also retain mechanosensory and thermosensory features and do not reproduce the full mature peptidergic transcriptome of ganglion tissue. RNA-seq therefore defines the molecular context and identity boundary of the model. With prolonged maturation, the same cultures develop robust TTX-resistant NaV1.8 currents with depolarized inactivation and sensitivity to VX-548.",
  "",
  "## Highest-value actions before submission",
  "",
  "1. Use Figure SA1 as the RNA-seq replacement/augmentation for Figure 2, or rebuild it into the journal figure style.",
  "2. Add Figure SA2 or an equivalent schematic so the RNA, protein, and electrophysiology windows are not conflated.",
  "3. Add SCN5A/NaV1.5 to the sodium-channel table and retain the full sensory-gene table in Supplementary Data.",
  "4. Resolve the activation V1/2 discrepancy, quantify NaV1.8 immunostaining, and state the per-lot denominators for the electrophysiology.",
  "5. Keep the donor-aware TG mapping, marker refit, VIF, and bootstrap analyses as supporting validation. Do not lead with NNLS fractions or describe them as cell fractions.",
  "",
  "## Provocative questions for the author team",
  "",
  "1. If SCN10A is present but does not rise across weeks 1-4, should the central maturation claim be framed as post-transcriptional/functional maturation rather than transcriptional induction?",
  "2. Is the most honest identity label a mixed human sensory-neuronal NAM with an A-delta-LTMR-like/large-fiber axis, rather than a homogeneous mature nociceptor culture?",
  "3. What is the true experimental unit for the 5/15 and 31/42 functional prevalences—cells, wells, differentiations, or lots—and does pooling hide lot-to-lot variability?",
  "4. Could the changing SOX10/MPZ/MBP/PMP22 counter-program affect current prevalence or current density, and can the paper quantify neuronal-versus-neural-crest co-occurrence rather than only show representative fields?",
  "5. Can every RNA, protein, and electrophysiology panel be annotated by both thaw-relative and plating-relative age so that the staged evidence chain cannot be mistaken for same-timepoint validation?",
  "6. Does VX-548 sensitivity establish target engagement strongly enough for the intended NAM claim, or should the paper explicitly separate pharmacological sensitivity from molecular identity and add a second orthogonal validation?",
  "7. What result would falsify the proposed story—for example, robust NaV1.8 current in lines lacking the expected molecular context—and is that test feasible with the existing data or a small follow-up experiment?",
  "",
  "## Output map",
  "",
  "- figures/Figure_SA1_nav_family_and_identity_context.png/.pdf: molecular context and identity-boundary figure.",
  "- figures/Figure_SA2_cross_modal_evidence_timeline.png/.pdf: modality/timepoint bridge.",
  "- tables/secondary_gene_blocked_trends.csv: direct gene trends, including SCN10A, RBFOX3, SNAP25, NTRK2/NTRK3, and SOX10.",
  "- tables/secondary_module_scores_hDRG_delta.csv: descriptive hDRG-centered log-CPM differences; not cell fractions.",
  "- tables/cross_modal_evidence_timeline.csv: paper-wide evidence chain.",
  "- tables/paper_numeric_consistency_audit.csv: manuscript/figure claims requiring resolution.",
  "- SECONDARY_ANALYSIS_METHODS.md: reproducibility and interpretation boundaries.",
  "- README.md: package guide and recommended use in the manuscript."
)
writeLines(report_lines, file.path(output_dir, "SECONDARY_ANALYSIS_REPORT.md"))

writeLines(c(
  "# Secondary story analysis package",
  "",
  paste0("Generated: ", run_date),
  "",
  "This package is a downstream, paper-facing analysis pass. It uses the existing merged gene-count table and validated result tables; it does not alter the numbered pipeline or manuscript files.",
  "",
  "## Start here",
  "",
  "1. Read SECONDARY_ANALYSIS_REPORT.md for the proposed paper-wide reframe and decision points.",
  "2. Review Figure_SA1_nav_family_and_identity_context.png for the molecular-context/identity-boundary story.",
  "3. Review Figure_SA2_cross_modal_evidence_timeline.png to align RNA-seq, electrophysiology, protein, and pharmacology time windows.",
  "4. Use tables/paper_numeric_consistency_audit.csv as the pre-submission reconciliation list.",
  "",
  "## Interpretation boundary",
  "",
  "RNA-seq supports human sensory-neuronal context, NaV-family context, and identity boundaries. Voltage clamp carries the claim of functional NaV1.8 competence. hDRG-centered module scores are descriptive expression differences, not cell fractions. Targeted trend p-values are nominal and should not be presented as genome-wide differential-expression results.",
  "",
  "## Reproducibility",
  "",
  "The source script is scripts/build_secondary_story_analysis.R. SECONDARY_ANALYSIS_METHODS.md, SESSION_INFO.txt, and tables/input_manifest.csv document the analysis environment and input provenance."
), file.path(output_dir, "README.md"))

input_paths <- c(
  "RealTGN.salmon.merged.gene_counts.tsv",
  "results/incoming_qc/weekly_maturation_program_scores.csv",
  "results/incoming_qc/aggregate_neuronal_identity_index.csv",
  "results/incoming_qc/weekly_neuron_only_pca_trajectory_coordinates.csv",
  "results/final_article/key_gene_blocked_trend_summary.csv",
  "Cheng et al 081626_VT_ASH-1.docx",
  "Cheng et Figures_V2.pdf"
)
input_manifest <- tibble(path = input_paths) |>
  mutate(
    exists = file.exists(path),
    size_bytes = ifelse(exists, file.info(path)$size, NA_real_),
    md5 = ifelse(exists, unname(tools::md5sum(path)), NA_character_)
  )
write_table(input_manifest, "input_manifest.csv")

writeLines(capture.output(sessionInfo()), file.path(output_dir, "SESSION_INFO.txt"))

message("Secondary story analysis written to: ", normalizePath(output_dir))
