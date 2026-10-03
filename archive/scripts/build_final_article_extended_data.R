#!/usr/bin/env Rscript

# Render one Extended Data composite in a clean R process.  This avoids keeping
# the large reference object and accumulated ggplot grobs in the Quarto session.

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[[1]] else ""
valid_targets <- c("analysis-context", "reference-calibration")
if (!target %in% valid_targets) {
  stop("Usage: Rscript scripts/build_final_article_extended_data.R [",
       paste(valid_targets, collapse = "|"), "]")
}

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scales)

source("R/anatomic_theme.R")

output_dir <- "results/final_article"
publication_run <- Sys.getenv("PUBLICATION_RUN", "2026-08-04")
figure_dir <- file.path("publication_runs", publication_run, "final_article_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (target == "analysis-context") {
  sample_metadata <- read_csv("results/incoming_qc/incoming_sample_metadata.csv", show_col_types = FALSE)
  complete_lines <- sample_metadata |>
    filter(sample_class == "TG bulk", !is.na(maturation_week)) |>
    distinct(specimen_id, maturation_week) |>
    count(specimen_id, name = "n_sampled_weeks") |>
    filter(n_sampled_weeks == 4L) |>
    pull(specimen_id)
  stopifnot(length(complete_lines) == 3L)
  matched_metadata <- sample_metadata |>
    filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
    arrange(specimen_id, maturation_week)

  qc_flags <- read_csv("results/incoming_qc/incoming_sample_qc_flags.csv", show_col_types = FALSE)
  program_scores_wide <- read_csv("results/incoming_qc/weekly_maturation_program_scores.csv", show_col_types = FALSE)
  top_reference_matches <- read_csv("results/figures/RealTGN_bulk_top_reference_matches.csv", show_col_types = FALSE)

  p_qc <- ggplot(qc_flags, aes(x = total_counts, y = detected_genes, color = sample_class, shape = qc_status)) +
    geom_vline(xintercept = 5e6, linetype = "dashed", color = anatomic_brand[["red"]]) +
    geom_hline(yintercept = 1e4, linetype = "dashed", color = anatomic_brand[["red"]]) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
    labs(title = "All incoming libraries pass operational QC floors", subtitle = "Dashed lines: 5 million counts and 10,000 detected genes", x = "Total counts", y = "Detected genes", color = "Sample class", shape = "QC status") +
    theme_anatomic(base_size = 9) + theme(legend.position = "bottom")

  program_scores_long <- program_scores_wide |>
    filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
    pivot_longer(cols = -c(specimen_id, sample_class, maturation_week), names_to = "program", values_to = "signature_score")
  p_programs <- ggplot(program_scores_long, aes(x = maturation_week, y = signature_score, group = specimen_id, color = specimen_id)) +
    geom_hline(yintercept = 0, color = "grey80") + geom_line(linewidth = 0.55, alpha = 0.75) +
    geom_point(size = 1.5, alpha = 0.85) + facet_wrap(~ program, scales = "free_y", ncol = 5) +
    scale_x_continuous(breaks = 1:4) +
    labs(title = "Full reference-weighted program readout", subtitle = "All neuronal and off-target programs retained from the originating analysis", x = "Culture week", y = "Weighted signature score", color = "Biological line") +
    theme_anatomic(base_size = 7.5) + theme(legend.position = "bottom")

  bridge_matches <- top_reference_matches |>
    inner_join(matched_metadata |> select(sample_id, specimen_id, maturation_week), by = c("bulk_sample" = "sample_id")) |>
    group_by(bulk_sample) |> slice_max(median_correlation, n = 1, with_ties = FALSE) |> ungroup()
  p_bridge <- ggplot(bridge_matches, aes(x = maturation_week, y = reference_label, color = median_correlation)) +
    geom_point(size = 4, alpha = 0.9) + facet_wrap(~ specimen_id, nrow = 1) +
    scale_x_continuous(breaks = 1:4) + scale_color_viridis_c(option = "C", end = 0.9) +
    labs(title = "Top reference match across culture trajectories", subtitle = "Color = median bulk-to-reference Spearman correlation", x = "Culture week", y = "Top reference identity", color = "Correlation") +
    theme_anatomic(base_size = 8.5) + theme(legend.position = "bottom")

  figure <- (p_qc | p_bridge) / p_programs +
    plot_layout(heights = c(0.8, 1.2)) +
    plot_annotation(title = "Extended analysis context: sequencing quality, correlation-weighted reference bridge, and full programs", subtitle = "Supporting panels retain the breadth of the originating RNA-seq analysis while the main figures prioritize the primary biological story.", tag_levels = "A")
  png_path <- file.path(figure_dir, "Extended_Data_Figure_1_analysis_context.png")
  pdf_path <- file.path(figure_dir, "Extended_Data_Figure_1_analysis_context.pdf")
  ggsave(png_path, figure, width = 12, height = 8, dpi = 200, bg = "white")
  ggsave(pdf_path, figure, width = 12, height = 8, bg = "white", useDingbats = FALSE)
  writeLines(c("# Extended Data Figure 1", "", "Analytical context for the final maturity and composition figures.", "", "A. Library depth and detected-gene complexity for all incoming samples relative to operational QC floors.", "B. Top reference match for each matched TG culture sample; color encodes median bulk-to-reference Spearman correlation.", "C. Full reference-weighted neuronal and non-neuronal program trajectories."), file.path(output_dir, "EXTENDED_DATA_FIGURE_1_CAPTION.md"))
}

if (target == "reference-calibration") {
  library(Seurat)
  consensus_object <- readRDS("results/consensus_identity/GSE197289_human_TG_neurons_consensus.rds")
  reference_umap <- Embeddings(consensus_object, "umap_sct") |>
    as.data.frame() |> rownames_to_column("cell")
  colnames(reference_umap)[2:3] <- c("UMAP_1", "UMAP_2")
  reference_umap <- reference_umap |>
    left_join(consensus_object[[]] |> rownames_to_column("cell"), by = "cell")
  p_reference_umap <- ggplot(reference_umap, aes(x = UMAP_1, y = UMAP_2, color = literature_class)) +
    geom_point(size = 0.55, alpha = 0.6) +
    geom_point(data = filter(reference_umap, consensus_status != "Donor-supported anchor"), shape = 4, size = 1.1, color = anatomic_brand[["ink"]]) +
    coord_equal() +
    guides(color = guide_legend(override.aes = list(size = 3.2, alpha = 1), nrow = 2, byrow = TRUE)) +
    labs(tag = "A", title = "Human TG neuronal reference", subtitle = "Crosses: low-support intersections retained without a neuronal-subtype call", x = "UMAP 1", y = "UMAP 2", color = "Literature-defined class") +
    theme_anatomic(base_size = 9) +
    theme(legend.position = "bottom", legend.key.width = grid::unit(0.45, "cm"), legend.key.height = grid::unit(0.4, "cm"), plot.tag = element_text(face = "bold", color = anatomic_brand[["red"]], size = 15), plot.tag.position = c(0.01, 0.99))

  reference_signed_scores <- read_csv("results/marker_refinement/refined_reference_signed_scores.csv", show_col_types = FALSE)
  score_matrix <- reference_signed_scores |>
    group_by(true_label, target_label) |>
    summarise(mean_signed_score = mean(signed_score), .groups = "drop")
  p_signed_reference <- ggplot(score_matrix, aes(x = target_label, y = true_label, fill = mean_signed_score)) +
    geom_tile() + scale_fill_gradient2(low = "#2166AC", mid = "white", high = anatomic_brand[["red"]], midpoint = 0) +
    labs(tag = "B", title = "Positive-minus-anti-marker identity separation", subtitle = "Mean donor-pseudobulk signed scores; diagonal enrichment supports specific identity recovery", x = "Scored identity", y = "Reference identity", fill = "Mean signed score") +
    theme_anatomic(base_size = 7.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.tag = element_text(face = "bold", color = anatomic_brand[["red"]], size = 15), plot.tag.position = c(0.01, 0.99))

  vif_joint <- read_csv("results/deconv_rerun_expanded/vif_joint_space_comparison.csv", show_col_types = FALSE) |>
    pivot_longer(-reference_identity, names_to = "manifest", values_to = "vif") |>
    mutate(manifest = recode(manifest, vif_original = "Positive markers only", vif_expanded = "Positive + anti-marker matrix"))
  p_vif <- ggplot(vif_joint, aes(x = reorder(reference_identity, vif), y = vif, fill = manifest)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64) +
    geom_hline(yintercept = 10, linetype = "dashed", color = anatomic_brand[["red"]]) + coord_flip() +
    labs(tag = "C", title = "Positive-plus-anti-marker recalibration lowers VIF", subtitle = "Anti-markers are recomputed in donor pseudobulks; dashed line = VIF 10 severe-collinearity cutoff", x = NULL, y = "Variance inflation factor", fill = NULL) +
    theme_anatomic(base_size = 8) + theme(legend.position = "bottom", plot.tag = element_text(face = "bold", color = anatomic_brand[["red"]], size = 15), plot.tag.position = c(0.01, 0.99))
  figure <- p_reference_umap | p_signed_reference | p_vif +
    plot_annotation(title = "Reference calibration supports conservative culture-state interpretation", subtitle = "Identity calls are anchored to donor-supported reference states; positive and anti-marker evidence is recomputed to improve deconvolution conditioning.")
  png_path <- file.path(figure_dir, "Extended_Data_Figure_2_reference_calibration.png")
  pdf_path <- file.path(figure_dir, "Extended_Data_Figure_2_reference_calibration.pdf")
  ggsave(png_path, figure, width = 15, height = 7, dpi = 200, bg = "white")
  ggsave(pdf_path, figure, width = 15, height = 7, bg = "white", useDingbats = FALSE)
  writeLines(c("# Extended Data Figure 2", "", "Reference identity calibration.", "", "A. Human TG neuronal reference UMAP; crosses mark low-support intersections retained without a neuronal-subtype call. Enlarged legend markers denote literature-defined classes.", "B. Mean donor-pseudobulk signed scores calculated from positive-marker evidence minus anti-marker evidence; diagonal enrichment supports specific identity recovery.", "C. VIF after recalculating the marker matrix in donor pseudobulks. The comparison contrasts a simple positive-marker-only matrix with the positive-plus-anti-marker matrix; the dashed line marks the VIF 10 severe-collinearity cutoff."), file.path(output_dir, "EXTENDED_DATA_FIGURE_2_CAPTION.md"))
}

message("Rendered extended data target: ", target)
