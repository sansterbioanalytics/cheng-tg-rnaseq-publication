#!/usr/bin/env Rscript

# Render the two high-resolution main composites in a clean R process.  Quarto
# retains plot grobs for HTML rendering; this script prevents those grobs from
# accumulating in the notebook process while preserving identical PNG/PDF assets.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(edgeR)
  library(scales)
})
source("R/anatomic_theme.R")

publication_run <- Sys.getenv("PUBLICATION_RUN", "2026-08-04")
figure_dir <- file.path("publication_runs", publication_run, "final_article_figures")
panel_tag_theme <- theme(
  plot.tag = element_text(face = "bold", color = anatomic_brand[["red"]], size = 15),
  plot.tag.position = c(0.01, 0.99)
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

summarise_replicates <- function(data, value_column, grouping) {
  data |>
    group_by(across(all_of(grouping))) |>
    summarise(
      n_lines = sum(!is.na(.data[[value_column]])),
      mean_value = mean(.data[[value_column]], na.rm = TRUE),
      sd_value = sd(.data[[value_column]], na.rm = TRUE),
      se_value = sd_value / sqrt(n_lines),
      t_critical = qt(0.975, df = pmax(n_lines - 1, 1)),
      ci_low = mean_value - t_critical * se_value,
      ci_high = mean_value + t_critical * se_value,
      .groups = "drop"
    )
}

plot_replicate_trajectory <- function(raw_data, summary_data, y_column, y_label, title, subtitle = NULL) {
  ggplot(raw_data, aes(x = maturation_week, y = .data[[y_column]], group = specimen_id, color = specimen_id)) +
    geom_line(linewidth = 0.7, alpha = 0.65) +
    geom_point(size = 2, alpha = 0.8) +
    geom_ribbon(data = summary_data, aes(x = maturation_week, ymin = ci_low, ymax = ci_high, group = 1), inherit.aes = FALSE, fill = anatomic_brand[["red"]], alpha = 0.12) +
    geom_line(data = summary_data, aes(x = maturation_week, y = mean_value, group = 1), inherit.aes = FALSE, linewidth = 1.15, color = anatomic_brand[["red"]]) +
    geom_point(data = summary_data, aes(x = maturation_week, y = mean_value), inherit.aes = FALSE, size = 2.7, shape = 21, fill = "white", color = anatomic_brand[["red"]], stroke = 0.9) +
    scale_x_continuous(breaks = 1:4) +
    labs(title = title, subtitle = subtitle, x = "Culture week", y = y_label, color = "Biological line") +
    theme_anatomic() + theme(legend.position = "bottom")
}

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

incoming_counts <- read_tsv("RealTGN.salmon.merged.gene_counts.tsv", show_col_types = FALSE, name_repair = "minimal")
incoming_count_table <- incoming_counts |>
  select(gene_name, all_of(matched_metadata$sample_id)) |>
  mutate(gene_name = toupper(gene_name))
incoming_count_matrix <- as.matrix(incoming_count_table[, matched_metadata$sample_id, drop = FALSE])
rownames(incoming_count_matrix) <- incoming_count_table[["gene_name"]]
storage.mode(incoming_count_matrix) <- "numeric"
incoming_count_matrix <- rowsum(incoming_count_matrix, group = rownames(incoming_count_matrix), reorder = FALSE)
incoming_logcpm <- cpm(calcNormFactors(DGEList(counts = incoming_count_matrix), method = "TMM"), log = TRUE, prior.count = 2)
key_gene_dictionary <- tribble(
  ~gene, ~program,
  "RBFOX3", "Pan-neuronal", "SNAP25", "Pan-neuronal", "SCN10A", "Nav1.8 nociceptor program",
  "CASQ2", "C-LTMR-like", "P2RY1", "C-LTMR-like", "NTRK3", "A-beta-LTMR-like",
  "PVALB", "A-beta-LTMR-like", "NTRK2", "A-delta-LTMR-like", "PIEZO2", "A-delta-LTMR-like",
  "SOX10", "Schwann counter-program"
) |>
  filter(gene %in% rownames(incoming_logcpm))
key_gene_expression <- incoming_logcpm[key_gene_dictionary$gene, , drop = FALSE] |>
  as.data.frame() |> rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "log_cpm") |>
  left_join(key_gene_dictionary, by = "gene") |>
  left_join(matched_metadata |> select(sample_id, specimen_id, maturation_week), by = "sample_id")

maturity_data <- read_csv("results/incoming_qc/aggregate_neuronal_identity_index.csv", show_col_types = FALSE) |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  arrange(specimen_id, maturation_week)
maturity_summary <- summarise_replicates(maturity_data, "neuronal_identity_index", "maturation_week")
pca_data <- read_csv("results/incoming_qc/weekly_neuron_only_pca_trajectory_coordinates.csv", show_col_types = FALSE) |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  arrange(specimen_id, maturation_week)

p_maturity <- plot_replicate_trajectory(maturity_data, maturity_summary, "neuronal_identity_index", "Aggregate neuronal identity index", "Matched-line neuronal identity trajectory") +
  labs(tag = "A") + panel_tag_theme
p_key_genes <- ggplot(key_gene_expression, aes(x = maturation_week, y = log_cpm, group = specimen_id, color = specimen_id)) +
  geom_line(linewidth = 0.55, alpha = 0.7) + geom_point(size = 1.7, alpha = 0.85) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1, color = anatomic_brand[["red"]], inherit.aes = TRUE) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 21, size = 2.3, fill = "white", color = anatomic_brand[["red"]]) +
  facet_wrap(~ program + gene, scales = "free_y", ncol = 3) + scale_x_continuous(breaks = 1:4) +
  labs(tag = "B", title = "Direct expression of neuronal and identity genes", subtitle = "Red = matched-line mean; colored paths = biological lines", x = "Culture week", y = "TMM log-CPM", color = "Biological line") +
  theme_anatomic(base_size = 9) + theme(legend.position = "bottom") + panel_tag_theme
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = specimen_id, group = specimen_id)) +
  geom_path(linewidth = 0.85, arrow = grid::arrow(length = grid::unit(0.09, "inches"))) +
  geom_point(aes(fill = maturation_week), shape = 21, size = 3.4, color = anatomic_brand[["ink"]], stroke = 0.45) +
  geom_point(data = pca_data |> group_by(maturation_week) |> summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop"), aes(x = PC1, y = PC2), inherit.aes = FALSE, shape = 23, size = 4.5, fill = anatomic_brand[["red"]], color = anatomic_brand[["ink"]]) +
  scale_fill_viridis_c(option = "C", end = 0.9, guide = guide_colorbar(title = "Week")) +
  labs(tag = "C", title = "Movement through fixed neuronal reference space", subtitle = "Diamonds = week-specific biological-line centroid", x = "Reference PC1", y = "Reference PC2", color = "Biological line") +
  theme_anatomic(base_size = 9) + theme(legend.position = "bottom") + panel_tag_theme
figure_1 <- p_maturity / (p_key_genes | p_pca) +
  plot_layout(heights = c(0.72, 1.28), widths = c(1.45, 1)) +
  plot_annotation(title = "Developing TG cultures show coordinated neuronal-program and reference-space changes", subtitle = "Three matched biological cell-line trajectories (weeks 1–4); one RNA-seq library per line at each week.")
ggsave(file.path(figure_dir, "Figure_1_neuronal_maturity.png"), figure_1, width = 15, height = 10, dpi =  200, bg = "white")
ggsave(file.path(figure_dir, "Figure_1_neuronal_maturity.pdf"), figure_1, width = 15, height = 10, bg = "white", useDingbats = FALSE)

composition_ci <- read_csv("results/deconv_rerun_expanded/weekly_nnls_deconvolution_with_bootstrap_ci.csv", show_col_types = FALSE)
combined_reference_metadata <- read_csv("results/pseudobulk/GSE197289_TG_combined_reference_metadata.csv", show_col_types = FALSE)
identity_lookup <- combined_reference_metadata |> distinct(reference_identity = reference_label, reference_level)
dominant_identity_labels <- c("cLTMR__SC01", "NF2__SC03", "NF3__SC02")
composition_matched <- composition_ci |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4, manifest == "Signed refined (788 genes)") |>
  left_join(identity_lookup, by = "reference_identity") |>
  mutate(composition_class = case_when(
    reference_identity == "cLTMR__SC01" ~ "C-LTMR-like (CASQ2/P2RY1)",
    reference_identity == "NF2__SC03" ~ "A-beta-LTMR-like (NTRK3/PVALB)",
    reference_identity == "NF3__SC02" ~ "A-delta-LTMR-like (NTRK2/PIEZO2)",
    reference_level == "Neuronal subtype" ~ "Other neuronal identities",
    TRUE ~ "Non-neuronal / off-target"))
composition_class <- composition_matched |>
  group_by(specimen_id, maturation_week, composition_class) |>
  summarise(estimated_fraction = sum(estimated_fraction), .groups = "drop")
composition_summary <- summarise_replicates(composition_class, "estimated_fraction", c("maturation_week", "composition_class"))
signed_scores <- read_csv("results/deconv_rerun_expanded/weekly_signed_refined_identity_scores.csv", show_col_types = FALSE)
signed_class_scores <- signed_scores |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4, refined_identity %in% dominant_identity_labels) |>
  mutate(composition_class = recode(refined_identity, "cLTMR__SC01" = "C-LTMR-like (CASQ2/P2RY1)", "NF2__SC03" = "A-beta-LTMR-like (NTRK3/PVALB)", "NF3__SC02" = "A-delta-LTMR-like (NTRK2/PIEZO2)"))
signed_summary <- summarise_replicates(signed_class_scores, "signed_score", c("maturation_week", "composition_class"))
composition_order <- c("C-LTMR-like (CASQ2/P2RY1)", "A-beta-LTMR-like (NTRK3/PVALB)", "A-delta-LTMR-like (NTRK2/PIEZO2)", "Other neuronal identities", "Non-neuronal / off-target")
p_composition <- ggplot(composition_class |> mutate(composition_class = factor(composition_class, levels = composition_order)), aes(x = maturation_week, y = estimated_fraction, group = specimen_id, color = specimen_id)) +
  geom_line(linewidth = 0.55, alpha = 0.65) + geom_point(size = 1.8, alpha = 0.8) +
  geom_ribbon(data = composition_summary |> mutate(composition_class = factor(composition_class, levels = composition_order)), aes(x = maturation_week, ymin = ci_low, ymax = ci_high, group = 1), inherit.aes = FALSE, fill = anatomic_brand[["red"]], alpha = 0.12) +
  geom_line(data = composition_summary |> mutate(composition_class = factor(composition_class, levels = composition_order)), aes(x = maturation_week, y = mean_value, group = 1), inherit.aes = FALSE, color = anatomic_brand[["red"]], linewidth = 1) +
  facet_wrap(~ composition_class, scales = "free_y", ncol = 3) + scale_x_continuous(breaks = 1:4) + scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(tag = "A", title = "Relative composition across matched biological lines", subtitle = "Red band = 95% CI across lines; colored paths = biological lines", x = "Culture week", y = "Relative signature contribution", color = "Biological line") +
  theme_anatomic(base_size = 9) + theme(legend.position = "bottom") + panel_tag_theme
p_signed <- ggplot(signed_class_scores, aes(x = maturation_week, y = signed_score, group = specimen_id, color = specimen_id)) +
  geom_hline(yintercept = 0, color = "grey75") + geom_line(linewidth = 0.6, alpha = 0.7) + geom_point(size = 1.8, alpha = 0.85) +
  geom_ribbon(data = signed_summary, aes(x = maturation_week, ymin = ci_low, ymax = ci_high, group = 1), inherit.aes = FALSE, fill = anatomic_brand[["red"]], alpha = 0.12) +
  geom_line(data = signed_summary, aes(x = maturation_week, y = mean_value, group = 1), inherit.aes = FALSE, color = anatomic_brand[["red"]], linewidth = 1) +
  facet_wrap(~ composition_class, ncol = 3) + scale_x_continuous(breaks = 1:4) +
  labs(tag = "B", title = "Signed evidence for culture-associated identity axes", subtitle = "Positive markers increase the score; validated anti-markers decrease it", x = "Culture week", y = "Reference-standardized signed score", color = "Biological line") +
  theme_anatomic(base_size = 9) + theme(legend.position = "bottom") + panel_tag_theme
figure_2 <- p_composition / p_signed + plot_layout(heights = c(1.18, 0.82)) +
  plot_annotation(title = "Relative sensory-identity contributions are interpreted alongside direct signed-marker evidence", subtitle = "NNLS reports a relative transcriptional signature contribution, not a literal cell fraction; gene-bootstrap sensitivity is reported separately.")
ggsave(file.path(figure_dir, "Figure_2_relative_composition.png"), figure_2, width = 15, height = 10, dpi =  200, bg = "white")
ggsave(file.path(figure_dir, "Figure_2_relative_composition.pdf"), figure_2, width = 15, height = 10, bg = "white", useDingbats = FALSE)
message("Rendered final main figures")
