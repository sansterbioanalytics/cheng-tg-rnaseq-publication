#!/usr/bin/env Rscript

# Build a single full-page primary figure and a methods-forward Extended Data set.
# All panels are regenerated from exported analysis objects/tables; no composite
# screenshots are embedded. The primary figure follows the computational method
# order: reference labels -> marker programs -> donor-aware signature model ->
# fixed-reference culture mapping -> subtype evidence -> relative contributions.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(grid)
})

source("R/anatomic_theme.R")

run_name <- "2026-08-04-single-primary-figure"
run_dir <- file.path("publication_runs", run_name)
figure_dir <- file.path(run_dir, "figures")
pdf_dir <- file.path("output", "pdf")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

ink <- anatomic_brand[["ink"]]
red <- anatomic_brand[["red"]]
blue <- anatomic_brand[["blue"]]
grey_light <- anatomic_brand[["grey_light"]]

panel_theme <- theme_anatomic(base_size = 7.2) +
  theme(
    plot.title = element_text(size = 8.4, margin = margin(b = 2.5)),
    plot.subtitle = element_text(size = 6.6, lineheight = 0.95, margin = margin(b = 4)),
    plot.caption = element_text(size = 6.1, color = "#55525D", hjust = 0, margin = margin(t = 3)),
    axis.title = element_text(size = 6.8),
    axis.text = element_text(size = 6.2),
    legend.title = element_text(size = 6.6),
    legend.text = element_text(size = 6.1),
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.34, "cm"),
    strip.text = element_text(size = 6.6),
    plot.tag = element_text(face = "bold", size = 12, color = red),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(4, 4, 4, 4)
  )

ed_theme <- theme_anatomic(base_size = 8.2) +
  theme(
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 7.7, lineheight = 0.95),
    plot.caption = element_text(size = 7, color = "#55525D", hjust = 0),
    axis.title = element_text(size = 7.7),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 7.4),
    plot.tag = element_text(face = "bold", size = 14, color = red),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(5, 5, 5, 5)
  )

save_figure <- function(filename, figure, width, height, dpi = 300) {
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  ggsave(png_path, figure, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, figure, width = width, height = height, bg = "white", useDingbats = FALSE, limitsize = FALSE)
  file.copy(pdf_path, file.path(pdf_dir, paste0(filename, ".pdf")), overwrite = TRUE)
  invisible(c(png_path, pdf_path))
}

short_identity <- function(x) gsub("__", "-", x, fixed = TRUE)

matrix_long <- function(x) {
  output <- as.data.frame(as.table(x), stringsAsFactors = FALSE)
  names(output) <- c("gene", "reference_identity", "z")
  output
}

expected_order <- c("NF1", "NF2", "NF3", "cLTMR", "NP", "PEP", "SST", "TRPM8")
expected_colors <- setNames(anatomic_discrete_palette(length(expected_order)), expected_order)
line_order <- c("TGN1003A", "TGN1028B", "TGN1215AvX")
line_colors <- setNames(c("#0072B2", "#009E73", "#CC79A7"), line_order)
identity_colors <- c(
  "A-beta-LTMR-like" = "#0072B2",
  "C-LTMR-like" = "#D55E00",
  "A-delta-LTMR-like" = "#009E73",
  "Other neuronal" = "#9A9492",
  "Non-neuronal / off-target" = "#D8D4D2"
)

sample_metadata <- read_csv("results/incoming_qc/incoming_sample_metadata.csv", show_col_types = FALSE)
complete_lines <- sample_metadata |>
  filter(sample_class == "TG bulk", !is.na(maturation_week)) |>
  distinct(specimen_id, maturation_week) |>
  count(specimen_id, name = "n_weeks") |>
  filter(n_weeks == 4L) |>
  pull(specimen_id)
stopifnot(setequal(complete_lines, line_order))

combined_metadata <- read_csv(
  "results/pseudobulk/GSE197289_TG_combined_reference_metadata.csv",
  show_col_types = FALSE
)
combined_logcpm <- as.matrix(read.csv(
  "results/pseudobulk/GSE197289_TG_combined_reference_logCPM.csv",
  row.names = 1, check.names = FALSE
))
storage.mode(combined_logcpm) <- "numeric"

neuron_object <- readRDS("results/consensus_identity/GSE197289_human_TG_neurons_consensus.rds")
umap <- Embeddings(neuron_object, "umap_sct") |>
  as.data.frame() |>
  setNames(c("UMAP_1", "UMAP_2")) |>
  rownames_to_column("cell") |>
  left_join(neuron_object[[]] |> rownames_to_column("cell"), by = "cell") |>
  mutate(
    expected_type = factor(expected_type, levels = expected_order),
    final_label = if_else(
      consensus_status == "Donor-supported anchor",
      short_identity(consensus_identity),
      "Unresolved"
    )
  )

final_centroids <- umap |>
  filter(consensus_status == "Donor-supported anchor") |>
  group_by(consensus_identity, expected_type) |>
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop") |>
  mutate(label = sub(".*__", "", consensus_identity))

p_umap_author <- ggplot(umap, aes(UMAP_1, UMAP_2, color = expected_type)) +
  geom_point(size = 0.22, alpha = 0.58, stroke = 0) +
  coord_equal() +
  scale_color_manual(values = expected_colors, drop = FALSE) +
  labs(
    tag = "A", title = "Author types",
    subtitle = "Author-provided neuronal annotation",
    x = NULL, y = NULL, color = "Author type"
  ) +
  panel_theme +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
    panel.grid = element_blank(), legend.position = "bottom",
    plot.tag.position = c(-0.045, 1.01),
    plot.title = element_text(margin = margin(l = 8, b = 2.5)),
    legend.box.margin = margin(t = -2), legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.2, alpha = 1)))

p_umap_final <- ggplot(umap, aes(UMAP_1, UMAP_2, color = expected_type)) +
  geom_point(size = 0.22, alpha = 0.5, stroke = 0) +
  geom_point(
    data = filter(umap, consensus_status != "Donor-supported anchor"),
    shape = 4, size = 0.75, stroke = 0.35, color = ink
  ) +
  geom_label_repel(
    data = final_centroids,
    aes(label = label), size = 1.65, label.size = 0.12,
    label.padding = unit(0.08, "lines"), box.padding = 0.22,
    min.segment.length = 0, seed = 20260804, show.legend = FALSE,
    color = ink, fill = alpha("white", 0.84), max.overlaps = Inf
  ) +
  coord_equal() +
  scale_color_manual(values = expected_colors, drop = FALSE) +
  labs(
    title = "Donor-supported labels",
    subtitle = "Pseudobulk-defined SC labels; crosses are unresolved",
    x = NULL, y = NULL
  ) +
  panel_theme +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
    panel.grid = element_blank(), legend.position = "none"
  )

p_reference_labels <- p_umap_author | p_umap_final

marker_modules <- tribble(
  ~gene, ~module,
  "RBFOX3", "Pan-neuronal", "SNAP25", "Pan-neuronal", "SCN10A", "Pan-neuronal",
  "NTRK1", "Peptidergic", "TRPV1", "Peptidergic", "CALCA", "Peptidergic",
  "P2RX3", "Non-peptidergic", "RET", "Non-peptidergic", "MRGPRX1", "Non-peptidergic",
  "NEFH", "Mechanosensory", "NTRK3", "Mechanosensory", "PVALB", "Mechanosensory",
  "NTRK2", "Mechanosensory", "PIEZO2", "Mechanosensory",
  "CASQ2", "C-LTMR", "P2RY1", "C-LTMR",
  "TRPM8", "Cold / pruriceptor", "SST", "Cold / pruriceptor", "NMB", "Cold / pruriceptor",
  "SOX10", "Counter-programs", "S100B", "Counter-programs",
  "PTPRC", "Counter-programs", "VWF", "Counter-programs"
) |>
  filter(gene %in% rownames(combined_logcpm)) |>
  mutate(
    module = factor(module, levels = unique(module)),
    gene = factor(gene, levels = rev(unique(gene)))
  )

identity_order <- combined_metadata |>
  distinct(reference_level, reference_label, expected_type) |>
  mutate(
    expected_type = factor(expected_type, levels = expected_order),
    level_order = if_else(reference_level == "Neuronal subtype", 1L, 2L)
  ) |>
  arrange(level_order, expected_type, reference_label) |>
  pull(reference_label)

identity_expression <- sapply(identity_order, function(identity) {
  ids <- combined_metadata$pseudobulk_id[combined_metadata$reference_label == identity]
  rowMeans(combined_logcpm[as.character(marker_modules$gene), ids, drop = FALSE])
})
rownames(identity_expression) <- as.character(marker_modules$gene)
identity_z <- t(scale(t(identity_expression)))
identity_z[!is.finite(identity_z)] <- 0
identity_z[identity_z < -2.5] <- -2.5
identity_z[identity_z > 2.5] <- 2.5

marker_heatmap <- matrix_long(identity_z) |>
  left_join(marker_modules |> mutate(gene = as.character(gene)), by = "gene") |>
  mutate(
    reference_identity = factor(reference_identity, levels = identity_order),
    short_reference = factor(short_identity(reference_identity), levels = short_identity(identity_order)),
    gene = factor(gene, levels = levels(marker_modules$gene))
  )

p_marker_heatmap <- ggplot(marker_heatmap, aes(short_reference, gene, fill = z)) +
  geom_tile() +
  facet_grid(module ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, limits = c(-2.5, 2.5)) +
  labs(
    tag = "B", title = "Reference neuronal signature and counter-programs",
    subtitle = "Donor-pseudobulk identity means; gene-wise z-scores",
    x = NULL, y = NULL, fill = "z-score"
  ) +
  panel_theme +
  theme(
    axis.text.x = element_text(angle = 58, hjust = 1, vjust = 1, size = 5.4),
    axis.text.y = element_text(size = 5.5, face = "italic"),
    axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
    strip.placement = "outside", strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, size = 5.6, hjust = 1),
    legend.position = "bottom", legend.key.width = unit(0.65, "cm")
  ) +
  guides(fill = guide_colorbar(title.position = "top", barheight = unit(0.18, "cm")))

signed_scores <- read_csv("results/marker_refinement/refined_reference_signed_scores.csv", show_col_types = FALSE)
signed_calls <- read_csv("results/marker_refinement/refined_reference_signed_score_calls.csv", show_col_types = FALSE)
signed_order <- combined_metadata |>
  filter(reference_level == "Neuronal subtype") |>
  distinct(reference_label, expected_type) |>
  mutate(expected_type = factor(expected_type, levels = expected_order)) |>
  arrange(expected_type, reference_label) |>
  pull(reference_label)

signed_matrix <- signed_scores |>
  group_by(true_label, target_label) |>
  summarise(mean_signed_score = mean(signed_score), .groups = "drop") |>
  mutate(
    true_label = factor(true_label, levels = rev(signed_order)),
    target_label = factor(target_label, levels = signed_order),
    diagonal = as.character(true_label) == as.character(target_label)
  )

vif <- read_csv("results/marker_refinement/vif_before_after_expansion.csv", show_col_types = FALSE)
mean_vif_before <- mean(vif$vif_before)
mean_vif_after <- mean(vif$vif_after)
signed_accuracy <- mean(signed_calls$correct_label)

p_signed_matrix <- ggplot(signed_matrix, aes(target_label, true_label, fill = mean_signed_score)) +
  geom_tile() +
  geom_tile(data = filter(signed_matrix, diagonal), fill = NA, color = ink, linewidth = 0.28) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) +
  scale_x_discrete(labels = short_identity) +
  scale_y_discrete(labels = short_identity) +
  labs(
    tag = "C", title = "Donor-aware signatures separate refined identities",
    subtitle = "60 positive markers/identity enter NNLS; 10 anti-markers/identity inform signed scores only",
    caption = sprintf(
      "Original author-marker matrix mean VIF %.1f -> donor-validated positive-marker matrix %.2f; signed reference calls %d/%d (%.1f%%).",
      mean_vif_before, mean_vif_after, sum(signed_calls$correct_label), nrow(signed_calls), 100 * signed_accuracy
    ),
    x = "Scored identity", y = "True donor pseudobulk identity", fill = "Positive - anti"
  ) +
  panel_theme +
  theme(
    axis.text.x = element_text(angle = 58, hjust = 1, size = 5.5),
    axis.text.y = element_text(size = 5.5), axis.ticks = element_blank(),
    axis.line = element_blank(), panel.grid = element_blank(), legend.position = "bottom",
    legend.key.width = unit(0.75, "cm")
  ) +
  guides(fill = guide_colorbar(title.position = "top", barheight = unit(0.18, "cm")))

projection <- read_csv("results/figures/RealTGN_bulk_neuron_only_pca_coordinates.csv", show_col_types = FALSE)
projection_reference <- projection |>
  filter(sample_origin == "Reference") |>
  mutate(expected_type = factor(expected_type, levels = expected_order), donor = factor(donor))
projection_bulk <- projection |>
  filter(sample_origin == "Bulk", sample_class == "TG bulk") |>
  left_join(sample_metadata |> select(sample_id, maturation_week), by = "sample_id") |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  mutate(specimen_id = factor(specimen_id, levels = line_order))
projection_centroids <- projection_reference |>
  group_by(expected_type) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

p_projection_context <- ggplot() +
  geom_point(
    data = projection_reference,
    aes(PC1, PC2, color = expected_type, shape = donor), size = 1.4, alpha = 0.52
  ) +
  geom_text_repel(
    data = projection_centroids, aes(PC1, PC2, label = expected_type, color = expected_type),
    size = 1.8, seed = 20260804, min.segment.length = 0, box.padding = 0.15,
    show.legend = FALSE, max.overlaps = Inf
  ) +
  geom_path(
    data = projection_bulk, aes(PC1, PC2, group = specimen_id, color = specimen_id),
    linewidth = 0.8, arrow = arrow(length = unit(0.07, "inches"))
  ) +
  geom_point(
    data = projection_bulk, aes(PC1, PC2, fill = maturation_week),
    shape = 21, size = 2.5, color = ink, stroke = 0.35
  ) +
  scale_color_manual(values = c(expected_colors, line_colors), breaks = line_order) +
  scale_fill_viridis_c(option = "C", begin = 0.08, end = 0.92, breaks = 1:4) +
  labs(
    tag = "D", title = "Cultures traverse a fixed neuron-only reference",
    subtitle = "Donor pseudobulks define axes; arrows connect weeks 1-4 by biological line",
    x = "Reference PC1", y = "Reference PC2", color = "Biological line", fill = "Week", shape = "Donor"
  ) +
  panel_theme +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 1.1, size = 2)),
    fill = guide_colorbar(order = 2, title.position = "top", barheight = unit(0.18, "cm")),
    shape = "none"
  )

trajectory_pca <- read_csv(
  "results/incoming_qc/weekly_neuron_only_pca_trajectory_coordinates.csv",
  show_col_types = FALSE
) |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  mutate(specimen_id = factor(specimen_id, levels = line_order))

trajectory_centroids <- trajectory_pca |>
  group_by(maturation_week) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

p_projection <- ggplot(
  trajectory_pca,
  aes(PC1, PC2, color = specimen_id, group = specimen_id)
) +
  geom_path(
    linewidth = 0.95,
    arrow = arrow(length = unit(0.09, "inches"), type = "closed")
  ) +
  geom_point(aes(fill = maturation_week), shape = 21, size = 3.2, color = ink, stroke = 0.45) +
  geom_point(
    data = trajectory_centroids, aes(PC1, PC2), inherit.aes = FALSE,
    shape = 23, size = 4, fill = red, color = ink, stroke = 0.45
  ) +
  geom_text_repel(
    data = trajectory_centroids,
    aes(PC1, PC2, label = paste0("W", maturation_week)),
    inherit.aes = FALSE, size = 2.1, color = ink,
    box.padding = 0.25, point.padding = 0.55,
    min.segment.length = 0, seed = 20260805
  ) +
  scale_color_manual(values = line_colors, drop = FALSE) +
  scale_fill_viridis_c(option = "C", begin = 0.08, end = 0.92, breaks = 1:4) +
  labs(
    tag = "D", title = "Centroid movement in fixed neuron-only PCA space",
    subtitle = "Paths show biological lines; diamonds show week-specific centroids",
    x = "Reference PC1", y = "Reference PC2", color = "Biological line", fill = "Week"
  ) +
  panel_theme +
  theme(legend.position = "bottom")

dominant_lookup <- c(
  "NF2__SC03" = "A-beta-LTMR-like",
  "cLTMR__SC01" = "C-LTMR-like",
  "NF3__SC02" = "A-delta-LTMR-like"
)
dominant_order <- unname(dominant_lookup)
culture_signed <- read_csv(
  "results/deconv_rerun_expanded/weekly_signed_refined_identity_scores.csv",
  show_col_types = FALSE
) |>
  filter(
    specimen_id %in% complete_lines, maturation_week %in% 1:4,
    refined_identity %in% names(dominant_lookup)
  ) |>
  mutate(
    identity = recode(refined_identity, !!!dominant_lookup),
    identity = factor(identity, levels = dominant_order),
    specimen_id = factor(specimen_id, levels = line_order)
  )
culture_signed_mean <- culture_signed |>
  group_by(identity, maturation_week) |>
  summarise(mean_score = mean(signed_score), .groups = "drop")

p_subtype_scores <- ggplot(culture_signed, aes(maturation_week, signed_score, group = specimen_id)) +
  geom_hline(yintercept = 0, color = "#CFCBCA", linewidth = 0.35) +
  geom_line(aes(color = specimen_id), linewidth = 0.45, alpha = 0.55) +
  geom_point(aes(color = specimen_id), size = 1.2, alpha = 0.75) +
  geom_line(
    data = culture_signed_mean,
    aes(maturation_week, mean_score, color = identity, group = identity),
    linewidth = 1.05, inherit.aes = FALSE
  ) +
  geom_point(
    data = culture_signed_mean,
    aes(maturation_week, mean_score, fill = identity),
    shape = 21, size = 2.1, color = "white", stroke = 0.35, inherit.aes = FALSE
  ) +
  facet_wrap(~identity, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = 1:4) +
  scale_color_manual(values = c(line_colors, identity_colors[dominant_order]), breaks = line_order) +
  scale_fill_manual(values = identity_colors[dominant_order]) +
  labs(
    tag = "E", title = "Subtype-resolved positive-minus-anti-marker evidence",
    subtitle = "Thin paths are matched biological lines; bold colored paths are line means",
    x = "Culture week", y = "Signed score", color = "Biological line"
  ) +
  panel_theme +
  theme(legend.position = "bottom") +
  guides(fill = "none", color = guide_legend(nrow = 1, override.aes = list(linewidth = 0.9)))

composition <- read_csv(
  "results/deconv_rerun_expanded/weekly_nnls_deconvolution_with_bootstrap_ci.csv",
  show_col_types = FALSE
) |>
  filter(
    specimen_id %in% complete_lines, maturation_week %in% 1:4,
    manifest == "Signed refined (788 genes)"
  ) |>
  left_join(
    combined_metadata |> distinct(reference_identity = reference_label, reference_level),
    by = "reference_identity"
  ) |>
  mutate(
    identity = case_when(
      reference_identity == "NF2__SC03" ~ "A-beta-LTMR-like",
      reference_identity == "cLTMR__SC01" ~ "C-LTMR-like",
      reference_identity == "NF3__SC02" ~ "A-delta-LTMR-like",
      reference_level == "Neuronal subtype" ~ "Other neuronal",
      TRUE ~ "Non-neuronal / off-target"
    ),
    specimen_id = factor(specimen_id, levels = line_order)
  ) |>
  group_by(specimen_id, maturation_week, identity) |>
  summarise(contribution = sum(estimated_fraction), .groups = "drop") |>
  complete(
    specimen_id, maturation_week,
    identity = names(identity_colors), fill = list(contribution = 0)
  ) |>
  mutate(identity = factor(identity, levels = rev(names(identity_colors))))

week_means <- composition |>
  mutate(identity = as.character(identity)) |>
  group_by(maturation_week, identity) |>
  summarise(contribution = mean(contribution), .groups = "drop")
week4_abeta <- week_means$contribution[week_means$maturation_week == 4 & week_means$identity == "A-beta-LTMR-like"]
week4_cltmr <- week_means$contribution[week_means$maturation_week == 4 & week_means$identity == "C-LTMR-like"]

p_composition <- ggplot(composition, aes(factor(maturation_week), contribution, fill = identity)) +
  geom_col(width = 0.76, color = "white", linewidth = 0.12) +
  facet_wrap(~specimen_id, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = identity_colors, breaks = names(identity_colors)) +
  labs(
    tag = "F", title = "A-beta-LTMR-like signal dominates the refined NNLS profile",
    subtitle = sprintf(
      "Week 4 mean: %.1f%% A-beta-LTMR-like and %.1f%% C-LTMR-like",
      100 * week4_abeta, 100 * week4_cltmr
    ),
    caption = "Relative transcriptional signature contributions; not literal cell fractions.",
    x = "Culture week", y = "Relative contribution", fill = NULL
  ) +
  panel_theme +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

primary_row_1 <- wrap_plots(wrap_elements(full = p_reference_labels), p_marker_heatmap,
                            ncol = 2, widths = c(1.04, 0.96))
primary_row_2 <- wrap_plots(p_signed_matrix, p_projection, ncol = 2)
primary_row_3 <- wrap_plots(p_subtype_scores, p_composition, ncol = 2)

primary_figure <- wrap_plots(
  wrap_elements(full = primary_row_1), wrap_elements(full = primary_row_2),
  wrap_elements(full = primary_row_3), ncol = 1, heights = c(1.02, 1.03, 0.92)
) +
  plot_annotation(
    title = "A donor-aware human TG reference anchors sensory-neuronal identity in developing cultures",
    subtitle = paste0("Reference definition and signature construction precede projection and subtype interpretation.\n",
                      "NNLS uses positive markers; anti-markers are reserved for signed scoring."),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12.5, color = ink, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.5, color = "#55525D", margin = margin(b = 5)),
      plot.margin = margin(6, 7, 6, 7)
    )
  )

save_figure("Figure_1_primary_TG_neuronal_identity_story", primary_figure, 8.5, 11)

# -----------------------------------------------------------------------------
# Extended Data Figure 1: reference construction and complete marker programs
# -----------------------------------------------------------------------------

support <- read_csv(
  "results/pseudobulk/GSE197289_neuronal_refined_pseudobulk_manifest.csv",
  show_col_types = FALSE
) |>
  filter(refined_neuronal_label %in% signed_order) |>
  mutate(
    refined_neuronal_label = factor(refined_neuronal_label, levels = signed_order),
    donor = factor(donor),
    support = case_when(
      include_pseudobulk ~ "Eligible pseudobulk",
      donor_meets_cell_minimum ~ "Excluded for label support",
      TRUE ~ "<20 nuclei"
    )
  )

p_support <- ggplot(support, aes(donor, refined_neuronal_label, fill = n_nuclei)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_point(data = filter(support, include_pseudobulk), shape = 21, size = 2.2, fill = "white", color = ink) +
  geom_text(aes(label = n_nuclei), size = 2.2, color = ink) +
  scale_color_identity() +
  scale_fill_gradient(low = "#F1F5F8", high = blue, trans = "sqrt") +
  scale_y_discrete(labels = short_identity) +
  labs(
    tag = "A", title = "Donor support defines eligible neuronal pseudobulks",
    subtitle = "Numbers are nuclei per donor-by-refined-label intersection; circles mark included pseudobulks",
    x = "Donor", y = "Refined neuronal label", fill = "Nuclei"
  ) +
  ed_theme + theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")

full_manifest <- read_csv(
  "results/pseudobulk/GSE197289_TG_marker_heatmap_gene_manifest.csv",
  show_col_types = FALSE
) |>
  filter(present, gene %in% rownames(combined_logcpm))
full_gene_order <- rev(unique(full_manifest$gene))
full_expression <- sapply(identity_order, function(identity) {
  ids <- combined_metadata$pseudobulk_id[combined_metadata$reference_label == identity]
  rowMeans(combined_logcpm[full_manifest$gene, ids, drop = FALSE])
})
rownames(full_expression) <- full_manifest$gene
full_z <- t(scale(t(full_expression)))
full_z[!is.finite(full_z)] <- 0
full_z[full_z < -2.5] <- -2.5
full_z[full_z > 2.5] <- 2.5
full_heatmap <- matrix_long(full_z) |>
  left_join(full_manifest |> select(gene, marker_group), by = "gene") |>
  mutate(
    gene = factor(gene, levels = full_gene_order),
    reference_identity = factor(reference_identity, levels = identity_order)
  )

p_full_heatmap <- ggplot(full_heatmap, aes(reference_identity, gene, fill = z)) +
  geom_tile() +
  facet_grid(marker_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, limits = c(-2.5, 2.5)) +
  scale_x_discrete(labels = short_identity) +
  labs(
    tag = "B", title = "Complete prespecified marker-program heatmap",
    subtitle = "Identity means across donor pseudobulks; gene-wise z-scores",
    x = NULL, y = NULL, fill = "z-score"
  ) +
  ed_theme +
  theme(
    axis.text.x = element_text(angle = 58, hjust = 1, size = 6.1),
    axis.text.y = element_text(size = 5.4, face = "italic"),
    axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
    strip.placement = "outside", strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, size = 6.2), legend.position = "bottom"
  )

neuron_pca <- read_csv(
  "results/pseudobulk/GSE197289_neuron_only_refined_pca_coordinates.csv",
  show_col_types = FALSE
) |>
  mutate(expected_type = factor(expected_type, levels = expected_order), donor = factor(donor))

p_neuron_pca <- ggplot(neuron_pca, aes(PC1, PC2, color = expected_type, shape = donor)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text_repel(aes(label = short_identity(refined_neuronal_label)), size = 2.2, seed = 20260804,
                  min.segment.length = 0, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = expected_colors, drop = FALSE) +
  labs(
    tag = "C", title = "Eligible donor pseudobulks retain refined neuronal structure",
    subtitle = "Neuron-only PCA of donor-by-refined-label profiles",
    x = "Reference PC1", y = "Reference PC2", color = "Author type", shape = "Donor"
  ) +
  ed_theme + theme(legend.position = "bottom")

p_unresolved <- ggplot(umap, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "#CFCBCA", size = 0.35, alpha = 0.45) +
  geom_point(
    data = filter(umap, consensus_status != "Donor-supported anchor"),
    aes(color = expected_type), shape = 4, size = 1.2, stroke = 0.5
  ) +
  scale_color_manual(values = expected_colors, drop = FALSE) + coord_equal() +
  labs(
    tag = "D", title = "Low-support nuclei: unresolved",
    subtitle = sprintf("%d nuclei without a supported subtype call", sum(umap$consensus_status != "Donor-supported anchor")),
    x = "UMAP 1", y = "UMAP 2", color = "Author type"
  ) +
  ed_theme + theme(legend.position = "none", panel.grid = element_blank())

extended_1 <- (p_support | p_full_heatmap) / (p_neuron_pca | p_unresolved) +
  plot_layout(widths = c(0.8, 1.2), heights = c(1.15, 0.85)) +
  plot_annotation(
    title = "Extended Data Figure 1 | Donor support and complete human TG reference programs",
    theme = theme(plot.title = element_text(face = "bold", size = 14, color = ink), plot.margin = margin(6, 8, 6, 8))
  )
save_figure("Extended_Data_Figure_1_reference_construction", extended_1, 14, 10)

# -----------------------------------------------------------------------------
# Extended Data Figure 2: marker/signature validation
# -----------------------------------------------------------------------------

expanded_manifest <- read_csv(
  "results/marker_refinement/GSE197289_TG_expanded_marker_manifest.csv",
  show_col_types = FALSE
) |>
  filter(marker_group %in% signed_order, direction %in% c("positive", "negative"))
neuronal_metadata <- combined_metadata |> filter(reference_level == "Neuronal subtype")
neuronal_ids <- neuronal_metadata$pseudobulk_id
signed_genes <- intersect(unique(expanded_manifest$gene), rownames(combined_logcpm))
signed_expression <- combined_logcpm[signed_genes, neuronal_ids, drop = FALSE]
gene_center <- rowMeans(signed_expression)
gene_scale <- apply(signed_expression, 1, sd)
gene_scale[!is.finite(gene_scale) | gene_scale == 0] <- 1
signed_z <- sweep(signed_expression, 1, gene_center, "-")
signed_z <- sweep(signed_z, 1, gene_scale, "/")

component_scores <- bind_rows(lapply(signed_order, function(target) {
  panel <- expanded_manifest |> filter(marker_group == target, gene %in% rownames(signed_z))
  bind_rows(lapply(c("positive", "negative"), function(direction_value) {
    genes <- panel |> filter(direction == direction_value)
    score <- colSums(signed_z[genes$gene, , drop = FALSE] * genes$marker_weight) / sum(genes$marker_weight)
    tibble(pseudobulk_id = names(score), target_label = target, component = direction_value, score = score)
  }))
})) |>
  left_join(neuronal_metadata |> transmute(pseudobulk_id, true_label = reference_label), by = "pseudobulk_id") |>
  group_by(component, true_label, target_label) |>
  summarise(score = mean(score), .groups = "drop") |>
  mutate(
    true_label = factor(true_label, levels = rev(signed_order)),
    target_label = factor(target_label, levels = signed_order)
  )

component_limits <- max(abs(component_scores$score), na.rm = TRUE)
p_components <- ggplot(component_scores, aes(target_label, true_label, fill = score)) +
  geom_tile() +
  facet_wrap(~component, nrow = 1, labeller = as_labeller(c(
    positive = "Positive-marker evidence", negative = "Anti-marker expression"
  ))) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0,
                       limits = c(-component_limits, component_limits)) +
  scale_x_discrete(labels = short_identity) + scale_y_discrete(labels = short_identity) +
  labs(
    tag = "A", title = "Positive and anti-marker components are retained as separate evidence",
    subtitle = "Anti-markers are subtracted in signed scoring and are not appended to the NNLS design matrix",
    x = "Scored identity", y = "True donor pseudobulk identity", fill = "Mean score"
  ) +
  ed_theme +
  theme(
    axis.text.x = element_text(angle = 58, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6), axis.ticks = element_blank(),
    axis.line = element_blank(), panel.grid = element_blank(), legend.position = "bottom"
  )

vif_order <- vif |> arrange(vif_before) |> pull(reference_identity)
vif_long <- vif |>
  pivot_longer(c(vif_before, vif_after), names_to = "matrix", values_to = "vif") |>
  mutate(
    matrix = recode(matrix,
      vif_before = "Original author-marker matrix",
      vif_after = "Donor-validated positive-marker matrix"
    ),
    reference_identity = factor(reference_identity, levels = vif_order)
  )
p_vif <- ggplot(vif_long, aes(reference_identity, vif, fill = matrix)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_hline(yintercept = 10, color = red, linetype = "dashed", linewidth = 0.55) +
  coord_flip() +
  scale_fill_manual(values = c("#C8C3C1", "#0072B2")) +
  labs(
    tag = "B", title = "Donor-validated positive-marker expansion improves matrix conditioning",
    subtitle = sprintf("Mean VIF %.1f -> %.2f; dashed line marks VIF 10", mean_vif_before, mean_vif_after),
    x = NULL, y = "Variance inflation factor", fill = NULL
  ) +
  ed_theme + theme(legend.position = "bottom")

call_margin <- signed_scores |>
  group_by(pseudobulk_id, true_label, donor) |>
  arrange(desc(signed_score), .by_group = TRUE) |>
  summarise(
    top_call = first(target_label),
    top_score = first(signed_score),
    runner_up = nth(target_label, 2),
    runner_score = nth(signed_score, 2),
    margin = top_score - runner_score,
    .groups = "drop"
  ) |>
  mutate(correct = top_call == true_label)
p_call_margin <- ggplot(call_margin, aes(reorder(short_identity(true_label), margin), margin, color = donor)) +
  geom_hline(yintercept = 0, color = red, linetype = "dashed") +
  geom_point(aes(shape = correct), size = 2.7) +
  coord_flip() +
  labs(
    tag = "C", title = "Signed scoring correctly calls 38 of 39 donor pseudobulks",
    subtitle = "Margin is the top signed score minus the runner-up score",
    x = "True refined identity", y = "Top-call margin", color = "Donor", shape = "Correct"
  ) +
  ed_theme + theme(legend.position = "bottom")

recovery <- read_csv(
  "results/deconvolution_diagnostics/donor_held_out_mixture_recovery.csv",
  show_col_types = FALSE
) |>
  group_by(held_out_donor, scenario) |>
  summarise(
    target_recovery = sum(recovered_fraction[true_fraction > 0]),
    l1_error = sum(abs(recovered_fraction - true_fraction)) / 2,
    scenario_type = if (startsWith(first(scenario), "pure:")) "Pure identity" else "50:50 nearest-neighbor mixture",
    .groups = "drop"
  )
p_recovery <- ggplot(recovery, aes(target_recovery, 1 - l1_error, color = factor(held_out_donor), shape = scenario_type)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#AAA5A3") +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "#AAA5A3") +
  geom_point(size = 2.3, alpha = 0.8) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    tag = "D", title = "Held-out recovery remains incomplete",
    subtitle = "Pure identities and 50:50 mixtures from unseen donors",
    x = "Recovered mass on true component(s)", y = "One - L1 error", color = "Held-out donor", shape = "Scenario"
  ) +
  ed_theme + theme(legend.position = "bottom")

extended_2 <- (p_components | p_vif) / (p_call_margin | p_recovery) +
  plot_layout(widths = c(1.18, 0.82), heights = c(1.05, 0.95)) +
  plot_annotation(
    title = "Extended Data Figure 2 | Donor-aware signature construction and validation",
    theme = theme(plot.title = element_text(face = "bold", size = 14, color = ink), plot.margin = margin(6, 8, 6, 8))
  )
save_figure("Extended_Data_Figure_2_signature_validation", extended_2, 14, 10)

# -----------------------------------------------------------------------------
# Extended Data Figure 3: incoming context and projection sensitivity
# -----------------------------------------------------------------------------

qc_flags <- read_csv("results/incoming_qc/incoming_sample_qc_flags.csv", show_col_types = FALSE)
p_qc <- ggplot(qc_flags, aes(total_counts, detected_genes, color = sample_class, shape = qc_status)) +
  geom_vline(xintercept = 5e6, linetype = "dashed", color = red) +
  geom_hline(yintercept = 1e4, linetype = "dashed", color = red) +
  geom_point(size = 2.8, alpha = 0.82) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    tag = "A", title = "Incoming libraries pass operational QC floors",
    subtitle = "Dashed lines: 5 million counts and 10,000 detected genes",
    x = "Total counts", y = "Detected genes", color = "Sample class", shape = "QC status"
  ) + ed_theme + theme(legend.position = "bottom")

program_scores <- read_csv(
  "results/incoming_qc/weekly_maturation_program_scores.csv",
  show_col_types = FALSE
) |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  pivot_longer(-c(specimen_id, sample_class, maturation_week), names_to = "program", values_to = "score") |>
  mutate(specimen_id = factor(specimen_id, levels = line_order))
p_programs <- ggplot(program_scores, aes(maturation_week, score, group = specimen_id, color = specimen_id)) +
  geom_hline(yintercept = 0, color = "#D4D0CE") +
  geom_line(linewidth = 0.48, alpha = 0.75) + geom_point(size = 1.25) +
  facet_wrap(~program, scales = "free_y", ncol = 5) +
  scale_x_continuous(breaks = 1:4) + scale_color_manual(values = line_colors) +
  labs(
    tag = "B", title = "Complete reference-weighted program trajectories",
    subtitle = "Broad neuronal and counter-program views are retained as supporting context",
    x = "Culture week", y = "Program score", color = "Biological line"
  ) + ed_theme + theme(legend.position = "bottom")

correlation <- read_csv("results/figures/RealTGN_bulk_top_reference_matches.csv", show_col_types = FALSE) |>
  inner_join(
    sample_metadata |> filter(specimen_id %in% complete_lines) |> select(sample_id, specimen_id, maturation_week),
    by = c("bulk_sample" = "sample_id")
  ) |>
  group_by(bulk_sample) |>
  slice_max(median_correlation, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(specimen_id = factor(specimen_id, levels = line_order))
p_correlation <- ggplot(correlation, aes(maturation_week, reference_label, color = median_correlation)) +
  geom_point(size = 3.3) + facet_wrap(~specimen_id, nrow = 1) +
  scale_x_continuous(breaks = 1:4) + scale_color_viridis_c(option = "C", end = 0.9) +
  labs(
    tag = "C", title = "Top reference match is tracked independently of PCA",
    subtitle = "Color is median donor-level bulk-to-reference Spearman correlation",
    x = "Culture week", y = "Top reference identity", color = "Correlation"
  ) + ed_theme + theme(legend.position = "bottom")

sensitivity_metrics <- read_csv("results/consensus_identity/pc_sensitivity_metrics.csv", show_col_types = FALSE)
p_sensitivity <- sensitivity_metrics |>
  select(n_pcs, cumulative_identity_centroid_variance, umap_reference_knn_overlap, tsne_reference_knn_overlap) |>
  pivot_longer(-c(n_pcs, cumulative_identity_centroid_variance), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
    umap_reference_knn_overlap = "Reference-trained UMAP",
    tsne_reference_knn_overlap = "Exploratory joint t-SNE"
  )) |>
  ggplot(aes(n_pcs, value, color = metric)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  scale_x_continuous(breaks = sensitivity_metrics$n_pcs) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(
    tag = "D", title = "Nonlinear views remain sensitivity analyses",
    subtitle = "Mean 5-nearest-neighbor overlap across retained neuronal PC counts",
    x = "Retained neuronal PCs", y = "Neighborhood preservation", color = NULL
  ) + ed_theme + theme(legend.position = "bottom")

extended_3 <- (p_qc | p_correlation) / p_programs / p_sensitivity +
  plot_layout(heights = c(0.9, 1.3, 0.72), widths = c(1, 1)) +
  plot_annotation(
    title = "Extended Data Figure 3 | Incoming-sample context and reference-mapping sensitivity",
    theme = theme(plot.title = element_text(face = "bold", size = 14, color = ink), plot.margin = margin(6, 8, 6, 8))
  )
save_figure("Extended_Data_Figure_3_mapping_context", extended_3, 12, 12)

# -----------------------------------------------------------------------------
# Extended Data Figure 4: composition robustness and interpretation boundary
# -----------------------------------------------------------------------------

composition_all <- read_csv(
  "results/deconv_rerun_expanded/weekly_nnls_deconvolution_with_bootstrap_ci.csv",
  show_col_types = FALSE
) |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  left_join(
    combined_metadata |> distinct(reference_identity = reference_label, reference_level),
    by = "reference_identity"
  ) |>
  mutate(
    identity = case_when(
      reference_identity == "NF2__SC03" ~ "A-beta-LTMR-like",
      reference_identity == "cLTMR__SC01" ~ "C-LTMR-like",
      reference_identity == "NF3__SC02" ~ "A-delta-LTMR-like",
      reference_level == "Neuronal subtype" ~ "Other neuronal",
      TRUE ~ "Non-neuronal / off-target"
    ),
    manifest = recode(manifest,
      "Original (79 genes)" = "Original author markers",
      "Signed refined (788 genes)" = "Donor-validated positive markers"
    ),
    manifest = factor(manifest, levels = c("Original author markers", "Donor-validated positive markers")),
    specimen_id = factor(specimen_id, levels = line_order)
  )

before_after <- composition_all |>
  group_by(specimen_id, maturation_week, manifest, identity) |>
  summarise(contribution = sum(estimated_fraction), .groups = "drop") |>
  mutate(identity = factor(identity, levels = rev(names(identity_colors))))
p_before_after <- ggplot(before_after, aes(factor(maturation_week), contribution, fill = identity)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.12) +
  facet_grid(specimen_id ~ manifest) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = identity_colors) +
  labs(
    tag = "A", title = "Marker refinement changes the inferred contribution profile",
    subtitle = "This preserves the informative former Figure 4C comparison in Extended Data",
    x = "Culture week", y = "Relative signature contribution", fill = NULL
  ) + ed_theme + theme(legend.position = "bottom")

bootstrap_summary <- composition_all |>
  filter(manifest == "Donor-validated positive markers", identity %in% c("A-beta-LTMR-like", "C-LTMR-like")) |>
  group_by(specimen_id, maturation_week, identity) |>
  summarise(
    estimate = sum(estimated_fraction),
    ci_low = sum(ci_low), ci_high = sum(ci_high),
    .groups = "drop"
  ) |>
  mutate(identity = factor(identity, levels = c("A-beta-LTMR-like", "C-LTMR-like")))
p_bootstrap <- ggplot(bootstrap_summary, aes(maturation_week, estimate, group = specimen_id, color = specimen_id)) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), linewidth = 0.45, alpha = 0.55) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  facet_wrap(~identity, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = 1:4) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = line_colors) +
  labs(
    tag = "B", title = "Gene-resampling intervals quantify marker sensitivity",
    subtitle = "Intervals are not between-line confidence intervals",
    x = "Culture week", y = "Relative contribution", color = "Biological line"
  ) + ed_theme + theme(legend.position = "bottom")

all_signed_heatmap <- read_csv(
  "results/deconv_rerun_expanded/weekly_signed_refined_identity_scores.csv",
  show_col_types = FALSE
) |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  mutate(
    sample = paste0(sub("TGN", "", specimen_id), " W", maturation_week),
    sample = factor(sample, levels = unlist(lapply(line_order, function(id) paste0(sub("TGN", "", id), " W", 1:4)))),
    refined_identity = factor(refined_identity, levels = rev(signed_order))
  )
p_all_signed <- ggplot(all_signed_heatmap, aes(sample, refined_identity, fill = signed_score)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0) +
  scale_y_discrete(labels = short_identity) +
  labs(
    tag = "C", title = "All refined signed identity scores retain the broader context",
    subtitle = "The main figure displays the three prespecified culture-associated axes",
    x = NULL, y = "Refined identity", fill = "Signed score"
  ) + ed_theme +
  theme(axis.text.x = element_text(angle = 58, hjust = 1), axis.ticks = element_blank(),
        axis.line = element_blank(), panel.grid = element_blank(), legend.position = "bottom")

radar_data <- read_csv("results/consensus_identity/radar_identity_distribution.csv", show_col_types = FALSE)
radar_preview <- radar_data |>
  filter(specimen_id %in% complete_lines, maturation_week %in% 1:4) |>
  mutate(
    identity = factor(radar_identity, levels = unique(radar_identity)),
    specimen_id = factor(specimen_id, levels = line_order)
  )
p_radar_context <- ggplot(radar_preview, aes(maturation_week, estimated_fraction, color = identity, group = identity)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.7) +
  facet_wrap(~specimen_id, nrow = 1) +
  scale_x_continuous(breaks = 1:4) + scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    tag = "D", title = "Five-axis summaries are retained as orientation, not primary evidence",
    subtitle = "Line plots avoid the perceptual distortion of radar-area comparisons",
    x = "Culture week", y = "Relative contribution", color = "Identity axis"
  ) + ed_theme + theme(legend.position = "bottom")

extended_4 <- p_before_after / (p_bootstrap | p_all_signed) / p_radar_context +
  plot_layout(heights = c(1.15, 1, 0.82), widths = c(1, 1)) +
  plot_annotation(
    title = "Extended Data Figure 4 | Composition robustness and interpretation boundary",
    theme = theme(plot.title = element_text(face = "bold", size = 14, color = ink), plot.margin = margin(6, 8, 6, 8))
  )
save_figure("Extended_Data_Figure_4_composition_robustness", extended_4, 12, 13)

caption_lines <- c(
  "# Single-primary figure set",
  "",
  "## Figure 1 | A donor-aware human TG reference resolves the sensory-neuronal identity of developing cultures",
  "",
  "(A) The same neuronal UMAP colored by the author-provided neuronal type and annotated with the final donor-supported refined label. Crosses indicate low-support intersections retained as unresolved. (B) Selected neuronal marker programs and non-neuronal counter-programs in donor-pseudobulk identity means. (C) Mean positive-minus-anti-marker scores across donor pseudobulks. Positive donor-validated markers enter NNLS; anti-markers are retained only for complementary signed scoring. (D) Three matched TG culture lines projected through the fixed neuron-only reference PCA. (E) Longitudinal signed evidence for A-beta-LTMR-like, C-LTMR-like, and A-delta-LTMR-like axes. (F) Refined NNLS relative transcriptional signature contributions, showing a dominant A-beta-LTMR-like profile with a smaller C-LTMR-like contribution. NNLS values are not literal cell fractions.",
  "",
  "## Extended Data Figure 1 | Donor support and complete human TG reference programs",
  "",
  "Donor-by-label pseudobulk eligibility, the complete prespecified marker-program heatmap, neuron-only pseudobulk PCA, and the unresolved low-support intersections.",
  "",
  "## Extended Data Figure 2 | Donor-aware signature construction and validation",
  "",
  "Separate positive-marker and anti-marker components, per-identity VIF before and after donor-validated positive-marker expansion, signed-call margins, and held-out donor mixture recovery.",
  "",
  "## Extended Data Figure 3 | Incoming-sample context and reference-mapping sensitivity",
  "",
  "Incoming library QC, correlation-based top matches, complete broad program trajectories, and nonlinear projection sensitivity.",
  "",
  "## Extended Data Figure 4 | Composition robustness and interpretation boundary",
  "",
  "Original-versus-refined marker composition, gene-bootstrap sensitivity intervals, all signed refined-identity scores, and five-axis contribution context rendered as line plots."
)
writeLines(caption_lines, file.path(run_dir, "FIGURE_CAPTIONS.md"))
caption_lines <- c(
  "# Final figure captions",
  "",
  "## Figure 1 | A donor-aware human TG reference anchors sensory-neuronal identity in developing cultures",
  "",
  "A donor-grounded human trigeminal ganglion (TG) reference was used to define the neuronal state space before interpreting the culture RNA-seq data. The matched longitudinal unit is the biological cell line (n = 3 lines, one RNA-seq library per line at each of weeks 1-4).",
  "",
  "(A) The same single-nucleus UMAP is shown with the author-provided neuronal types (left) and the final donor-supported refined labels (right). Refined SC labels were retained only when eligible donor-by-label pseudobulks contained at least 20 nuclei and the label was supported by at least two donors. Crosses mark low-support intersections retained as unresolved rather than promoted to neuronal subtypes.",
  "(B) Selected neuronal identity programs and non-neuronal counter-programs across donor-pseudobulk identity means. Values are gene-wise z-scores; columns are the 18 donor-supported neuronal identities followed by broad non-neuronal reference classes.",
  "(C) Mean positive-minus-anti-marker scores across donor pseudobulks. For each refined neuronal identity, 60 donor-validated positive markers contribute to the positive-marker NNLS design, whereas 10 anti-markers are reserved for the complementary signed score. The expanded positive-marker matrix reduced mean variance inflation factor (VIF) from 40.18 to 7.73, and signed scoring correctly classified 38 of 39 donor pseudobulks (97.4%). Per-identity conditioning and held-out validation are shown in Extended Data Figure 2.",
  "(D) Matched culture trajectories in the fixed neuron-only PCA trained on donor pseudobulks. Colored paths and arrows connect weeks 1-4 within each biological line; red diamonds mark the week-specific centroid across the three lines. Off-target reference classes were excluded from axis construction.",
  "(E) Positive-minus-anti-marker evidence for the three prespecified culture-associated axes. Thin paths show individual biological lines, and bold paths show the line mean. These signed scores provide identity-specific evidence independently of the NNLS contribution estimates.",
  "(F) Relative transcriptional signature contributions from the donor-validated positive-marker NNLS model. The mean A-beta-LTMR-like contribution was 93.3% at week 1 and 85.7% at week 4; the mean C-LTMR-like contribution was 6.7% at week 1 and 13.8% at week 4. Values are relative reference-signature contributions, not literal histological cell fractions.",
  "",
  "## Extended Data Figure 1 | Donor support and complete human TG reference programs",
  "",
  "(A) Nuclei per donor-by-refined-label intersection. Open circles mark pseudobulks included under the prespecified nucleus and donor-support criteria. (B) Complete prespecified marker-program heatmap across donor-pseudobulk identity means, shown as gene-wise z-scores. (C) Neuron-only PCA of eligible donor pseudobulks, colored by author type and shaped by donor. (D) UMAP locations of 211 low-support nuclei retained without a donor-supported subtype call.",
  "",
  "## Extended Data Figure 2 | Donor-aware signature construction and validation",
  "",
  "(A) Positive-marker evidence and anti-marker expression shown separately across true and scored donor-pseudobulk identities. Anti-markers are subtracted only in signed scoring and are not appended to the NNLS design matrix. (B) Per-identity VIF for the original author-marker matrix and the donor-validated positive-marker matrix; the dashed line marks VIF 10. (C) Margin between the top signed call and runner-up for each donor pseudobulk; 38 of 39 top calls matched the true refined identity. (D) Recovery of pure identities and 50:50 nearest-neighbor mixtures from held-out donors. Residual recovery error defines the boundary for subtype-level interpretation.",
  "",
  "## Extended Data Figure 3 | Incoming-sample context and reference-mapping sensitivity",
  "",
  "(A) Library depth and detected-gene complexity relative to the operational QC floors of 5 million counts and 10,000 detected genes. (B) Complete broad neuronal and counter-program trajectories for the three matched biological lines. These supporting views are not used as the primary subtype result. (C) Highest-matching reference identity for each matched TG culture library, colored by median donor-level Spearman correlation. (D) Preservation of five-nearest-neighbor reference structure for nonlinear embeddings using 7, 10, or 15 neuronal PCs. Fixed-reference PCA remains the quantitative trajectory representation; UMAP and t-SNE are sensitivity views.",
  "",
  "## Extended Data Figure 4 | Composition robustness and interpretation boundary",
  "",
  "(A) Relative contribution profiles from the original author-marker NNLS model and the donor-validated positive-marker model, shown for every matched biological line and week. (B) Gene-bootstrap sensitivity intervals for the two nonzero dominant contribution axes. These intervals quantify sensitivity to marker-gene resampling and are not confidence intervals across biological lines. (C) Signed scores for all 18 refined neuronal identities across matched cultures. (D) Five-axis contribution summaries rendered as line plots rather than radar areas; this panel is retained for orientation and is not primary quantitative evidence."
)

writeLines(caption_lines, file.path(run_dir, "FIGURE_CAPTIONS.md"))

readme_lines <- c(
  "# Single-primary TG figure set",
  "",
  "This run replaces the nested four-main-figure package with one source-native, full-page primary figure and four Extended Data figures.",
  "",
  "Scientific allocation:",
  "",
  "1. The primary page follows method order from reference labels to culture interpretation.",
  "2. Broad maturation programs, QC, matrix conditioning, held-out recovery, projection sensitivity, and before/after composition are kept in Extended Data.",
  "3. NNLS is described as using donor-validated positive markers. Anti-markers are used only in the complementary signed score.",
  "4. NNLS outputs are relative transcriptional signature contributions, not literal cell fractions.",
  "",
  "Regenerate with:",
  "",
  "```bash",
  "Rscript scripts/build_single_primary_figure_set.R",
  "```"
)
writeLines(readme_lines, file.path(run_dir, "README.md"))

message("Wrote single-primary figure set to: ", normalizePath(run_dir))
