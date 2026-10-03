#!/usr/bin/env Rscript

# Build the final, direct-gene Cheng et al. insertion packet outputs.
# This script owns only the matched TGN bulk RNA-seq summaries and figures.

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
  library(png)
})

source("R/anatomic_theme.R")

out_dir <- file.path("publication_runs", "2026-08-27-cheng-figure2-powerpoint")
fig_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
asset_dir <- file.path(out_dir, "source_assets")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- "RealTGN.salmon.merged.gene_counts.tsv"
source_genotype <- "ANAT002"
source_sex <- "male"
lots <- c("TGN1003A", "TGN1028B", "TGN1215AvX")

gene_manifest <- tribble(
  ~gene, ~figure_block, ~manuscript_role,
  "SIX1", "Cranial sensory lineage", "trigeminal-lineage-associated developmental context",
  "EYA1", "Cranial sensory lineage", "SIX1 co-regulatory cranial sensory competence",
  "SOX10", "Cranial sensory lineage", "neural-crest-associated state and developmental heterogeneity",
  "ISL1", "Sensory neurogenesis", "pan-sensory neuronal differentiation",
  "POU4F1", "Sensory neurogenesis", "BRN3A; pan-sensory neuronal differentiation",
  "PRDM12", "Sensory neurogenesis", "nociceptor-lineage regulatory context",
  "RBFOX3", "Neuronal architecture", "neuronal identity",
  "TUBB3", "Neuronal architecture", "neuronal cytoskeletal identity",
  "SNAP25", "Neuronal architecture", "synaptic architecture",
  "SCN5A", "NaV channel context", "TTX-resistant developmental comparator NaV1.5",
  "SCN8A", "NaV channel context", "NaV1.6 comparator",
  "SCN9A", "NaV channel context", "sensory sodium channel NaV1.7",
  "SCN10A", "NaV channel context", "target-channel molecular context NaV1.8",
  "SCN11A", "NaV channel context", "sensory sodium channel NaV1.9",
  "TRPV1", "Sensory context and boundaries", "heat/nociceptor-associated context",
  "TRPM8", "Sensory context and boundaries", "cold-sensory context",
  "P2RX3", "Sensory context and boundaries", "sensory/nociceptor-associated context",
  "NTRK1", "Sensory context and boundaries", "TRKA nociceptor-lineage context",
  "NTRK2", "Sensory context and boundaries", "TRKB sensory context",
  "NTRK3", "Sensory context and boundaries", "TRKC sensory context",
  "PIEZO2", "Sensory context and boundaries", "mechanosensory context",
  "CALCA", "Sensory context and boundaries", "peptidergic boundary",
  "CALCB", "Sensory context and boundaries", "peptidergic boundary"
)

imported <- read_tsv(counts_path, show_col_types = FALSE, name_repair = "minimal")
stopifnot(all(c("gene_id", "gene_name") %in% colnames(imported)))
sample_ids <- setdiff(colnames(imported), c("gene_id", "gene_name"))
counts <- imported |>
  transmute(gene = toupper(trimws(gene_name)), across(all_of(sample_ids), as.numeric)) |>
  filter(!is.na(gene), gene != "") |>
  group_by(gene) |>
  summarise(across(all_of(sample_ids), ~ sum(.x, na.rm = TRUE)), .groups = "drop") |>
  column_to_rownames("gene") |>
  as.matrix()
storage.mode(counts) <- "numeric"
stopifnot(!anyDuplicated(rownames(counts)))

total_counts <- colSums(counts)
detected_genes <- colSums(counts > 0)
metadata <- tibble(sample_id = colnames(counts)) |>
  mutate(
    sample_label = sub("^[^_]+_", "", sample_id),
    sample_class = case_when(
      grepl("^hDRG", sample_label) ~ "hDRG comparator",
      grepl("^TGN", sample_label) ~ "TGN culture",
      grepl("^TAK", sample_label) ~ "TAK control",
      TRUE ~ "Other"
    ),
    differentiation_lot = case_when(
      grepl("^TGN.*W[0-9]+$", sample_label) ~ sub("W[0-9]+$", "", sample_label),
      TRUE ~ NA_character_
    ),
    rna_week_after_thaw = as.integer(str_remove(str_extract(sample_label, "W[0-9]+$"), "^W")),
    source_hiPSC_genotype = if_else(sample_class == "TGN culture", source_genotype, NA_character_),
    source_hiPSC_sex = if_else(sample_class == "TGN culture", source_sex, NA_character_),
    total_counts = unname(total_counts[sample_id]),
    detected_genes = unname(detected_genes[sample_id])
  )

matched_metadata <- metadata |>
  filter(sample_class == "TGN culture", differentiation_lot %in% lots, rna_week_after_thaw %in% 1:4) |>
  arrange(differentiation_lot, rna_week_after_thaw)
hdrg_metadata <- metadata |> filter(sample_class == "hDRG comparator")
stopifnot(
  nrow(matched_metadata) == 12L,
  n_distinct(matched_metadata$differentiation_lot) == 3L,
  all(table(matched_metadata$differentiation_lot) == 4L),
  all(matched_metadata$source_hiPSC_genotype == source_genotype),
  length(hdrg_metadata$sample_id) > 0L,
  all(gene_manifest$gene %in% rownames(counts))
)

tmm_logcpm <- function(x) {
  dge <- DGEList(counts = x) |> calcNormFactors(method = "TMM")
  cpm(dge, log = TRUE, prior.count = 2)
}

context_metadata <- bind_rows(matched_metadata, hdrg_metadata)
context_counts <- counts[, context_metadata$sample_id, drop = FALSE]
context_logcpm <- tmm_logcpm(context_counts)
expression_long <- as.data.frame(context_logcpm[gene_manifest$gene, , drop = FALSE]) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "tmm_log_cpm") |>
  left_join(context_metadata, by = "sample_id") |>
  left_join(gene_manifest, by = "gene") |>
  mutate(
    context = case_when(
      sample_class == "TGN culture" ~ paste0("TGN W", rna_week_after_thaw),
      sample_class == "hDRG comparator" ~ "hDRG",
      TRUE ~ sample_class
    )
  )

tgn_expression <- expression_long |> filter(sample_class == "TGN culture", differentiation_lot %in% lots, rna_week_after_thaw %in% 1:4)
tgn_week_summary <- tgn_expression |>
  group_by(gene, figure_block, manuscript_role, rna_week_after_thaw) |>
  summarise(
    mean_tmm_log_cpm = mean(tmm_log_cpm),
    sd_tmm_log_cpm = sd(tmm_log_cpm),
    n_differentiation_lots = n_distinct(differentiation_lot),
    .groups = "drop"
  ) |>
  mutate(context = paste0("TGN W", rna_week_after_thaw), n_hDRG_libraries = NA_integer_)
hdrg_summary <- expression_long |>
  filter(sample_class == "hDRG comparator") |>
  group_by(gene, figure_block, manuscript_role) |>
  summarise(
    mean_tmm_log_cpm = mean(tmm_log_cpm),
    sd_tmm_log_cpm = sd(tmm_log_cpm),
    n_differentiation_lots = NA_integer_,
    n_hDRG_libraries = n(),
    .groups = "drop"
  ) |>
  mutate(context = "hDRG")
heatmap_summary <- bind_rows(tgn_week_summary, hdrg_summary) |>
  select(gene, figure_block, manuscript_role, context, mean_tmm_log_cpm, sd_tmm_log_cpm, n_differentiation_lots, n_hDRG_libraries)

blocked_trend <- function(df) {
  fit <- lm(tmm_log_cpm ~ rna_week_after_thaw + differentiation_lot, data = df)
  co <- summary(fit)$coefficients["rna_week_after_thaw", , drop = FALSE]
  df_resid <- df.residual(fit)
  est <- unname(co[1, "Estimate"])
  se <- unname(co[1, "Std. Error"])
  tibble(
    blocked_slope_per_week = est,
    blocked_slope_ci_low = est - qt(0.975, df_resid) * se,
    blocked_slope_ci_high = est + qt(0.975, df_resid) * se,
    blocked_slope_p_value = unname(co[1, "Pr(>|t|)"]),
    n_differentiation_lots = n_distinct(df$differentiation_lot),
    n_libraries = nrow(df)
  )
}
trends <- tgn_expression |>
  group_by(gene, figure_block, manuscript_role) |>
  group_modify(~ blocked_trend(.x)) |>
  ungroup() |>
  mutate(blocked_slope_fdr_bh = p.adjust(blocked_slope_p_value, method = "BH"))

write_csv(metadata, file.path(table_dir, "owned_sample_metadata.csv"), na = "NA")
write_csv(gene_manifest |> mutate(present_in_owned_counts = gene %in% rownames(counts)), file.path(table_dir, "main_figure_gene_manifest.csv"), na = "NA")
write_csv(expression_long, file.path(table_dir, "owned_23_gene_expression_logcpm.csv"), na = "NA")
write_csv(heatmap_summary, file.path(table_dir, "owned_23_gene_week_summary.csv"), na = "NA")
write_csv(trends, file.path(table_dir, "owned_23_gene_blocked_trends.csv"), na = "NA")

gene_order <- gene_manifest$gene
block_order <- unique(gene_manifest$figure_block)
heat_plot_data <- heatmap_summary |>
  mutate(
    context = factor(context, levels = c(paste0("TGN W", 1:4), "hDRG")),
    gene = factor(gene, levels = rev(gene_order)),
    figure_block = factor(figure_block, levels = block_order),
    label_color = if_else(mean_tmm_log_cpm >= 7.5, "white", "#2F2E2E")
  )

p_heat <- ggplot(heat_plot_data, aes(context, gene, fill = mean_tmm_log_cpm)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = number(mean_tmm_log_cpm, accuracy = 0.1), color = label_color), size = 2.1, fontface = "bold") +
  facet_grid(rows = vars(figure_block), scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient(low = "#F2F2F2", high = "#1677B8", name = "TMM\nlog-CPM") +
  scale_color_identity() +
  labs(title = "23-gene sensory-neuronal signature", subtitle = "Means across three matched ANAT002 lots; hDRG is contextual", x = NULL, y = NULL) +
  theme_anatomic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 7.2, color = "#55525D"),
    axis.text.x = element_text(face = "bold", size = 7), axis.text.y = element_text(face = "italic", size = 6.8),
    strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, face = "bold", size = 6.7),
    panel.grid = element_blank(), legend.position = "right", plot.margin = margin(5, 4, 5, 4)
  )

trend_plot_data <- trends |>
  mutate(gene = factor(gene, levels = rev(gene_order)), figure_block = factor(figure_block, levels = block_order))
p_trend <- ggplot(trend_plot_data, aes(blocked_slope_per_week, gene, color = figure_block)) +
  geom_vline(xintercept = 0, color = "#8C8C8C", linewidth = 0.4) +
  geom_segment(aes(x = blocked_slope_ci_low, xend = blocked_slope_ci_high, y = gene, yend = gene), linewidth = 0.55) +
  geom_point(size = 1.7) +
  facet_grid(rows = vars(figure_block), scales = "free_y", space = "free_y", switch = "y") +
  scale_color_manual(values = c("#4E79A7", "#59A14F", "#9C755F", "#E15759", "#B07AA1"), guide = "none") +
  labs(title = "Blocked temporal remodeling", subtitle = "Slope per week (95% CI); no significance calls", x = "log-CPM change per week", y = NULL) +
  theme_anatomic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 7.2, color = "#55525D"),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(),
    panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.25), plot.margin = margin(5, 5, 5, 2)
  )

ggsave(file.path(fig_dir, "Figure_2_RNAseq_heatmap.png"), p_heat, width = 7.2, height = 8.2, dpi = 400, bg = "white")
ggsave(file.path(fig_dir, "Figure_2_RNAseq_heatmap.pdf"), p_heat, width = 7.2, height = 8.2, bg = "white")
ggsave(file.path(fig_dir, "Figure_2_RNAseq_slope.png"), p_trend, width = 4.8, height = 8.2, dpi = 400, bg = "white")
ggsave(file.path(fig_dir, "Figure_2_RNAseq_slope.pdf"), p_trend, width = 4.8, height = 8.2, bg = "white")

read_asset <- function(path) {
  if (!file.exists(path)) return(NULL)
  grid::rasterGrob(readPNG(path, native = FALSE), interpolate = TRUE)
}
asset_plot <- function(path, label) {
  ggplot() + annotation_custom(read_asset(path), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.1, vjust = 1.25, fontface = "bold", size = 5, color = "#2F2E2E") +
    coord_cartesian(expand = FALSE) + theme_void() + theme(plot.margin = margin(1, 1, 1, 1))
}
phase_paths <- c(file.path(asset_dir, "source_phase_D7.png"), file.path(asset_dir, "source_phase_D28.png"))
ihc_paths <- c(file.path(asset_dir, "source_IHC_ISL1_TUJ1_BRN3A.png"), file.path(asset_dir, "source_IHC_SIX1_TUJ1.png"))


sample_plot_data <- tgn_expression |>
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    differentiation_lot = factor(differentiation_lot, levels = lots),
    rna_week_after_thaw = factor(rna_week_after_thaw, levels = 1:4)
  )
sample_heat <- sample_plot_data |>
  group_by(gene) |>
  mutate(z = as.numeric(scale(tmm_log_cpm))) |>
  ungroup() |>
  mutate(sample_key = factor(paste0(differentiation_lot, " W", rna_week_after_thaw), levels = as.vector(outer(lots, 1:4, paste, sep = " W"))))
p_s2_heat <- ggplot(sample_heat, aes(sample_key, gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(rows = vars(figure_block), scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D7301F", midpoint = 0, name = "Gene-centered\nz-score") +
  labs(title = "Supplementary Figure S2A. Lot- and week-resolved signature", subtitle = "Each gene is centered across the 12 matched TGN libraries", x = NULL, y = NULL) +
  theme_anatomic(base_size = 8) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), axis.text.y = element_text(size = 6.8), strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, face = "bold", size = 6.5), panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 7.2, color = "#55525D"), plot.margin = margin(5, 5, 5, 5))
p_s2_lines <- ggplot(sample_plot_data, aes(as.numeric(as.character(rna_week_after_thaw)), tmm_log_cpm, group = differentiation_lot, color = differentiation_lot)) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.1) +
  facet_wrap(~ gene, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = 1:4) + scale_color_manual(values = c("#0072B2", "#009E73", "#CC79A7"), name = "Lot") +
  labs(title = "Supplementary Figure S2B. Per-lot expression trajectories", subtitle = "Raw TMM log-CPM for all 23 prespecified genes", x = "Week after thawing", y = "TMM log-CPM") +
  theme_anatomic(base_size = 8) + theme(strip.text = element_text(face = "italic", size = 7), plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 7.2, color = "#55525D"), legend.position = "bottom", panel.grid.minor = element_blank(), plot.margin = margin(5, 5, 5, 5))
ggsave(file.path(fig_dir, "Figure_S2_direct_gene_heatmap.png"), p_s2_heat, width = 11, height = 8.3, dpi = 350, bg = "white")
ggsave(file.path(fig_dir, "Figure_S2_direct_gene_heatmap.pdf"), p_s2_heat, width = 11, height = 8.3, bg = "white")
ggsave(file.path(fig_dir, "Figure_S2_direct_gene_trajectories.png"), p_s2_lines, width = 11, height = 9.5, dpi = 350, bg = "white")
ggsave(file.path(fig_dir, "Figure_S2_direct_gene_trajectories.pdf"), p_s2_lines, width = 11, height = 9.5, bg = "white")

writeLines(c(
  "# Final Cheng et al. RNA-seq insertion packet",
  "",
  "This packet is a separate writing and figure deliverable; the Cheng manuscript is not edited.",
  "",
  "## Placement",
  "",
  "Replace the current Figure 2 transcriptomic subsection and RNA-seq Methods placeholder with the supplied Results and Methods text. Use the Figure 2 and Figure S2 legends verbatim unless author-supplied image timing requires a wording correction.",
  "",
  "## Provenance",
  "",
  paste0("Counts source: ", counts_path, ". Gene-level counts were upper-cased, duplicate symbols summed, and normalized with edgeR TMM; expression is reported as log2 counts per million with prior count 2."),
  "The longitudinal unit is one library per lot per week for three ANAT002 lots (12 libraries). hDRG is contextual only and is not part of the blocked models.",
  "",
  "## Scope boundary",
  "",
  "This packet contains direct gene-expression evidence only. Reference-subtype and composition analyses are not part of this submission story."
), file.path(out_dir, "PLACEMENT_AND_PROVENANCE.md"))

message("Final analysis and figures written to ", normalizePath(out_dir))
