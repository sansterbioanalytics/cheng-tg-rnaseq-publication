#!/usr/bin/env Rscript

# Build the manuscript-aligned RNA-seq visual abstract and literature-curated
# developmental scaffold. Regulatory arrows are curated from cited literature;
# they are not inferred from the bulk RNA-seq data.

source("scripts/manuscript_aligned/00_shared.R")

heatmap_path <- file.path(manuscript_table_dir, "owned_main_heatmap_summary.csv")
trend_path <- file.path(manuscript_table_dir, "owned_direct_gene_blocked_trends.csv")
assert_input(heatmap_path)
assert_input(trend_path)

heatmap_data <- read_csv(heatmap_path, show_col_types = FALSE)
direct_trends <- read_csv(trend_path, show_col_types = FALSE)

grn_nodes <- tribble(
  ~node_id, ~x, ~y, ~node_label, ~node_class,
  "border", 1.0, 0.55, "Cranial ectoderm /\nneural plate border\nTFAP2A", "early competence",
  "placode", 2.1, 0.95, "Preplacodal / cranial\nsensory competence\nSIX1 · EYA1", "cranial sensory",
  "crest", 2.1, 0.15, "Neural-crest-associated\nstate\nSOX10 · FOXD3", "neural crest",
  "neurogenesis", 3.3, 0.55, "Sensory neurogenesis\nNEUROG1 · NEUROD1\nISL1 · POU4F1", "sensory neuron",
  "channel", 4.6, 0.55, "Nociceptor / channel context\nPRDM12 · NTRK1\nSCN9A · SCN10A · TRPV1 · TRPM8", "NaV1.8 context"
)

grn_edges <- tribble(
  ~from, ~to, ~edge_class, ~evidence_scope, ~citation,
  "border", "placode", "lineage transition", "Conceptual cranial sensory lineage transition; not a direct regulatory edge", "Hovland et al. 2020, doi:10.1002/wsbm.1468",
  "border", "crest", "lineage transition", "Conceptual neural crest lineage transition; not a direct regulatory edge", "Hovland et al. 2020, doi:10.1002/wsbm.1468",
  "placode", "neurogenesis", "curated regulation", "SIX1/EYA1-dependent cranial sensory neurogenesis targets include NEUROG1 and POU4F1", "Riddiford and Schlosser 2016, doi:10.7554/eLife.17666",
  "crest", "neurogenesis", "lineage transition", "Neural crest progenitors contribute sensory neurons; not a single direct edge", "Hovland et al. 2022, doi:10.1016/j.devcel.2022.09.006",
  "neurogenesis", "channel", "curated regulation", "PRDM12 is required for the nociceptive lineage and downstream NEUROD1/BRN3A/ISL1 program", "Bartesaghi et al. 2019, doi:10.1016/j.celrep.2019.02.098"
) |>
  left_join(grn_nodes |> select(from = node_id, x_from = x, y_from = y), by = "from") |>
  left_join(grn_nodes |> select(to = node_id, x_to = x, y_to = y), by = "to")

write_manuscript_table(grn_nodes, "literature_curated_grn_nodes.csv")
write_manuscript_table(
  grn_edges |> select(from, to, edge_class, evidence_scope, citation),
  "literature_curated_grn_edges.csv"
)

node_colors <- c(
  "early competence" = "#E9E6E3",
  "cranial sensory" = "#D8EAF4",
  "neural crest" = "#EADBE5",
  "sensory neuron" = "#DDEFE8",
  "NaV1.8 context" = "#F2D9D5"
)

p_grn <- ggplot() +
  geom_curve(
    data = grn_edges,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linetype = edge_class),
    curvature = 0.08,
    color = anatomic_brand[["grey_mid"]],
    linewidth = 0.65,
    arrow = grid::arrow(length = grid::unit(0.08, "inches"), type = "closed")
  ) +
  geom_label(
    data = grn_nodes,
    aes(x = x, y = y, label = node_label, fill = node_class),
    color = anatomic_brand[["ink"]],
    size = 3.0,
    fontface = "bold",
    lineheight = 0.92,
    label.padding = grid::unit(0.25, "lines"),
    label.size = 0.25
  ) +
  scale_fill_manual(values = node_colors, guide = "none") +
  scale_linetype_manual(
    values = c("lineage transition" = "dashed", "curated regulation" = "solid"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(0.55, 5.05), ylim = c(-0.08, 1.22), clip = "off") +
  labs(
    title = "Developmental scaffold",
    subtitle = "Literature-curated lineage and regulatory context; not inferred from bulk RNA-seq",
    caption = "SIX1 and SOX10 capture complementary developmental evidence; neither alone establishes mature trigeminal nociceptor identity."
  ) +
  theme_void(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", color = anatomic_brand[["ink"]], size = 11),
    plot.subtitle = element_text(color = "#55525D", size = 8.5, margin = margin(b = 8)),
    plot.caption = element_text(color = "#55525D", hjust = 0, size = 7.5, margin = margin(t = 8)),
    legend.position = "bottom",
    legend.text = element_text(size = 7.5),
    plot.margin = margin(8, 8, 8, 8)
  )

figure_block_order <- c(
  "Cranial sensory lineage",
  "Sensory neurogenesis",
  "Neuronal architecture",
  "NaV channel context",
  "Sensory context and boundaries"
)
gene_order <- c(
  "SIX1", "EYA1", "SOX10",
  "ISL1", "POU4F1", "PRDM12",
  "RBFOX3", "TUBB3", "SNAP25",
  "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A",
  "TRPV1", "TRPM8", "P2RX3", "NTRK1", "NTRK2", "NTRK3", "PIEZO2", "CALCA", "CALCB"
)

heatmap_plot_data <- heatmap_data |>
  mutate(
    context = factor(context, levels = c(paste0("TGN W", 1:4), "hDRG")),
    figure_block = factor(figure_block, levels = figure_block_order),
    gene = factor(gene, levels = rev(gene_order)),
    label_color = if_else(mean_tmm_log_cpm >= 7.5, "white", anatomic_brand[["ink"]])
  )

p_heatmap <- ggplot(heatmap_plot_data, aes(x = context, y = gene, fill = mean_tmm_log_cpm)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = number(mean_tmm_log_cpm, accuracy = 0.1), color = label_color),
    size = 2.55,
    fontface = "bold"
  ) +
  facet_grid(
    rows = vars(figure_block),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_gradient(
    low = anatomic_brand[["grey_light"]],
    high = anatomic_brand[["blue"]],
    name = "Mean TMM\nlog-CPM"
  ) +
  scale_color_identity() +
  labs(
    title = "Owned RNA-seq: target context and identity boundaries",
    subtitle = "Three ANAT002 differentiation lots per week; hDRG is a pooled-tissue expression comparator",
    x = NULL,
    y = NULL,
    caption = "Low peptidergic markers and high TRPM8 remain visible by design; the heatmap is contextual rather than a purity or cell-fraction assay."
  ) +
  manuscript_theme(base_size = 8.5) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic", size = 7.2),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 7.2, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

evidence_stages <- tribble(
  ~stage, ~window, ~heading, ~result, ~interpretation,
  1, "W1-W4 after thaw", "RNA-seq", "SCN10A and a sensory-neuronal program are present", "Developmental and molecular context",
  2, "6-7 weeks after plating", "Early function", "5/15 NaV1.8-positive; 2/15 above 1 nA", "Functional emergence",
  3, "10-11 weeks after plating", "Robust function", "31/42 NaV1.8-positive; 20/42 above 1 nA", "Reproducible late function across two differentiations",
  4, "10-13 weeks after plating", "Target engagement and protein", "30 nM VX-548 sensitivity; representative NaV1.8 immunoreactivity", "Pharmacological and protein corroboration"
)

p_bridge <- ggplot(evidence_stages, aes(x = stage, y = 0.5)) +
  annotate(
    "segment", x = 1, xend = 4, y = 0.5, yend = 0.5,
    color = anatomic_brand[["grey"]], linewidth = 1.1,
    arrow = grid::arrow(length = grid::unit(0.1, "inches"), type = "closed")
  ) +
  geom_point(
    aes(color = heading),
    size = 5
  ) +
  geom_text(aes(y = 0.87, label = heading), fontface = "bold", size = 3.1) +
  geom_text(aes(y = 0.70, label = window), size = 2.75) +
  geom_text(aes(y = 0.15, label = str_wrap(result, width = 31)), size = 2.75, lineheight = 0.9) +
  geom_text(
    aes(y = -0.42, label = str_wrap(interpretation, width = 31)),
    size = 2.55,
    lineheight = 0.9,
    color = "#55525D"
  ) +
  scale_color_manual(
    values = c(
      "RNA-seq" = anatomic_brand[["blue"]],
      "Early function" = anatomic_brand[["red"]],
      "Robust function" = anatomic_brand[["red_dark"]],
      "Target engagement and protein" = anatomic_okabe_ito[[4]]
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0.6, 4.4), ylim = c(-0.68, 1.03), clip = "off") +
  labs(
    title = "Evidence chain",
    subtitle = "The molecular program is established early; robust, pharmacologically tractable NaV1.8 emerges with prolonged culture",
    caption = "RNA-seq and functional assays are staged measurements, not synchronized same-timepoint validation."
  ) +
  theme_void(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 11, color = anatomic_brand[["ink"]]),
    plot.subtitle = element_text(size = 8.5, color = "#55525D", margin = margin(b = 8)),
    plot.caption = element_text(size = 7.5, color = "#55525D", hjust = 0, margin = margin(t = 8)),
    plot.margin = margin(8, 8, 8, 8)
  )

figure_2_visual_abstract <- (p_grn | p_heatmap) / p_bridge +
  plot_layout(widths = c(0.92, 1.08), heights = c(1.22, 0.78)) +
  plot_annotation(
    title = "Figure 2C. Developmental and molecular context for a maturing NaV1.8 platform",
    subtitle = "RNA-seq was performed on three independent differentiation lots from one ANAT002 genotype during weeks 1-4 after thawing.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, color = anatomic_brand[["ink"]]),
      plot.subtitle = element_text(size = 10, color = "#55525D")
    )
  )

save_manuscript_figure(
  figure_2_visual_abstract,
  "Figure_2C_RNAseq_visual_abstract_manuscript_aligned",
  width = 16,
  height = 11.5
)

standalone_grn <- p_grn +
  plot_annotation(
    title = "Literature-curated developmental scaffold for the TGN platform",
    subtitle = "A conceptual gene-regulatory and lineage model; measured expression does not establish the direction of these edges.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = anatomic_brand[["ink"]]),
      plot.subtitle = element_text(size = 9.5, color = "#55525D")
    )
  )
save_manuscript_figure(
  standalone_grn,
  "Figure_S1_literature_curated_lineage_GRN",
  width = 13,
  height = 6.5
)

scn10a_trend <- direct_trends |> filter(gene == "SCN10A")
sox10_trend <- direct_trends |> filter(gene == "SOX10")

caption_lines <- c(
  "# Manuscript-aligned figure captions",
  "",
  "## Figure 2C. Developmental and molecular context for a maturing NaV1.8 platform",
  "",
  "The left schematic is a literature-curated developmental scaffold connecting cranial sensory competence, neural-crest-associated state, sensory neurogenesis, and nociceptor/channel context. Dashed arrows denote conceptual lineage transitions; solid arrows denote curated regulatory relationships. The network is not inferred from the bulk RNA-seq data. The heatmap reports TMM-normalized log2 counts per million for lineage, neuronal, sodium-channel, and sensory-context genes. TGN values are means across three independent ANAT002 differentiation lots at each week after thawing; hDRG is a pooled-tissue expression comparator. The lower evidence chain places the early RNA-seq window in temporal relation to later functional, pharmacological, and protein measurements.",
  "",
  paste0(
    "SCN10A was detected throughout weeks 1-4 but did not show a conclusive linear increase (blocked slope ",
    number(scn10a_trend$blocked_slope_per_week, accuracy = 0.001),
    " log-CPM/week; 95% CI ",
    number(scn10a_trend$blocked_slope_ci_low, accuracy = 0.001),
    " to ",
    number(scn10a_trend$blocked_slope_ci_high, accuracy = 0.001),
    "; nominal p=", number(scn10a_trend$blocked_slope_p_value, accuracy = 0.001), ")."
  ),
  "",
  "## Figure S1. Literature-curated developmental scaffold",
  "",
  "Conceptual lineage and regulatory relationships relevant to cranial sensory-neuronal and nociceptor development. This figure visualizes a literature model and overlays no inferred regulatory edges from the present bulk RNA-seq experiment.",
  "",
  paste0(
    "SOX10 changed across the four-week RNA window in the blocked model (slope ",
    number(sox10_trend$blocked_slope_per_week, accuracy = 0.001),
    "; nominal p=", number(sox10_trend$blocked_slope_p_value, accuracy = 0.001),
    "). This targeted result should receive one bounded sentence, not a broad maturation interpretation."
  )
)
writeLines(caption_lines, file.path(manuscript_text_dir, "FIGURE_CAPTIONS.md"))

message("Manuscript-aligned visual abstract written to: ", normalizePath(manuscript_figure_dir))
