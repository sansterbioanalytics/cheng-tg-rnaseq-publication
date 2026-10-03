#!/usr/bin/env Rscript

# Assemble a supplemental-only view from pre-existing reference/deconvolution
# outputs. This script imports those artifacts verbatim and does not rerun or
# claim authorship of the upstream single-nucleus or deconvolution analyses.

source("scripts/manuscript_aligned/00_shared.R")

composition_path <- "results/final_article/relative_composition_replicate_summary.csv"
vif_path <- "results/deconv_rerun_expanded/vif_joint_space_comparison.csv"
signed_validation_path <- "results/marker_refinement/refined_reference_signed_score_calls.csv"
bootstrap_path <- "results/deconv_rerun_expanded/weekly_nnls_deconvolution_with_bootstrap_ci.csv"

vapply(
  c(composition_path, vif_path, signed_validation_path, bootstrap_path),
  assert_input,
  character(1)
)

composition <- read_csv(composition_path, show_col_types = FALSE)
vif <- read_csv(vif_path, show_col_types = FALSE)
signed_validation <- read_csv(signed_validation_path, show_col_types = FALSE)

validation_accuracy <- mean(signed_validation$correct_label, na.rm = TRUE)
validation_n <- sum(!is.na(signed_validation$correct_label))
validation_correct <- sum(signed_validation$correct_label, na.rm = TRUE)

composition_palette <- c(
  "A-beta-LTMR-like (NTRK3/PVALB)" = anatomic_okabe_ito[[1]],
  "A-delta-LTMR-like (NTRK2/PIEZO2)" = anatomic_okabe_ito[[2]],
  "C-LTMR-like (CASQ2/P2RY1)" = anatomic_okabe_ito[[3]],
  "Other neuronal identities" = anatomic_okabe_ito[[4]],
  "Non-neuronal / off-target" = anatomic_brand[["grey_mid"]]
)

p_composition <- ggplot(
  composition,
  aes(x = maturation_week, y = mean_value, color = composition_class, fill = composition_class)
) +
  geom_ribbon(aes(ymin = pmax(ci_low, 0), ymax = pmin(ci_high, 1)), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.0) +
  scale_color_manual(values = composition_palette) +
  scale_fill_manual(values = composition_palette) +
  scale_x_continuous(breaks = 1:4) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Imported relative reference-signature contributions",
    subtitle = "Supplemental context only; values are not cell fractions",
    x = "Week after thawing/plating",
    y = "Relative signature contribution",
    color = NULL,
    fill = NULL
  ) +
  manuscript_theme(base_size = 9) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 7.2),
    legend.key.width = grid::unit(0.9, "lines"),
    axis.title.x = element_text(margin = margin(t = 7))
  )

vif_long <- vif |>
  pivot_longer(
    cols = c(vif_original, vif_expanded),
    names_to = "marker_panel",
    values_to = "vif"
  ) |>
  mutate(
    marker_panel = recode(
      marker_panel,
      vif_original = "Original panel",
      vif_expanded = "Refined panel"
    )
  )

vif_summary <- vif_long |>
  group_by(marker_panel) |>
  summarise(
    mean_vif = mean(vif, na.rm = TRUE),
    median_vif = median(vif, na.rm = TRUE),
    max_vif = max(vif, na.rm = TRUE),
    .groups = "drop"
  )
write_manuscript_table(vif_summary, "imported_vif_summary.csv")

p_vif <- ggplot(vif_long, aes(x = marker_panel, y = vif, fill = marker_panel)) +
  geom_hline(yintercept = 10, color = anatomic_brand[["red"]], linetype = "dashed", linewidth = 0.55) +
  geom_boxplot(width = 0.56, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.10, height = 0, size = 1.15, alpha = 0.55, color = anatomic_brand[["ink"]]) +
  scale_fill_manual(values = c("Original panel" = anatomic_brand[["grey"]], "Refined panel" = anatomic_brand[["blue"]])) +
  labs(
    title = "Imported reference-model identifiability checks",
    subtitle = paste0(
      "Signed-marker validation: ", validation_correct, "/", validation_n,
      " correct (", percent(validation_accuracy, accuracy = 0.1),
      "); residual collinearity remains"
    ),
    x = NULL,
    y = "Variance inflation factor",
    fill = NULL,
    caption = "Dashed line: VIF = 10 diagnostic threshold"
  ) +
  manuscript_theme(base_size = 9) +
  coord_cartesian(ylim = c(0, max(vif_long$vif, na.rm = TRUE) * 1.04), clip = "off") +
  theme(
    legend.position = "none",
    plot.caption = element_text(hjust = 0, color = "#55525D", size = 7.2)
  )

supplemental_reference_figure <- (p_composition | p_vif) +
  plot_layout(widths = c(1.2, 0.8)) +
  plot_annotation(
    tag_levels = "A",
    title = "Supplemental reference mapping and model diagnostics",
    subtitle = "All values are imported from the existing validated repository outputs; no upstream model was rerun here.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = anatomic_brand[["ink"]]),
      plot.subtitle = element_text(size = 9.5, color = "#55525D"),
      plot.tag = element_text(face = "bold", size = 12, color = anatomic_brand[["red"]])
    )
  )

save_manuscript_figure(
  supplemental_reference_figure,
  "Figure_S2_imported_reference_mapping_and_diagnostics",
  width = 15,
  height = 8
)

import_provenance <- tribble(
  ~artifact, ~source_path, ~upstream_analysis, ~treatment_here, ~manuscript_location, ~required_language,
  "Relative reference-signature summary", composition_path, "Existing donor-aware reference and NNLS workflow", "Imported and replotted only", "Supplement", "Relative signature contribution; not cell fraction",
  "Joint-space VIF comparison", vif_path, "Existing marker-refinement workflow", "Imported and summarized only", "Supplementary Methods", "Residual collinearity remains",
  "Signed-marker validation calls", signed_validation_path, "Existing donor-aware marker-validation workflow", "Imported and counted only", "Supplementary Methods", "Reference classification validation; not culture identity accuracy",
  "Gene-bootstrap intervals", bootstrap_path, "Existing 500-resample NNLS sensitivity workflow", "Provenance only; not regenerated", "Supplementary Methods/Data", "Sensitivity interval; not biological confidence interval"
) |>
  mutate(source_md5 = vapply(source_path, md5_or_na, character(1)))
write_manuscript_table(import_provenance, "imported_reference_provenance.csv")

writeLines(c(
  "# Supplemental reference/deconvolution methods boundary",
  "",
  "The donor-aware human trigeminal-ganglion reference, marker refinement, signed-score validation, NNLS fitting, VIF calculations, and gene-bootstrap sensitivity analysis were completed in the existing numbered analysis workflow. The manuscript-aligned scripts do not reconstruct those analyses. They import versioned output tables, verify file provenance, and place the resulting material in Supplementary Methods and figures.",
  "",
  paste0("The imported signed-marker classifier correctly labeled ", validation_correct, " of ", validation_n, " donor-level neuronal reference pseudobulks (", percent(validation_accuracy, accuracy = 0.1), "). This validates discrimination inside the reference; it does not prove literal subtype composition in bulk TGN cultures."),
  "",
  "NNLS values are reported as relative reference-signature contributions. They must not be described as cell fractions, percentages of neurons, or histological purity. Gene-bootstrap intervals quantify sensitivity to signature-gene selection and are not biological confidence intervals.",
  "",
  "Primary source: Yang et al., Human and mouse trigeminal ganglia cell atlas implicates multiple cell types in migraine, Neuron (2022), doi:10.1016/j.neuron.2022.03.003."
), file.path(manuscript_text_dir, "SUPPLEMENTAL_REFERENCE_METHODS_BOUNDARY.md"))

message("Imported reference supplement written to: ", normalizePath(manuscript_figure_dir))
