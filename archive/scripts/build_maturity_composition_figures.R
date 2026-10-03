#!/usr/bin/env Rscript

# Assemble high-impact manuscript figures from the validated review panels.
# This script does not recompute the RNA-seq analysis. It adds consistent
# panel labels and a restrained journal-style hierarchy to the existing plots.

suppressPackageStartupMessages({
  library(grid)
  library(png)
})

run_dir <- file.path("publication_runs", "2026-07-24")
figure_dir <- file.path(run_dir, "figures")
output_dir <- file.path(run_dir, "assembled_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_path <- function(notebook_dir, filename) {
  path <- file.path(figure_dir, notebook_dir, filename)
  if (!file.exists(path)) stop("Missing panel: ", path, call. = FALSE)
  path
}

read_panel <- function(path) {
  rasterGrob(readPNG(path), interpolate = TRUE)
}

draw_panel <- function(panel, label, descriptor, row, col, layout, label_size = 14) {
  pushViewport(viewport(
    layout.pos.row = row, layout.pos.col = col,
    layout = layout
  ))
  grid.rect(gp = gpar(fill = "white", col = "#D9D6D5", lwd = 0.7))
  pushViewport(viewport(
    x = unit(0.5, "npc"), y = unit(0.48, "npc"),
    width = unit(0.96, "npc"), height = unit(0.86, "npc"),
    just = c("center", "center")
  ))
  grid.draw(panel)
  popViewport()
  grid.roundrect(
    x = unit(0.055, "npc"), y = unit(0.94, "npc"),
    width = unit(0.09, "npc"), height = unit(0.09, "npc"),
    r = unit(0.018, "snpc"),
    gp = gpar(fill = "#A9252A", col = NA)
  )
  grid.text(
    label, x = unit(0.055, "npc"), y = unit(0.94, "npc"),
    gp = gpar(col = "white", fontsize = label_size, fontface = "bold")
  )
  grid.text(
    descriptor, x = unit(0.12, "npc"), y = unit(0.94, "npc"),
    just = "left", gp = gpar(col = "#242323", fontsize = 10.5, fontface = "bold")
  )
  popViewport()
}

draw_figure <- function(filename, title, subtitle, panels, nrow, ncol,
                        width, height, widths = NULL, heights = NULL) {
  if (is.null(widths)) widths <- rep(1, ncol)
  if (is.null(heights)) heights <- rep(1, nrow)
  output_png <- file.path(output_dir, paste0(filename, ".png"))
  output_pdf <- file.path(output_dir, paste0(filename, ".pdf"))

  draw <- function() {
    grid.newpage()
    page_layout <- grid.layout(
      nrow = nrow + 1, ncol = ncol,
      heights = c(unit(0.8, "in"), unit(heights, "null")),
      widths = unit(widths, "null")
    )
    pushViewport(viewport(layout = page_layout))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:ncol))
    grid.text(title, x = unit(0, "npc"), y = unit(0.78, "npc"), just = "left",
              gp = gpar(fontsize = 18, fontface = "bold", col = "#242323"))
    grid.text(subtitle, x = unit(0, "npc"), y = unit(0.30, "npc"), just = "left",
              gp = gpar(fontsize = 10.5, col = "#4E4B4B"))
    popViewport()

    panel_layout <- grid.layout(nrow, ncol, widths = unit(widths, "null"),
                                heights = unit(heights, "null"))
    for (panel in panels) {
      draw_panel(
        panel = panel$image, label = panel$label, descriptor = panel$descriptor,
        row = panel$row + 1, col = panel$col, layout = page_layout
      )
    }
    popViewport()
  }

  png(output_png, width = width * 300, height = height * 300, res = 300, type = "cairo")
  draw()
  dev.off()
  pdf(output_pdf, width = width, height = height, useDingbats = FALSE)
  draw()
  dev.off()
  invisible(c(png = output_png, pdf = output_pdf))
}

# Figure 1: compact calibration figure. Suitable as main Figure 1 or Extended Data
# Figure 1 when word count requires a biology-first opening figure.
figure_1_panels <- list(
  list(
    image = read_panel(panel_path("07_TG_consensus_identity_refinement", "038_fig-consensus-neuronal-umap-1.png")),
    label = "A", descriptor = "Donor-supported neuronal reference", row = 1, col = 1
  ),
  list(
    image = read_panel(panel_path("01_TG_sc_reprocessing", "014_combined-pseudobulk-marker-heatmap-1.png")),
    label = "B", descriptor = "Gene programs define neuronal and off-target states", row = 1, col = 2
  ),
  list(
    image = read_panel(panel_path("05_TG_single_cell_marker_refinement", "032_signed-marker-reference-validation-1.png")),
    label = "C", descriptor = "Donor-held-out signed-marker validation", row = 2, col = 1
  ),
  list(
    image = read_panel(panel_path("05_TG_single_cell_marker_refinement", "033_fig-refit-vif-comparison-1.png")),
    label = "D", descriptor = "Refined markers improve subtype identifiability", row = 2, col = 2
  )
)
draw_figure(
  "Figure_1_reference_and_marker_calibration",
  "A donor-aware TG reference supports gene-grounded culture-state inference",
  "Consensus identities are supported across donors; direct marker programs and signed scores define the readout used for maturity and composition.",
  figure_1_panels, nrow = 2, ncol = 2, width = 14, height = 10
)

# Figure 2: primary maturity figure. The narrative follows gene programs ->
# aggregate maturity -> fixed-reference trajectory, all with specimen identity retained.
figure_2_panels <- list(
  list(
    image = read_panel(panel_path("03_TG_incoming_sample_QC_maturation", "022_fig-composite-signature-trajectory-1.png")),
    label = "A", descriptor = "Gene-program trajectories across culture weeks", row = 1, col = 1
  ),
  list(
    image = read_panel(panel_path("03_TG_incoming_sample_QC_maturation", "023_fig-aggregate-neuronal-index-1.png")),
    label = "B", descriptor = "Aggregate neuronal identity index", row = 1, col = 2
  ),
  list(
    image = read_panel(panel_path("03_TG_incoming_sample_QC_maturation", "024_fig-pca-projection-trajectory-1.png")),
    label = "C", descriptor = "Fixed-reference neuronal trajectory", row = 2, col = 1
  ),
  list(
    image = read_panel(panel_path("07_TG_consensus_identity_refinement", "041_fig-selected-pc-trajectory-comparison-1.png")),
    label = "D", descriptor = "Trajectory pattern persists across embedding views", row = 2, col = 2
  )
)
draw_figure(
  "Figure_2_reference_anchored_neuronal_maturity",
  "TG cultures acquire a sensory-neuronal transcriptomic state across the time course",
  "Maturity is supported by coordinated gene programs and fixed-reference projection; nonlinear embeddings are sensitivity views, not the quantitative endpoint.",
  figure_2_panels, nrow = 2, ncol = 2, width = 16, height = 11
)

# Figure 3: primary composition figure. This deliberately combines a direct
# signed-score readout, uncertainty-aware NNLS estimates, and the final concise
# identity summary; conditioning is included so the limitations travel with the result.
figure_3_panels <- list(
  list(
    image = read_panel(panel_path("06_TG_deconvolution_rerun_expanded_markers", "034_signed-refined-identity-trajectories-1.png")),
    label = "A", descriptor = "Signed identity evidence from positive and anti-markers", row = 1, col = 1
  ),
  list(
    image = read_panel(panel_path("06_TG_deconvolution_rerun_expanded_markers", "037_fig-composition-ci-trajectory-1.png")),
    label = "B", descriptor = "Relative contributions with gene-bootstrap intervals", row = 1, col = 2
  ),
  list(
    image = read_panel(panel_path("07_TG_consensus_identity_refinement", "043_fig-radar-identity-distribution-1.png")),
    label = "C", descriptor = "Accessible summary of culture-associated identities", row = 2, col = 1
  ),
  list(
    image = read_panel(panel_path("06_TG_deconvolution_rerun_expanded_markers", "035_fig-vif-recheck-1.png")),
    label = "D", descriptor = "Improved but not eliminated signature collinearity", row = 2, col = 2
  )
)
draw_figure(
  "Figure_3_estimated_relative_composition",
  "Cultures show shifting sensory-identity signatures with uncertainty-aware composition estimates",
  "Composition is interpreted as a relative transcriptional signature contribution, supported by signed scores and bounded by gene-bootstrap sensitivity intervals.",
  figure_3_panels, nrow = 2, ncol = 2, width = 16, height = 11
)

message("Wrote assembled figures to: ", normalizePath(output_dir))
