#!/usr/bin/env Rscript

# Assemble the selected index-gallery validation panels into an explicit
# methods-oriented Extended Data figure for the refined publication run.

suppressPackageStartupMessages({
  library(grid)
  library(png)
})

source("R/anatomic_theme.R")

publication_run <- Sys.getenv("PUBLICATION_RUN", "2026-08-04")
figure_dir <- file.path("publication_runs", publication_run, "final_article_figures")
output_dir <- "results/final_article"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_paths <- c(
  "publication_runs/2026-07-24/figures/02_TG_reviewer_figures/019_fig-bulk-neuron-only-projection-1.png",
  "publication_runs/2026-07-24/figures/05_TG_single_cell_marker_refinement/032_signed-marker-reference-validation-1.png",
  "publication_runs/2026-07-24/figures/05_TG_single_cell_marker_refinement/033_fig-refit-vif-comparison-1.png",
  "publication_runs/2026-07-24/figures/04_TG_deconvolution_diagnostics/027_donor-held-out-recovery-test-1.png"
)
if (!all(file.exists(panel_paths))) stop("A gallery panel required for Extended Data Figure 3 is missing.")

read_panel <- function(path) rasterGrob(readPNG(path), interpolate = TRUE)

draw_panel <- function(panel, label, descriptor, row, col, layout) {
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = col, layout = layout))
  grid.rect(gp = gpar(fill = "white", col = "#D9D6D5", lwd = 0.7))
  pushViewport(viewport(
    x = unit(0.5, "npc"), y = unit(0.46, "npc"),
    width = unit(0.95, "npc"), height = unit(0.82, "npc")
  ))
  grid.draw(panel)
  popViewport()
  grid.roundrect(
    x = unit(0.055, "npc"), y = unit(0.94, "npc"),
    width = unit(0.09, "npc"), height = unit(0.10, "npc"),
    r = unit(0.018, "snpc"),
    gp = gpar(fill = anatomic_brand[["red"]], col = NA)
  )
  grid.text(
    label, x = unit(0.055, "npc"), y = unit(0.94, "npc"),
    gp = gpar(col = "white", fontsize = 14, fontface = "bold")
  )
  grid.text(
    descriptor, x = unit(0.12, "npc"), y = unit(0.94, "npc"), just = "left",
    gp = gpar(col = anatomic_brand[["ink"]], fontsize = 10.5, fontface = "bold")
  )
  popViewport()
}

panels <- lapply(panel_paths, read_panel)
draw <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 3, ncol = 2,
    heights = unit(c(0.75, 1, 1), c("in", "null", "null"))
  )
  pushViewport(viewport(layout = layout))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid.text(
    "Extended validation of reference mapping and marker-matrix refitting",
    x = unit(0, "npc"), y = unit(0.78, "npc"), just = "left",
    gp = gpar(fontsize = 18, fontface = "bold", col = anatomic_brand[["ink"]])
  )
  grid.text(
    "Selected index-gallery panels show the reference projection, donor-aware signed-marker selection, refit conditioning, and donor-held-out recovery test.",
    x = unit(0, "npc"), y = unit(0.30, "npc"), just = "left",
    gp = gpar(fontsize = 10.5, col = "#4E4B4B")
  )
  popViewport()
  draw_panel(panels[[1]], "A", "Neuron-only reference projection", 2, 1, layout)
  draw_panel(panels[[2]], "B", "Donor-held-out signed-marker validation", 2, 2, layout)
  draw_panel(panels[[3]], "C", "Positive-plus-anti-marker refit reduces VIF", 3, 1, layout)
  draw_panel(panels[[4]], "D", "Held-out donor mixture recovery", 3, 2, layout)
  popViewport()
}

png_path <- file.path(figure_dir, "Extended_Data_Figure_3_method_validation_gallery.png")
pdf_path <- file.path(figure_dir, "Extended_Data_Figure_3_method_validation_gallery.pdf")
png(png_path, width = 16 * 200, height = 12 * 200, res = 200, type = "cairo")
draw()
dev.off()
pdf(pdf_path, width = 16, height = 12, useDingbats = FALSE)
draw()
dev.off()

writeLines(c(
  "# Extended Data Figure 3", "",
  "Method validation for the human TG reference and marker-matrix refit.", "",
  "A. Incoming cultures projected into the neuron-only human TG reference space, providing the reference-mapping context used for culture-state interpretation.",
  "B. Donor-held-out signed-marker validation: positive markers are required to remain enriched and anti-markers depleted relative to nearby alternative identities in held-out donor pseudobulk profiles.",
  "C. Refit VIF comparison showing the gain in signature-matrix conditioning after the positive-plus-anti-marker panel is constructed.",
  "D. Donor-held-out pure and nearest-neighbor 50:50 mixture recovery; the dashed line marks maximum absolute fraction error of 0.10. This test evaluates donor generalization rather than self-reconstruction."
), file.path(output_dir, "EXTENDED_DATA_FIGURE_3_CAPTION.md"))

message("Rendered Extended Data Figure 3 for publication run: ", publication_run)
