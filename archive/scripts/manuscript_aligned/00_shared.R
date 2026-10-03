# Shared configuration and helpers for the manuscript-aligned NaV1.8 platform analysis.

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

manuscript_run_date <- "2026-08-23"
manuscript_run_id <- paste0(manuscript_run_date, "-manuscript-aligned-nav18-platform")
manuscript_output_dir <- file.path("publication_runs", manuscript_run_id)
manuscript_figure_dir <- file.path(manuscript_output_dir, "figures")
manuscript_table_dir <- file.path(manuscript_output_dir, "tables")
manuscript_text_dir <- file.path(manuscript_output_dir, "manuscript_text")

dir.create(manuscript_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_text_dir, recursive = TRUE, showWarnings = FALSE)

owned_counts_path <- "RealTGN.salmon.merged.gene_counts.tsv"
manuscript_docx_path <- "Cheng et al 081626_VT_ASH-1.docx"
manuscript_figure_pdf_path <- "Cheng et al Figures_V2.pdf"

source_genotype <- "ANAT002"
source_genotype_sex <- "male"
complete_differentiation_lots <- c("TGN1003A", "TGN1028B", "TGN1215AvX")
tgn_line_palette <- setNames(anatomic_okabe_ito[1:3], complete_differentiation_lots)

write_manuscript_table <- function(x, filename) {
  path <- file.path(manuscript_table_dir, filename)
  write_csv(x, path, na = "NA")
  path
}

save_manuscript_figure <- function(plot, filename, width, height, dpi = 300) {
  png_path <- file.path(manuscript_figure_dir, paste0(filename, ".png"))
  pdf_path <- file.path(manuscript_figure_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white", device = "png")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white", device = "pdf")
  invisible(c(png = png_path, pdf = pdf_path))
}

assert_input <- function(path) {
  if (!file.exists(path)) {
    stop("Required input is missing: ", path, call. = FALSE)
  }
  invisible(path)
}

parse_manuscript_metadata <- function(sample_ids, counts) {
  total_counts <- colSums(counts)
  detected_genes <- colSums(counts > 0)

  tibble(sample_id = sample_ids) |>
    mutate(
      sequencing_cohort = sub("^([^_]+)_.*$", "\\1", sample_id),
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
      hDRG_library_label = if_else(sample_class == "hDRG comparator", sample_label, NA_character_),
      source_hiPSC_genotype = if_else(sample_class == "TGN culture", source_genotype, NA_character_),
      source_hiPSC_sex = if_else(sample_class == "TGN culture", source_genotype_sex, NA_character_),
      total_counts = unname(total_counts[sample_id]),
      detected_genes = unname(detected_genes[sample_id])
    )
}

load_owned_gene_counts <- function(path = owned_counts_path) {
  assert_input(path)
  imported <- read_tsv(path, show_col_types = FALSE, name_repair = "minimal")
  stopifnot(all(c("gene_id", "gene_name") %in% colnames(imported)))

  sample_ids <- setdiff(colnames(imported), c("gene_id", "gene_name"))
  counts_tbl <- imported |>
    transmute(
      gene = toupper(trimws(gene_name)),
      across(all_of(sample_ids), as.numeric)
    ) |>
    filter(!is.na(gene), gene != "") |>
    group_by(gene) |>
    summarise(across(all_of(sample_ids), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  counts <- counts_tbl |>
    column_to_rownames("gene") |>
    as.matrix()
  storage.mode(counts) <- "numeric"
  stopifnot(!anyDuplicated(rownames(counts)))

  list(
    counts = counts,
    metadata = parse_manuscript_metadata(colnames(counts), counts)
  )
}

tmm_logcpm <- function(counts) {
  dge <- DGEList(counts = counts) |>
    calcNormFactors(method = "TMM")
  cpm(dge, log = TRUE, prior.count = 2)
}

blocked_lot_trend <- function(data, value_col) {
  model_data <- data |>
    transmute(
      differentiation_lot,
      rna_week_after_thaw,
      value = .data[[value_col]]
    ) |>
    filter(
      differentiation_lot %in% complete_differentiation_lots,
      rna_week_after_thaw %in% 1:4,
      is.finite(value)
    )

  if (nrow(model_data) != 12L || n_distinct(model_data$differentiation_lot) != 3L) {
    stop("Blocked trend requires the 12 matched libraries from three differentiation lots.", call. = FALSE)
  }

  fit <- lm(value ~ rna_week_after_thaw + differentiation_lot, data = model_data)
  coefficient_table <- summary(fit)$coefficients
  slope <- coefficient_table["rna_week_after_thaw", , drop = FALSE]
  residual_df <- df.residual(fit)
  estimate <- unname(slope[1, "Estimate"])
  standard_error <- unname(slope[1, "Std. Error"])

  tibble(
    blocked_slope_per_week = estimate,
    blocked_slope_ci_low = estimate - qt(0.975, residual_df) * standard_error,
    blocked_slope_ci_high = estimate + qt(0.975, residual_df) * standard_error,
    blocked_slope_p_value = unname(slope[1, "Pr(>|t|)"]),
    n_differentiation_lots = n_distinct(model_data$differentiation_lot),
    n_libraries = nrow(model_data)
  )
}

md5_or_na <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  unname(tools::md5sum(path))
}

manuscript_theme <- function(base_size = 10) {
  theme_anatomic(base_size = base_size) +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(color = "#55525D", hjust = 0, size = rel(0.78)),
      panel.spacing = grid::unit(0.7, "lines")
    )
}
