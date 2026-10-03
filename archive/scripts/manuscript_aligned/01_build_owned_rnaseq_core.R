#!/usr/bin/env Rscript

# Rebuild only the RNA-seq analyses supported by the in-house merged count table.
# This script does not reconstruct the human TG reference or deconvolution model.

source("scripts/manuscript_aligned/00_shared.R")

owned <- load_owned_gene_counts()
counts <- owned$counts
metadata <- owned$metadata

matched_metadata <- metadata |>
  filter(
    sample_class == "TGN culture",
    differentiation_lot %in% complete_differentiation_lots,
    rna_week_after_thaw %in% 1:4
  ) |>
  arrange(differentiation_lot, rna_week_after_thaw)

stopifnot(
  nrow(metadata) == 19L,
  nrow(matched_metadata) == 12L,
  n_distinct(matched_metadata$differentiation_lot) == 3L,
  all(table(matched_metadata$differentiation_lot) == 4L),
  all(unique(matched_metadata$source_hiPSC_genotype) == source_genotype)
)

manuscript_scope <- tribble(
  ~item, ~manuscript_aligned_value, ~interpretation_boundary,
  "Source hiPSC genotype", source_genotype, "Single male genotype; no cross-genotype generalization",
  "RNA-seq biological unit", "Three independent differentiation lots", "One bulk library per lot per week",
  "RNA-seq longitudinal window", "Weeks 1-4 after thawing/plating", "Not weeks after induction; Day-7 neurons were cryopreserved before this window",
  "RNA-seq longitudinal libraries", "12", "Three lots multiplied by four weekly timepoints",
  "hDRG comparator", "Four libraries from two commercial pooled-donor RNA lots", "Exact library-to-lot mapping should be confirmed before final Methods",
  "Primary paper claim", "NaV1.8 functional platform", "Broad sensory subtype labels are contextual and supplemental",
  "Primary RNA-seq role", "Molecular and developmental context", "RNA is not a proxy for functional channel competence"
)
write_manuscript_table(manuscript_scope, "manuscript_analysis_scope.csv")
write_manuscript_table(metadata, "owned_sample_metadata.csv")

main_gene_manifest <- tribble(
  ~gene, ~figure_block, ~manuscript_role,
  "SIX1", "Cranial sensory lineage", "trigeminal-lineage-associated developmental context",
  "EYA1", "Cranial sensory lineage", "SIX1 co-regulatory cranial sensory competence",
  "SOX10", "Cranial sensory lineage", "neural-crest-associated state and purity boundary",
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
  "TRPM8", "Sensory context and boundaries", "cold-sensory context retained in the main heatmap",
  "P2RX3", "Sensory context and boundaries", "sensory/nociceptor-associated context",
  "NTRK1", "Sensory context and boundaries", "TRKA nociceptor-lineage context",
  "NTRK2", "Sensory context and boundaries", "TRKB sensory context",
  "NTRK3", "Sensory context and boundaries", "TRKC sensory context",
  "PIEZO2", "Sensory context and boundaries", "mechanosensory context",
  "CALCA", "Sensory context and boundaries", "peptidergic boundary",
  "CALCB", "Sensory context and boundaries", "peptidergic boundary"
) |>
  mutate(present_in_owned_counts = gene %in% rownames(counts))

if (!all(main_gene_manifest$present_in_owned_counts)) {
  stop(
    "Required manuscript genes missing from owned count table: ",
    paste(main_gene_manifest$gene[!main_gene_manifest$present_in_owned_counts], collapse = ", "),
    call. = FALSE
  )
}
write_manuscript_table(main_gene_manifest, "main_figure_gene_manifest.csv")

context_metadata <- metadata |>
  filter(
    sample_class == "hDRG comparator" |
      differentiation_lot %in% complete_differentiation_lots
  )
context_counts <- counts[, context_metadata$sample_id, drop = FALSE]
context_logcpm <- tmm_logcpm(context_counts)

context_expression <- as.data.frame(context_logcpm[main_gene_manifest$gene, , drop = FALSE]) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "tmm_log_cpm") |>
  left_join(context_metadata, by = "sample_id") |>
  left_join(main_gene_manifest, by = "gene")
write_manuscript_table(context_expression, "owned_main_gene_expression_logcpm.csv")

tgn_week_summary <- context_expression |>
  filter(
    sample_class == "TGN culture",
    differentiation_lot %in% complete_differentiation_lots,
    rna_week_after_thaw %in% 1:4
  ) |>
  group_by(gene, figure_block, manuscript_role, rna_week_after_thaw) |>
  summarise(
    mean_tmm_log_cpm = mean(tmm_log_cpm),
    sd_tmm_log_cpm = sd(tmm_log_cpm),
    n_differentiation_lots = n_distinct(differentiation_lot),
    .groups = "drop"
  ) |>
  mutate(context = paste0("TGN W", rna_week_after_thaw))

hdrg_summary <- context_expression |>
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

main_heatmap_summary <- bind_rows(
  tgn_week_summary |> mutate(n_hDRG_libraries = NA_integer_),
  hdrg_summary
) |>
  select(
    gene, figure_block, manuscript_role, context,
    mean_tmm_log_cpm, sd_tmm_log_cpm,
    n_differentiation_lots, n_hDRG_libraries
  )
write_manuscript_table(main_heatmap_summary, "owned_main_heatmap_summary.csv")

matched_counts <- counts[, matched_metadata$sample_id, drop = FALSE]
matched_logcpm <- tmm_logcpm(matched_counts)
matched_expression <- as.data.frame(matched_logcpm[main_gene_manifest$gene, , drop = FALSE]) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "tmm_log_cpm") |>
  left_join(matched_metadata, by = "sample_id") |>
  left_join(main_gene_manifest, by = "gene")

direct_trends <- matched_expression |>
  group_by(gene, figure_block, manuscript_role) |>
  group_modify(~ blocked_lot_trend(.x, "tmm_log_cpm")) |>
  ungroup()
write_manuscript_table(direct_trends, "owned_direct_gene_blocked_trends.csv")

purity_evidence_map <- tribble(
  ~evidence_domain, ~current_manuscript_evidence, ~denominator_status, ~safe_language, ~needed_source_material,
  "Manufacturing QC", "Pass criterion states neuronal purity greater than 99%", "Lot-level observed values not present in this repository", "Manufacturing lots met predefined QC criteria", "Per-lot purity, yield, viability, and post-thaw recovery table",
  "Immunostaining", "Representative neuronal and sensory-marker images", "Cells, fields, and lots not quantified", "Representative marker expression", "Blinded counts by field and lot with prespecified denominator",
  "Transcriptomics", "Neuronal and sensory programs with explicit non-neuronal counter-programs", "Bulk RNA-seq cannot provide cell counts", "Transcriptionally enriched for neuronal programs", "Keep similarity language; do not convert scores to purity percentages"
)
write_manuscript_table(purity_evidence_map, "purity_evidence_integration_map.csv")

owned_provenance <- tribble(
  ~artifact, ~source_path, ~analysis_action, ~ownership_boundary,
  "Merged gene counts", owned_counts_path, "Imported and re-normalized with edgeR TMM", "Reproducible from the in-house count table",
  "Direct gene heatmap", owned_counts_path, "Recomputed from matched TGN lots and hDRG libraries", "Owned downstream summary",
  "Blocked gene trends", owned_counts_path, "Recomputed with week plus differentiation-lot model", "Owned downstream inference",
  "Functional milestones", manuscript_docx_path, "Transcribed as manuscript claims only", "Not reanalyzed; electrophysiology source data unavailable"
) |>
  mutate(source_md5 = vapply(source_path, md5_or_na, character(1)))
write_manuscript_table(owned_provenance, "owned_analysis_provenance.csv")

message("Owned RNA-seq core written to: ", normalizePath(manuscript_output_dir))
