# Cheng TG RNA-seq publication code export

The [Quarto analysis notebook](Cheng_TG_RNAseq_publication.qmd) is the curated publication-facing code path for Figure 2C–F and Supplemental Figure 1A–F and detachable reference/sample keys. It uses the corrected paired-control matrix and frozen reference sources, with PDF/SVG/600-dpi PNG files for each export. The final manuscript labels these supplemental RNA-seq panels Supplemental Figure 1A–E; the older `s2*` filenames in the code are historical export names. A complete render writes its panels and numerical tables to `output/`. The manuscript figure assembly includes non-RNA-seq images supplied by the authors.

Author-confirmed identity mapping: source `hDRG` = hDRG 1; source `TAK2204464A` = hDRG 2. The code and tables retain the source identifiers and provide `sample_label_key.csv` for plot labels. C1/C2 are technical sequencing runs within each of these two source samples, not four biological replicates.

## Reproduce this notebook

1. Download the corrected `salmon.merged.gene_counts.tsv` from the study data deposit (accession pending) to a local path. Its expected SHA-256 is pinned in the notebook.
2. From the repository root, set `TG_CORRECTED_COUNTS` to the local absolute path and run `quarto render publication_notebook/Cheng_TG_RNAseq_publication.qmd --to html`.
3. Review `publication_notebook/output/` for input/panel manifests, normalization and label tables, underlying value/coordinate tables, comparisons with legacy results, the individual panels, and both detachable keys. Generated results are deliberately excluded from Git by the existing ignore rules. Share data/results separately from code.

Dependencies: Quarto, R, `edgeR`, `ggplot2`, `digest`, `knitr`, `Seurat`, `ggrepel`, `nnls`, `Rtsne`, and `scales`. The notebook does not install software or fetch protected inputs during rendering. It fails if the corrected matrix hash, gene IDs, gene names, marker-set sizes, or biological sample IDs differ from the pinned run. In a restricted environment, set `XDG_CACHE_HOME` and `DENO_DIR` to writable temporary directories before invoking Quarto.

## Code and result boundary

| Role | Current source | Export decision |
| --- | --- | --- |
| Figure 2C corrected bulk heatmap | `publication_notebook/Cheng_TG_RNAseq_publication.qmd` | Generated for two hDRG reference samples; all 322 plotted values match the paired-control export exactly. Author approves final layout and legend. |
| Figure 2D UMAP, Supplemental Figure 1A marker UMAP, and Supplemental Figure 1B reference expression | `publication_notebook/Cheng_TG_RNAseq_publication.qmd` | Regenerated from GSE197289-only sources. Supplemental Figure 1B has genes in columns and hierarchically clustered identities in rows; values are mean log2 CPM and color is centered within each gene. Supplemental Figure 1A and Supplemental Figure 1B are both 6.8 × 7.3 inches. |
| Figure 2E t-SNE and detachable shared key | `publication_notebook/Cheng_TG_RNAseq_publication.qmd` | Seven-PC t-SNE includes explicit hDRG 1/2 and W1–W4 labels; `figure2e_nearby_reference_label_key.csv` records the direct-labeled source pseudobulks. Reference colors and short labels match the UMAP. |
| Figure 2F refined signatures, Supplemental Figure 1C conditioning, and Supplemental Figure 1E original/refined comparison | `publication_notebook/Cheng_TG_RNAseq_publication.qmd` | Figure 2F is the refined row of Supplemental Figure 1E: the same 0–100% stacked bars and colors, with the three lots as columns. Supplemental Figure 1C retains per-identity 79/788-gene VIF bars; Supplemental Figure 1E adds the original fit above the refined fit. These are relative fitted weights, not cell fractions. |
| Supplemental Figure 1D expression trajectories | `publication_notebook/Cheng_TG_RNAseq_publication.qmd` | All seven gene-set scores are displayed, including negative satellite-glia and Schwann-cell scores. The Aβ marker-contrast row is removed from the figure. All ten calculated axes remain in the 120-row values table. Absolute program scores differ materially from the legacy table; review Results/legend. |
| Supplemental Figure 1F refined cluster marker expression | `publication_notebook/Cheng_TG_RNAseq_publication.qmd` | SCT assay `data` values for eight paragraph markers across 18 donor-supported clusters; [legend and interpretation](Supplemental_Figure_1F_marker_expression_legend.md). |
| Marker refinement | `5_TG_single_cell_marker_refinement.qmd` | Preserve as reference-method provenance; freeze the chosen marker manifest before a joint rerun. |

The signed-refined 788-gene rule and seven-PC t-SNE method were recovered from the archived rendered analysis record and are now implemented in this notebook with frozen input checksums. Older exploratory notebooks remain available as methodological history, not the generation entry point. The final manuscript uses Figure 2C–F and Supplemental Figure 1A–F for these outputs. See the final manuscript for the approved figure legends. Before submission, reconcile the new Supplemental Figure 1D score scale and all affected captions/Results, inspect the final page assembly with Vince/Patrick's IHC assets, and confirm the hDRG source-label mapping in submission metadata. Keep source tables/results in a separate distribution; code export should remain lean and independently inspectable.

## Signature comparison

`s2e_original_vs_refined_weights.{pdf,svg,png}` compares the original 79-gene and refined 788-gene fits for the same 12 matched lot-week samples. Its stacked bars use a common percentage scale. The VIF comparison in Supplemental Figure 1C is the separate evidence for improved numerical conditioning; fitted-weight changes alone do not show biological accuracy or subtype conversion.
