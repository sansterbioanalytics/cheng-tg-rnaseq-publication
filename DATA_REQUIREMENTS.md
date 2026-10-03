# Data required for a full render

The 15-sample Salmon estimated gene-count matrix and matching sample metadata are planned for NCBI GEO under BioProject PRJNA1538692. The GEO accession is pending. That public matrix is not the same file as the corrected paired-control input pinned by the publication notebook. A full render also needs the frozen reference objects and other inputs listed below; their public locations are pending.

Extract `Cheng_TGN_RNAseq_analysis_inputs.zip` under the repository root to place the following files at these relative paths. They are deliberately excluded from Git; their public deposit locations are still pending. The separate `Cheng_TGN_RNAseq_hDRG_paired_controls.tsv` is specified by `TG_CORRECTED_COUNTS` and must match the SHA-256 in the notebook.

| Relative path | Local size |
| --- | ---: |
| `RealTGN.salmon.merged.gene_counts.tsv` | 5,589,051 bytes |
| `results/consensus_identity/GSE197289_human_TG_neurons_consensus.rds` | 123,198,948 bytes |
| `results/pseudobulk/GSE197289_TG_neuron_donor_subtype_logCPM.csv` | 17,223,757 bytes |
| `results/deconv_rerun_expanded/ensembl_to_symbol_id_map.csv` | 3,700,452 bytes |
| `results/pseudobulk/GSE197289_TG_combined_reference_edgeR_dge.rds` | 2,065,409 bytes |
| `results/pseudobulk/GSE197289_TG_combined_reference_metadata.csv` | 5,551 bytes |
| `results/marker_refinement/GSE197289_TG_expanded_marker_manifest.csv` | 133,587 bytes |
| `results/pseudobulk/GSE197289_TG_combined_marker_manifest.csv` | 2,630 bytes |
| `results/pseudobulk/GSE197289_TG_neuron_donor_subtype_edgeR_dge.rds` | 1,452,612 bytes |
| `results/pseudobulk/GSE197289_TG_neuron_donor_subtype_metadata.csv` | 3,111 bytes |
| `results/pseudobulk/GSE197289_TG_combined_reference_logCPM.csv` | 24,005,896 bytes |

The five small CSV files in `publication_notebook/inputs/` are committed with the notebook. The 123 MB Seurat reference object exceeds GitHub’s 100 MB per-file limit, so it belongs in a data repository or release asset with a stable identifier.

The corrected paired-control count matrix was downloaded from the prior analysis run and verified against the notebook SHA-256. The notebook rendered end to end locally. The generated panels still require comparison with the exact final manuscript assembly before publication.
