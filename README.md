# Cheng et al. TG RNA-seq analysis code

This repository contains analysis code for the RNA-seq panels in Cheng et al., “hiPSC-derived trigeminal sensory neurons” (final manuscript draft dated 2026-09-28). The publication entry point is [`publication_notebook/Cheng_TG_RNAseq_publication.qmd`](publication_notebook/Cheng_TG_RNAseq_publication.qmd). Figures 2C–F and Supplemental Figure 1A–F are generated there. The notebook retains `s2*` output filenames from an earlier figure numbering scheme; in the final manuscript these are Supplemental Figure 1 panels. Supplemental Figure 2 is antibody validation and is outside this RNA-seq code.

## Reproduce

The notebook requires the corrected paired-control Salmon gene-count matrix, the earlier TGN count matrix, the GSE197289 derived reference objects and manifests, and the small sample/gene tables cited in the notebook. These inputs are not committed to this Git repository. The corrected count matrix must match the SHA-256 pinned in the notebook. See [input and run instructions](publication_notebook/README.md) and [the required data manifest](DATA_REQUIREMENTS.md). Once the data are available at the documented relative paths, set `TG_CORRECTED_COUNTS` to the corrected matrix and run:

```sh
quarto render publication_notebook/Cheng_TG_RNAseq_publication.qmd --to html
```

This command needs Quarto, R, and the R packages listed in the notebook README. The code stops if a required input or the pinned corrected matrix hash is missing. Generated outputs go to `publication_notebook/output/` and are excluded from Git.

## Release scope

Earlier notebooks are in [`archive/notebooks/`](archive/notebooks/) as exploratory history. They do not generate the final publication panels. The final figure code is in `publication_notebook/`. Raw reads, derived analysis objects, and final panel files require a separate public data deposit before the analysis is independently reproducible. The manuscript's data and code availability statement still contains repository placeholders; it must be updated with the final accession and release DOI.

The [SRA and BioSample submission tables](sra_submission/README.md) share the publication sample IDs used in the processed count matrix. These tables are drafts for depositor review and have not been submitted to NCBI.
