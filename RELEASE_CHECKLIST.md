# Public release gate

- [x] Identify the final RNA-seq code entry point: `publication_notebook/Cheng_TG_RNAseq_publication.qmd`.
- [x] Map code outputs to final manuscript: Figure 2C–F and Supplemental Figure 1A–E. Historical `s2*` filenames are documented.
- [x] Exclude manuscript drafts, slide decks, working runs, credentials, and rendered outputs from Git.
- [x] Stage the corrected paired-control matrix, reference objects, and public 15-sample matrix in the private S3 release candidate.
- [ ] Confirm public redistribution rights and deposit those inputs with stable public identifiers.
- [x] Render the notebook end to end with the validated corrected matrix; internal Figure 2C values match 322/322 comparison rows.
- [ ] Compare rendered panel files with the exact final manuscript assembly.
- [ ] Review the final manuscript's figure legend and data/code availability placeholders with the authors.
- [ ] Review the exact Git diff and repository visibility before making the GitHub repository public.
- [ ] Create a version tag and GitHub release; then archive that release in Zenodo and record its DOI in the repository and manuscript.

A GitHub or Zenodo release should wait until the data inputs and author review are complete.

- [x] Confirm hDRG comparators are two separate supplier RNA pools, each sequenced twice; author reports no human-material privacy concern.
- [ ] Review remaining BioSample Human 1.0 missing-value fields before NCBI submission.
