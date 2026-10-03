# NCBI submission review

The SRA and Human 1.0 BioSample tables each contain 15 biological samples, with sample names identical to the public processed count matrix. Each of the 34 FASTQ filenames appears exactly once in the SRA table. C1/C2 are technical sequencing runs of each hDRG comparator and are listed together in one SRA row per comparator.

The Human 1.0 table uses `not collected` for mandatory source attributes that are not available in the draft metadata. Lot and week columns distinguish the TGN cultures. The author confirmed that hDRG comparator 2 is a separate pooled RNA source from the same supplier, sequenced in C1 and C2. Donor age/sex composition and collection details for comparator 2 remain unreported and are marked with controlled missing-value terms. Its source ID `TAK2204464A` is preserved only in the private source mapping; the public sample name is `Cheng_hDRG2`. BioProject/SAMN accessions are not yet assigned.

The SRA TSV is the submission table; the XLSX is a review copy of the supplied SRA template. The Human 1.0 TSV is the BioSample attributes submission table; its XLSX is a review copy. Use the 34 publication-named FASTQs shown in the checksum manifest for upload. These files have not been submitted to NCBI.
