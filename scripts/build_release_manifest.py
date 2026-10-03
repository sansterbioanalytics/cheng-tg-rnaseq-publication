#!/usr/bin/env python3
"""Build an MD5/SHA-256 release manifest without private source S3 paths."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def digest(path):
    md5 = hashlib.md5(usedforsecurity=False)
    sha = hashlib.sha256()
    size = 0
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            md5.update(chunk)
            sha.update(chunk)
            size += len(chunk)
    return size, md5.hexdigest(), sha.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--copy-plan", type=Path, required=True)
    parser.add_argument("--fastq-checksums", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.copy_plan.read_text())["objects"]
    with args.fastq_checksums.open(newline="") as f:
        fastq = {r["filename"]: r for r in csv.DictReader(f, delimiter="\t")}
    if len(plan) != len(fastq) or len(fastq) != 34:
        raise ValueError("Expected exactly 34 FASTQ checksums")
    rows = []
    for item in plan:
        checksum = fastq[item["filename"]]
        if int(checksum["bytes"]) != item["size_bytes"]:
            raise ValueError(f"FASTQ size mismatch: {item['filename']}")
        rows.append([f"raw_fastq/{item['filename']}", "raw FASTQ", checksum["bytes"],
                     checksum["md5"], checksum["sha256"]])
    paths = {
        "processed/Cheng_TGN_RNAseq_salmon_gene_counts.tsv": "release_data/processed/Cheng_TGN_RNAseq_salmon_gene_counts.tsv",
        "processed/Cheng_TGN_RNAseq_sample_key.tsv": "release_data/processed/Cheng_TGN_RNAseq_sample_key.tsv",
        "processed/Cheng_TGN_RNAseq_counts_manifest.json": "release_data/processed/Cheng_TGN_RNAseq_counts_manifest.json",
        "processed/Cheng_TGN_RNAseq_hDRG_paired_controls.tsv": "release_data/source/paired_controls.salmon.merged.gene_counts.tsv",
        "processed/Cheng_TGN_RNAseq_analysis_inputs.zip": "release_data/processed/Cheng_TGN_RNAseq_analysis_inputs.zip",
        "processed/Cheng_TGN_RNAseq_analysis_input_checksums.tsv": "release_data/processed/Cheng_TGN_RNAseq_analysis_input_checksums.tsv",
        "processed/Cheng_TGN_RNAseq_DATA_README.txt": "release_data/processed/Cheng_TGN_RNAseq_DATA_README.txt",
        "code/Cheng_TGN_RNAseq_code.zip": "release_data/processed/Cheng_TGN_RNAseq_code.zip",
        "submission/Cheng_TGN_SRA_metadata.tsv": "sra_submission/Cheng_TGN_SRA_metadata.tsv",
        "submission/Cheng_TGN_SRA_metadata.xlsx": "sra_submission/Cheng_TGN_SRA_metadata.xlsx",
        "submission/Cheng_TGN_BioSample_Human1.0.tsv": "sra_submission/Cheng_TGN_BioSample_Human1.0.tsv",
        "submission/Cheng_TGN_BioSample_Human1.0.xlsx": "sra_submission/Cheng_TGN_BioSample_Human1.0.xlsx",
        "submission/Cheng_TGN_FASTQ_checksums.tsv": "sra_submission/Cheng_TGN_FASTQ_checksums.tsv",
    }
    for relative, source in paths.items():
        size, md5, sha = digest(Path(source))
        rows.append([relative, "release artifact", size, md5, sha])
    if len({r[0] for r in rows}) != len(rows):
        raise ValueError("Duplicate release filename")
    with args.out.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "role", "bytes", "md5", "sha256"])
        writer.writerows(sorted(rows))
    print(f"wrote {len(rows)} MD5/SHA-256 entries")


if __name__ == "__main__":
    main()
