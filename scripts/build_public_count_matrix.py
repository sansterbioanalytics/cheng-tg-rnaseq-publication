#!/usr/bin/env python3
"""Assemble a publication Salmon estimated gene-count matrix from two validated runs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from public_sample_names import public_name

CORRECTED_SHA256 = "8e2a974c81b17c8a4988dfda79ec27a74d7ce8bb060219add3ec53688a9b24ad"
GENE_ROWS = 62757


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def numeric(value, label):
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"Invalid count for {label}: {value!r}") from exc
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"Invalid count for {label}: {value!r}")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy", type=Path, required=True)
    parser.add_argument("--corrected", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if sha256(args.corrected) != CORRECTED_SHA256:
        raise ValueError("Corrected matrix SHA-256 does not match the validated run")
    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = args.outdir / "Cheng_TGN_RNAseq_salmon_gene_counts.tsv"
    mapping = args.outdir / "Cheng_TGN_RNAseq_sample_key.tsv"
    with args.legacy.open(newline="") as old_file, args.corrected.open(newline="") as new_file, matrix.open("w", newline="") as out_file:
        old = csv.DictReader(old_file, delimiter="\t")
        new = csv.DictReader(new_file, delimiter="\t")
        tgn = [c for c in old.fieldnames if c.startswith(("C1_TGN", "C2_TGN"))]
        if len(tgn) != 13 or len({c[3:] for c in tgn}) != 13:
            raise ValueError(f"Expected 13 unique TGN libraries; found {tgn}")
        if new.fieldnames != ["gene_id", "gene_name", "hDRG", "TAK2204464A"]:
            raise ValueError(f"Unexpected corrected matrix columns: {new.fieldnames}")
        output_names = [public_name(c) for c in tgn] + ["Cheng_hDRG1", "Cheng_hDRG2"]
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "gene_name", *output_names])
        count = 0
        for old_row, new_row in zip(old, new, strict=True):
            if (old_row["gene_id"], old_row["gene_name"]) != (new_row["gene_id"], new_row["gene_name"]):
                raise ValueError(f"Gene identity mismatch at row {count + 1}")
            values = [numeric(old_row[c], c) for c in tgn]
            values.extend([numeric(new_row[c], c) for c in ("hDRG", "TAK2204464A")])
            writer.writerow([old_row["gene_id"], old_row["gene_name"], *values])
            count += 1
        if count != GENE_ROWS:
            raise ValueError(f"Expected {GENE_ROWS} genes, found {count}")
    with mapping.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["output_sample", "source_matrix", "source_column", "role"])
        for col in tgn:
            writer.writerow([public_name(col), "legacy", col, "TGN culture" if col != "C1_TGN251028BW1" else "pilot TGN culture"])
        writer.writerow(["Cheng_hDRG1", "corrected", "hDRG", "hDRG comparator"])
        writer.writerow(["Cheng_hDRG2", "corrected", "TAK2204464A", "hDRG comparator"])
    manifest = {
        "matrix": matrix.name,
        "matrix_sha256": sha256(matrix),
        "gene_rows": count,
        "biological_samples": len(output_names),
        "source_map": mapping.name,
        "source_map_sha256": sha256(mapping),
        "legacy_matrix_sha256": sha256(args.legacy),
        "corrected_matrix_sha256": sha256(args.corrected),
        "count_type": "Salmon estimated gene counts; non-integer estimates, not raw read counts",
    }
    (args.outdir / "Cheng_TGN_RNAseq_counts_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
