#!/usr/bin/env python3
"""Stream an S3 FASTQ manifest and record actual content MD5 and SHA-256."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import subprocess
from pathlib import Path


def digest_object(bucket, item):
    uri = f"s3://{bucket}/{item['source_key']}"
    process = subprocess.Popen(
        ["aws", "s3", "cp", uri, "-", "--only-show-errors"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    size = 0
    while True:
        chunk = process.stdout.read(4 * 1024 * 1024)
        if not chunk:
            break
        md5.update(chunk)
        sha256.update(chunk)
        size += len(chunk)
    stderr = process.stderr.read().decode(errors="replace")
    code = process.wait()
    if code:
        raise RuntimeError(f"Download failed for {uri}: {stderr}")
    if size != item["size_bytes"]:
        raise ValueError(f"Size mismatch for {uri}: {size} != {item['size_bytes']}")
    return {"source_key": item["source_key"], "bytes": size,
            "md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    data = json.loads(args.manifest.read_text())
    output = args.output
    existing = {}
    if output.exists():
        with output.open(newline="") as f:
            existing = {r["source_key"]: r for r in csv.DictReader(f, delimiter="\t")}
    items = [i for i in data["objects"] if i["source_key"] not in existing]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(digest_object, data["bucket"], item): item for item in items}
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            result = future.result()
            existing[result["source_key"]] = result
            with output.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["source_key", "bytes", "md5", "sha256"], delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for key in sorted(existing):
                    writer.writerow(existing[key])
            print(f"{i}/{len(items)} MD5 {result['source_key']} {result['md5']}", flush=True)
    if len(existing) != len(data["objects"]):
        raise ValueError("Incomplete input digest table")


if __name__ == "__main__":
    main()
