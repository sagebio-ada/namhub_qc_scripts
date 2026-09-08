#!/usr/bin/env python3
"""QC pipeline for single-cell RNA-seq files under a Synapse project/folder.

Walks a Synapse project tree, classifies every file as a raw sequencing run
(often a link to an SRA accession rather than an uploaded FASTQ) or a
CellRanger-style processed sample (barcodes/features/matrix.mtx), then runs:

  - Raw QC (fastq_qc.py): resolve the real FASTQ behind Synapse's SRA/ENA
    link (or use an uploaded FASTQ directly), run FastQC, and check the
    downloaded data's integrity/read-count/length against ENA's own file
    report — external data checked against its own external registry.
  - GEO annotation QC (geo_metadata_qc.py): for GEO-derived raw runs, check
    that curated Synapse annotations (platform, referenceSet, library prep)
    match what GEO's own SOFT record says for that sample. This is distinct
    from curator-side schema validation (which only checks a value is
    *legal*, not that it's *correct* for this sample).
  - Processed QC (matrix_qc.py): load the CellRanger matrix.mtx.gz +
    barcodes/features, verify dimensions agree, and report per-cell
    UMI/gene/mito metrics.

This pipeline intentionally does NOT re-validate Synapse annotations against
the curation schema (required fields, controlled-vocabulary membership) --
that's already enforced by the curator tooling at entry time.

Usage
-----
    python run_scrnaseq_qc.py --synapse-id syn73675040 --output-dir output/qc

    # Also cross-check curated annotations against the GEO record:
    python run_scrnaseq_qc.py --synapse-id syn73675040 --geo-accession GSE293390

    # Fast pass, no FastQC / no matrix download:
    python run_scrnaseq_qc.py --synapse-id syn73675040 --skip-fastqc --skip-matrix

    # Spot-check raw reads without downloading full (multi-GB) fastq files:
    python run_scrnaseq_qc.py --synapse-id syn73675040 --max-download-mb 25

Requires: synapseclient, numpy, scipy, requests, and (for raw FASTQ QC)
the `fastqc` binary on PATH and the geo-synapse package:
    pip install git+https://github.com/sagebio-ada/geo_dataset_creation.git
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

_HERE = Path(__file__).parent
sys.path.insert(0, str(_HERE))

import fastq_qc
import geo_metadata_qc
import matrix_qc
from qc_common import Finding, first_value, make_finding, summarize, write_report_csv

try:
    import synapseclient as _synapse
except ImportError:
    print("synapseclient is required: pip install synapseclient", file=sys.stderr)
    raise

_SRA_RUN_URL_RE = re.compile(r"/sra/([A-Z]{2,3}\d+)/", re.IGNORECASE)


def _synapse_login():
    """Log into Synapse using SYNAPSE_AUTH_TOKEN, falling back to a cached login."""
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    syn = _synapse.Synapse()
    syn.login(authToken=token, silent=True)
    return syn


def walk_file_entities(syn, root_id: str) -> List[dict]:
    """Recursively list every FileEntity under a Synapse project/folder."""
    out: List[dict] = []
    stack = [root_id]
    while stack:
        parent = stack.pop()
        for child in syn.getChildren(parent, includeTypes=["folder", "file"]):
            if child["type"].endswith("Folder"):
                stack.append(child["id"])
            elif child["type"].endswith("FileEntity"):
                out.append(child)
    return out


def classify_entity(ent) -> str:
    """Return one of: raw_sra, raw_fastq_direct, processed_barcodes,
    processed_features, processed_matrix, other."""
    name = ent.name
    file_handle = getattr(ent, "_file_handle", None) or {}
    ext_url = file_handle.get("externalURL") or ""
    if _SRA_RUN_URL_RE.search(ext_url):
        return "raw_sra"
    if re.search(r"barcodes\.tsv", name, re.IGNORECASE):
        return "processed_barcodes"
    if re.search(r"(features|genes)\.tsv", name, re.IGNORECASE):
        return "processed_features"
    if re.search(r"matrix\.mtx", name, re.IGNORECASE):
        return "processed_matrix"
    if name.lower().endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")):
        return "raw_fastq_direct"
    return "other"


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--synapse-id", required=True, help="Synapse project or folder ID to scan (e.g. syn73675040)")
    p.add_argument("--geo-accession", default=None,
                   help="GEO series accession (e.g. GSE293390) to cross-check curated annotations against. "
                        "If omitted, that check is skipped.")
    p.add_argument("--output-dir", "-d", default="output", help="Directory to write the QC report CSV to")
    p.add_argument("--work-dir", default=None,
                   help="Directory for downloaded fastq/matrix files (default: <output-dir>/scrnaseq_qc_work)")
    p.add_argument("--max-download-mb", type=float, default=None,
                   help="Cap raw fastq downloads to this many MB for a fast spot-check "
                        "(default: download the full file). Strongly recommended when "
                        "--kraken2-db is set -- classification is slow on full files.")
    p.add_argument("--kraken2-db", default=None,
                   help="Path to an extracted Kraken2 database directory. If given, raw reads "
                        "are taxonomically classified to screen for species mismatches/contamination "
                        "(see https://benlangmead.github.io/aws-indexes/k2 for pre-built databases). "
                        "Skipped entirely if omitted.")
    p.add_argument("--kraken2-threads", type=int, default=4, help="Threads for Kraken2 classification")
    p.add_argument("--skip-fastqc", action="store_true", help="Skip raw-read FastQC QC entirely")
    p.add_argument("--skip-matrix", action="store_true", help="Skip processed CellRanger matrix QC entirely")
    p.add_argument("--keep-downloads", action="store_true",
                   help="Keep downloaded fastq/matrix files instead of deleting them after QC")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    out_dir = Path(args.output_dir)
    work_dir = Path(args.work_dir) if args.work_dir else out_dir / "scrnaseq_qc_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    print("Logging into Synapse …")
    syn = _synapse_login()

    print(f"Scanning Synapse tree under {args.synapse_id} …")
    children = walk_file_entities(syn, args.synapse_id)
    print(f"  {len(children)} file(s) found.")

    findings: List[Finding] = []
    raw_runs: Dict[str, dict] = {}
    processed_groups: Dict[str, Dict[str, str]] = defaultdict(dict)
    processed_meta: Dict[str, dict] = {}

    for child in children:
        entity_id = child["id"]
        ent = syn.get(entity_id, downloadFile=False)
        annotations = dict(ent.annotations)
        role = classify_entity(ent)

        if role == "raw_sra":
            srr = _SRA_RUN_URL_RE.search(ent._file_handle.get("externalURL", "")).group(1).upper()
            raw_runs[entity_id] = dict(name=ent.name, annotations=annotations, srr=srr)
        elif role == "raw_fastq_direct":
            raw_runs[entity_id] = dict(name=ent.name, annotations=annotations, srr=None)
        elif role.startswith("processed_"):
            sample_id = first_value(annotations.get("SampleId")) or ent.name
            sub_role = role.split("_", 1)[1]
            processed_groups[sample_id][sub_role] = entity_id
            processed_meta[sample_id] = annotations

    print(f"  {len(raw_runs)} raw run(s), {len(processed_groups)} processed sample group(s).")

    geo_by_gsm: Dict[str, dict] = {}
    if args.geo_accession:
        print(f"Fetching GEO metadata for {args.geo_accession} …")
        try:
            geo_by_gsm = geo_metadata_qc.fetch_geo_metadata_by_gsm(args.geo_accession)
            print(f"  {len(geo_by_gsm)} sample(s) found in GEO.")
        except Exception as exc:
            print(f"  Could not fetch GEO metadata: {exc}")
            geo_by_gsm = {}
        for entity_id, info in raw_runs.items():
            geo_row = geo_by_gsm.get(info["name"])
            if geo_row:
                findings.extend(geo_metadata_qc.check_geo_annotations(
                    entity_id, info["name"], info["annotations"], geo_row))
            else:
                findings.append(make_finding(entity_id, info["name"], "geo_lookup", "WARN",
                                              f"{info['name']} not found in GEO series {args.geo_accession}"))
    else:
        print("No --geo-accession given; skipping Synapse-annotations-vs-GEO check.")

    if args.skip_fastqc:
        print("Skipping raw FastQC step (--skip-fastqc).")
    else:
        for entity_id, info in raw_runs.items():
            print(f"Running FastQC QC for {info['name']} …")
            geo_row = geo_by_gsm.get(info["name"], {})
            expected_species = str(geo_row.get("sample_organism_ch1") or "") or None
            if info["srr"]:
                findings.extend(fastq_qc.qc_raw_run(
                    entity_id, info["name"], info["srr"],
                    work_dir=work_dir, max_download_mb=args.max_download_mb,
                    keep_files=args.keep_downloads,
                    kraken2_db=args.kraken2_db, expected_species=expected_species,
                    kraken2_threads=args.kraken2_threads,
                ))
            else:
                findings.extend(fastq_qc.qc_direct_fastq(
                    syn, entity_id, info["name"],
                    work_dir=work_dir, keep_files=args.keep_downloads,
                    kraken2_db=args.kraken2_db, expected_species=expected_species,
                    kraken2_threads=args.kraken2_threads,
                ))

    if args.skip_matrix:
        print("Skipping processed matrix QC step (--skip-matrix).")
    else:
        for sample_id, files in processed_groups.items():
            print(f"Running matrix QC for sample {sample_id} …")
            findings.extend(matrix_qc.qc_processed_sample(
                syn, sample_id, files, work_dir,
                cellranger_output_class=first_value(processed_meta[sample_id].get("CellrangerOutputClass")),
            ))

    report_path = out_dir / "scrnaseq_qc_report.csv"
    write_report_csv(findings, report_path)

    print(f"\nQC complete: {summarize(findings)}")
    print(f"Report written to {report_path}")

    return 1 if any(f["status"] == "FAIL" for f in findings) else 0


if __name__ == "__main__":
    sys.exit(main())
