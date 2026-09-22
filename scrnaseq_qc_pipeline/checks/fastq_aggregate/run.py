#!/usr/bin/env python3
"""Dockerizable check: post-loop aggregate checks over every mate file's
FastQC summary vs. the ENA registry (read count, R1/R2 parity, mean length).
Wraps fastq_qc.aggregate_fastq_checks(). Pure computation over small JSON
artifacts -- no tool dependency, no network.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from fastq_qc import aggregate_fastq_checks  # noqa: E402
from qc_common import write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--entity-id", required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--fastqc-summary", action="append", required=True, dest="fastqc_summaries",
                    help="Path to a fastqc_summary.json (repeatable, one per mate file)")
    p.add_argument("--ena-metadata-json", required=True, help="ena_metadata.json produced by fastq_resolution")
    p.add_argument("--any-capped", action="store_true",
                    help="Set if any mate file's download was byte-capped for spot-QC")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    per_file_summaries = []
    for path in args.fastqc_summaries:
        with open(path) as fh:
            per_file_summaries.append(json.load(fh))

    with open(args.ena_metadata_json) as fh:
        ena_metadata = json.load(fh)

    findings = aggregate_fastq_checks(
        args.entity_id, args.sample_name, per_file_summaries, ena_metadata.get("links", []), args.any_capped,
    )
    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
