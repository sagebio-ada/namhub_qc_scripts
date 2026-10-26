#!/usr/bin/env python3
"""Post-loop aggregate checks over every mate file's FastQC summary vs. the
ENA registry (read count, R1/R2 parity, mean length). Wraps
fastq_qc.aggregate_fastq_checks(). Pure computation over small JSON
artifacts already on disk -- no tool dependency, no network, no credentials.

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: the comparison ran. There is no verdict attached -- every one of this
     check's findings is INFO (no validated threshold), so it only ever
     resolves to pass or error, same as kraken2_qc/multiqc.
  2: the check could not run (a summary/registry file was missing or
     malformed).
(Exit code 1 is defined by the contract but genuinely unused -- there is no
condition under which this reports FAIL.)

Metrics (one JSON object on stdout): the three aggregate comparisons, where
each applies (read_count_vs_ena/paired_fastq_parity/mean_length_vs_registry).

KNOWN DEVIATION: takes a variable number of FastQC-summary input files (one
per mate), beyond py-dcqc's currently-tested Target cardinality
(SingleTarget=1, PairedTarget=2) -- same deviation already noted for
matrix_content_qc.

Invoked with the ena_metadata.json filename and one or more fastqc_summary
filenames, all relative to the current directory (nf-dcqc would stage these
as symlinks, same as any other input file) -- do not pass directories.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from fastq_qc import aggregate_fastq_checks  # noqa: E402
from qc_common import DCQC_ERROR, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ena-metadata-json", required=True, help="ena_metadata.json filename produced by fastq_resolution")
    p.add_argument("--any-capped", action="store_true",
                    help="Set if any mate file's download was byte-capped for spot-QC")
    p.add_argument("fastqc_summaries", nargs="+", help="fastqc_summary.json filename(s), one per mate file")
    args = p.parse_args()

    try:
        per_file_summaries = []
        for path in args.fastqc_summaries:
            with open(path) as fh:
                per_file_summaries.append(json.load(fh))
        with open(args.ena_metadata_json) as fh:
            ena_metadata = json.load(fh)
    except Exception as exc:
        emit_dcqc_reason(f"Could not read input: {exc}")
        return DCQC_ERROR

    findings = aggregate_fastq_checks(
        "aggregate", "aggregate", per_file_summaries, ena_metadata.get("links", []), args.any_capped,
    )
    by_check = {f["check"]: f["metric"] for f in findings}

    metrics = {}
    for check in ("read_count_vs_ena", "paired_fastq_parity", "mean_length_vs_registry"):
        if check in by_check:
            metrics[check] = by_check[check]
    emit_dcqc_metrics(metrics)
    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
