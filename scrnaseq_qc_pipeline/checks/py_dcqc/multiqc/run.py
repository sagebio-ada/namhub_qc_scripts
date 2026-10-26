#!/usr/bin/env python3
"""Aggregate FastQC/Kraken2 output already present in the working directory
into one MultiQC HTML report. Wraps multiqc_qc.run_multiqc(). Produces no
verdict of its own -- MultiQC "doesn't compute any new QC metrics" (see
multiqc_qc.py's docstring).

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: MultiQC ran and wrote a report. There is no verdict attached -- this
     check only ever resolves to pass or error, same as kraken2_qc and
     fastq_aggregate.
  2: MultiQC itself could not run (crash, missing binary).
(Exit code 1 is defined by the contract but genuinely unused -- there is no
condition under which this reports FAIL.)

Metrics (one JSON object on stdout): the report filename.

KNOWN DEVIATIONS from a strict py-dcqc external test:
  - It needs every FastQC/Kraken2 output file already staged in the working
    directory, not a single file or pair -- py-dcqc's Target model only
    supports SingleTarget (1 file) and PairedTarget (2); there is no real
    precedent for an arbitrary-N-file target like this one.
  - It writes its own output (the HTML report) into the working directory,
    not just the three reserved output files -- the report *is* the point
    of running this check at all, unlike a validator whose only job is a
    verdict.
Treat this as a working approximation of the contract (exit code + stdout/
stderr discipline), not a literal, submittable py-dcqc test.

Invoked with no positional file argument: it scans the current directory
(where the driver stages every upstream check's FastQC/Kraken2 output) for
input, and writes its report there too.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from multiqc_qc import run_multiqc  # noqa: E402
from scrnaseq_common import DCQC_ERROR, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--report-name", default="multiqc_report")
    args = p.parse_args()

    try:
        report_path = run_multiqc(Path("."), Path("."), report_name=args.report_name)
    except Exception as exc:
        emit_dcqc_reason(f"MultiQC failed: {exc}")
        return DCQC_ERROR

    emit_dcqc_metrics({"report_file": report_path.name})
    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
