#!/usr/bin/env python3
"""Run FastQC on one file and judge it by its own per-module verdicts.
Wraps fastq_qc.run_fastqc_and_summarize() -- one of the two tool-native
sources in this pipeline allowed to carry a real verdict.

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: FastQC ran and no module reported FAIL (a WARN-only or clean result
     both count as pass here -- py-dcqc has no WARN tier, and FastQC's own
     WARN threshold isn't a hard failure).
  1: at least one FastQC module reported FAIL -- a genuine tool-native
     verdict.
  2: FastQC itself could not run (missing binary, crash) -- not a verdict
     about the data.

Metrics (one JSON object on stdout): basic stats plus every non-PASS
module's status and detail. The failure reason on exit 1 lists which
module(s) failed. The working directory stays clean -- FastQC's own output
(HTML/zip/data files) is written to a throwaway temp directory, not the
task's own working directory, and discarded after the summary is parsed.

Invoked with just the fastq's filename, relative to the current directory
(nf-dcqc stages the file there as a symlink) -- do not pass a directory.

--summary-out is an extension beyond the strict py-dcqc contract: py-dcqc's
own model has no mechanism for passing data between independent external
tests, but this pipeline's fastq_aggregate check needs each mate file's
FastQC summary to compute its own cross-file comparisons. Omit this flag
for a literal py-dcqc submission; keep it for this repo's own local driver.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from fastq_qc import FASTQC_INSTALL_HINT, fastqc_available, run_fastqc_and_summarize  # noqa: E402
from scrnaseq_common import DCQC_ERROR, DCQC_FAIL, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq_path", help="Fastq filename, relative to the current directory")
    p.add_argument("--summary-out", default=None,
                    help="[Extension beyond py-dcqc] write fastqc_summary.json here for fastq_aggregate")
    args = p.parse_args()

    fastq_path = Path(args.fastq_path)

    if not fastqc_available():
        emit_dcqc_reason(FASTQC_INSTALL_HINT)
        return DCQC_ERROR

    with tempfile.TemporaryDirectory() as tmp:
        try:
            findings, summary = run_fastqc_and_summarize(
                fastq_path.name, fastq_path.name, fastq_path, Path(tmp), label=fastq_path.name,
            )
        except Exception as exc:
            emit_dcqc_reason(f"{fastq_path.name}: FastQC failed: {exc}")
            return DCQC_ERROR

    if summary is None:
        # run_fastqc_and_summarize() already turned the underlying failure
        # into a FAIL finding (fastqc itself crashed on this file) -- that's
        # an operational failure, not a module verdict.
        reason = next((f["metric"] for f in findings if f["check"] == "fastqc"), "FastQC failed")
        emit_dcqc_reason(reason)
        return DCQC_ERROR

    modules = {
        f["check"].split(":", 1)[1]: f["status"] for f in findings if f["check"].startswith("fastqc_module:")
    }
    module_details = {
        f["check"].split(":", 1)[1]: f["metric"] for f in findings if f["check"].startswith("fastqc_module:")
    }
    metrics = {
        "fname": fastq_path.name,
        "total_sequences": summary.get("total_sequences"),
        "sequence_length": summary.get("sequence_length"),
        "pct_gc": summary.get("pct_gc"),
        "modules": modules,
        "module_details": module_details,
    }
    emit_dcqc_metrics(metrics)

    if args.summary_out:
        out_path = Path(args.summary_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as fh:
            json.dump(summary, fh, indent=2)

    failed_modules = [name for name, status in modules.items() if status == "FAIL"]
    if failed_modules:
        details = "; ".join(
            f["metric"] for f in findings if f["check"].startswith("fastqc_module:") and f["status"] == "FAIL"
        )
        emit_dcqc_reason(f"{fastq_path.name}: FAIL module(s) {', '.join(failed_modules)} -- {details}")
        return DCQC_FAIL

    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
