#!/usr/bin/env python3
"""Check that a file is structurally valid FASTQ, via `fq lint`
(stjude-rust-labs/fq). Wraps fq_lint_qc.lint_fastq().

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: fq lint found no structural issues.
  1: fq lint found real structural issues -- a genuine tool-native verdict.
  2: the check itself could not run (the fq binary is missing, or it
     crashed) -- not a verdict about the data.

Metrics (one JSON object on stdout): {"fname": str}. The failure/error
reason goes to stderr only for exit codes 1 and 2.

Invoked with just the fastq's filename, relative to the current directory
(nf-dcqc stages the file there as a symlink) -- do not pass a directory.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from fq_lint_qc import FQ_INSTALL_HINT, fq_available, lint_fastq  # noqa: E402
from qc_common import DCQC_ERROR, DCQC_FAIL, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq_path", help="Fastq filename, relative to the current directory")
    args = p.parse_args()

    fastq_path = Path(args.fastq_path)

    if not fq_available():
        emit_dcqc_reason(FQ_INSTALL_HINT)
        return DCQC_ERROR

    # entity_id/sample_name only label the Finding dict's unused fields here
    # (the dcqc contract only inspects the one status/metric this produces).
    findings = lint_fastq(fastq_path.name, fastq_path.name, fastq_path)
    finding = findings[0]

    emit_dcqc_metrics({"fname": fastq_path.name})
    if finding["status"] == "PASS":
        return DCQC_PASS

    # fq_available() already ruled out the "binary missing" case above, so a
    # FAIL finding here unambiguously means fq lint found a real structural
    # problem, not an operational failure.
    emit_dcqc_reason(finding["metric"])
    return DCQC_FAIL


if __name__ == "__main__":
    sys.exit(main())
