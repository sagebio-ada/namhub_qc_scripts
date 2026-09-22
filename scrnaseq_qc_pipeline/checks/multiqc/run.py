#!/usr/bin/env python3
"""Dockerizable step: aggregate FastQC/Kraken2 output into one MultiQC HTML
report. Wraps multiqc_qc.run_multiqc(). Produces no Findings of its own --
MultiQC "doesn't compute any new QC metrics" (see multiqc_qc.py's docstring)
-- so --output is still required for a consistent contract, but is always an
empty JSON array unless the run itself crashes.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from multiqc_qc import run_multiqc  # noqa: E402
from qc_common import make_finding, write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work-dir", required=True, help="Shared work directory to scan for FastQC/Kraken2 output")
    p.add_argument("--output-dir", required=True, help="Where to write the MultiQC HTML report")
    p.add_argument("--report-name", default="multiqc_report")
    p.add_argument("--output", required=True, help="Findings JSON path (always empty on success)")
    args = p.parse_args()

    try:
        run_multiqc(Path(args.work_dir), Path(args.output_dir), report_name=args.report_name)
        findings = []
    except Exception as exc:
        findings = [make_finding("multiqc", "multiqc", "multiqc", "FAIL", f"MultiQC failed: {exc}")]

    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
