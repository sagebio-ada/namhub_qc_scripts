#!/usr/bin/env python3
"""Dockerizable check: run FastQC and emit its basic-stats + per-module
findings. Wraps fastq_qc.run_fastqc_and_summarize() -- the other of the two
tool-native sources allowed to carry a real PASS/WARN/FAIL verdict (FastQC's
own published module thresholds).

Also writes the fastqc_summary.json artifact (--summary-out) that
fastq_aggregate needs to compute the post-loop ENA-registry comparisons.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from fastq_qc import run_fastqc_and_summarize  # noqa: E402
from qc_common import write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--entity-id", required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--fastq-path", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--label", default=None,
                    help="Mate filename to prefix messages with, for a multi-file (R1/R2) run; "
                         "omit for a single-file run")
    p.add_argument("--output", required=True)
    p.add_argument("--summary-out", default=None, help="Path to write fastqc_summary.json (if FastQC succeeded)")
    args = p.parse_args()

    findings, summary = run_fastqc_and_summarize(
        args.entity_id, args.sample_name, Path(args.fastq_path), Path(args.outdir), label=args.label,
    )
    write_findings_json(findings, Path(args.output))

    if summary is not None and args.summary_out:
        out_path = Path(args.summary_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as fh:
            json.dump(summary, fh, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
