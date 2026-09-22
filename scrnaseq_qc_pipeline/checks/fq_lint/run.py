#!/usr/bin/env python3
"""Dockerizable check: FASTQ structural integrity via `fq lint`
(stjude-rust-labs/fq). Wraps fq_lint_qc.lint_fastq() -- one of the two
tool-native sources in this pipeline allowed to carry a real PASS/FAIL
verdict (fq lint's own exit code).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from fq_lint_qc import lint_fastq  # noqa: E402
from qc_common import write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--entity-id", required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--fastq-path", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    findings = lint_fastq(args.entity_id, args.sample_name, Path(args.fastq_path))
    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
