"""Run and parse `fq lint` (stjude-rust-labs/fq) output.

Shared by every pipeline that validates raw FASTQ structure:
scrnaseq_qc_pipeline's dockerized `fq_lint` check (via fq_lint_qc.py) and
the nf-core-summarizing pipelines (via qc_common.tools/fastq_sections).

`--lint-mode log` reports every structural issue found instead of stopping
at the first one (fq lint's own default, "panic", exits as soon as it hits
one). `--disable-validator P001` turns off the paired-read name-correspondence
validator; it's a no-op on a single-file call (no paired validators run
without an R2) and is otherwise desirable so a mismatched read-name
convention between mates doesn't get reported as a structural defect.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Union

ANSI = re.compile(r"\x1b\[[0-9;]*m")
_ARGS = ["--lint-mode", "log", "--disable-validator", "P001"]


def run_fq_lint(fastqs: List[Union[str, Path]], fq: str = "fq") -> Tuple[str, int]:
    """Run `fq lint` on one FASTQ (single-end) or a pair.

    Returns (combined stdout+stderr text, exit code).
    """
    result = subprocess.run([fq, "lint", *_ARGS, *map(str, fastqs)], capture_output=True, text=True)
    return result.stdout + result.stderr, result.returncode


def parse_fq_lint_output(text: str) -> Dict[str, object]:
    """fq lint's log-mode output -> structural facts about the run.

    "records" is None if fq lint never reached its final summary line (e.g.
    it crashed before finishing).
    """
    text = ANSI.sub("", text)
    lines = text.splitlines()
    errors = [ln for ln in lines if re.search(r"\bERROR\b", ln)]
    warnings = [ln for ln in lines if re.search(r"\bWARN\b", ln)]
    recs = re.findall(r"read (\d+)(?: \* \d+)? records", text)
    return {
        "errors": errors,
        "warnings": warnings,
        "records": int(recs[-1]) if recs else None,
        "completed": "fq-lint end" in text,
    }
