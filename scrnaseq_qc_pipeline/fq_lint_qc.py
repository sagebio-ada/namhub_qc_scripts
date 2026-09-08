"""FASTQ structural integrity via `fq lint` (stjude-rust-labs/fq).

FastQC tolerates and analyzes malformed input rather than flagging
structural problems outright -- it'll happily report quality metrics on a
file with mismatched read/quality-line lengths or duplicate read names.
`fq lint` is a dedicated validator for exactly that class of issue, so it
complements FastQC rather than duplicating it.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import List

from qc_common import Finding, make_finding

FQ_INSTALL_HINT = (
    "fq was not found on PATH. Download a prebuilt binary for your platform from "
    "https://github.com/stjude-rust-labs/fq/releases and place it on PATH "
    "(e.g. /opt/homebrew/bin on an Apple Silicon Mac -- the x86_64-apple-darwin "
    "build runs fine there under Rosetta 2)."
)

_MAX_ERRORS_SHOWN = 10


def fq_available() -> bool:
    return shutil.which("fq") is not None


def lint_fastq(entity_id: str, entity_name: str, fastq_path: Path) -> List[Finding]:
    findings: List[Finding] = []

    def add(check: str, status: str, detail: str) -> None:
        findings.append(make_finding(entity_id, entity_name, check, status, detail))

    if not fq_available():
        add("fq_lint", "FAIL", FQ_INSTALL_HINT)
        return findings

    result = subprocess.run(
        ["fq", "lint", "--lint-mode", "log", str(fastq_path)],
        capture_output=True, text=True,
    )
    # fq logs INFO lines (startup, validator list, "read N records") even on
    # success -- only ERROR-level lines are actual validation failures.
    error_lines = [ln for ln in (result.stdout + result.stderr).splitlines() if "ERROR" in ln]

    if result.returncode == 0:
        add("fq_lint", "PASS", f"{fastq_path.name}: no structural issues found")
    else:
        shown = "; ".join(error_lines[:_MAX_ERRORS_SHOWN]) or f"fq lint exited {result.returncode}"
        more = f" (+{len(error_lines) - _MAX_ERRORS_SHOWN} more)" if len(error_lines) > _MAX_ERRORS_SHOWN else ""
        add("fq_lint", "FAIL", f"{fastq_path.name}: {shown}{more}")

    return findings
