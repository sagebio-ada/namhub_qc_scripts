"""FASTQ structural integrity via `fq lint` (stjude-rust-labs/fq).

FastQC tolerates and analyzes malformed input rather than flagging
structural problems outright -- it'll happily report quality metrics on a
file with mismatched read/quality-line lengths or duplicate read names.
`fq lint` is a dedicated validator for exactly that class of issue, so it
complements FastQC rather than duplicating it.

The actual `fq lint` invocation and output parsing live in
qc_common.fq_lint, shared with the nf-core-summarizing pipelines
(rnaseq_qc_pipeline/atacseq_qc_pipeline); only the PASS/FAIL judgment
below (by exit code) is specific to this pipeline.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path
from typing import List

# In a repo checkout this file lives at scrnaseq_qc_pipeline/fq_lint_qc.py,
# one level below the repo root that holds qc_common/. In a Docker image
# (see checks/*/Dockerfile) this file is flattened to /app/fq_lint_qc.py
# alongside a /app/qc_common/ copied in beside it, so the two layouts need
# different sys.path entries.
_THIS_DIR = Path(__file__).resolve().parent
if (_THIS_DIR / "qc_common").is_dir():
    sys.path.insert(0, str(_THIS_DIR))
else:
    sys.path.insert(0, str(_THIS_DIR.parent))
from qc_common.fq_lint import parse_fq_lint_output, run_fq_lint  # noqa: E402

from scrnaseq_common import Finding, make_finding

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

    text, returncode = run_fq_lint([fastq_path])
    # fq logs INFO lines (startup, validator list, "read N records") even on
    # success -- only ERROR-level lines are actual validation failures.
    error_lines = parse_fq_lint_output(text)["errors"]

    if returncode == 0:
        add("fq_lint", "PASS", f"{fastq_path.name}: no structural issues found")
    else:
        shown = "; ".join(error_lines[:_MAX_ERRORS_SHOWN]) or f"fq lint exited {returncode}"
        more = f" (+{len(error_lines) - _MAX_ERRORS_SHOWN} more)" if len(error_lines) > _MAX_ERRORS_SHOWN else ""
        add("fq_lint", "FAIL", f"{fastq_path.name}: {shown}{more}")

    return findings
