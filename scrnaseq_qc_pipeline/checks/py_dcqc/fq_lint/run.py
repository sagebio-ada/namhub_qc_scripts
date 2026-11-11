#!/usr/bin/env python3
"""Check that each FASTQ file is structurally valid, via `fq lint`
(stjude-rust-labs/fq).

Takes 1 file (SingleTarget) or 2 (a PairedTarget's mates). Each file is
linted on its own -- a separate single-file `fq lint` call -- so the paired
read-name check (fq's P001) is deliberately never done.

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: fq lint found no structural issues in any file.
  1: fq lint found real structural issues in at least one file -- a
     genuine tool-native verdict.
  2: the check itself could not run on at least one file (the fq binary is
     missing or can't execute, or a file is missing/unreadable) -- not a
     verdict about the data. Takes precedence over 1.

`fq lint` itself exits 1 both for a structurally bad file and for an
operational error (e.g. a missing file), so this script checks that each
file exists and is readable before calling fq (missing/unreadable -> 2),
and treats fq failing before it logs its own startup as 2 too. Once fq has
started on a readable file, any nonzero exit is a verdict about the
file (1) -- validator ERROR lines, or a file fq can't parse at all (e.g. a
truncated gzip).

Metrics (one JSON object on stdout, written on exit 0 and 1):
{"files": [{"fname", "status", "records", "n_errors", "n_warnings",
"errors"}, ...]} -- "status" is "pass"/"fail" per file, "errors" holds at
most the first 10 ERROR lines. The failure/error reason goes to stderr only
for exit codes 1 and 2.

Invoked with fastq filenames relative to the current directory (nf-dcqc
stages the files there as symlinks) -- do not pass a directory.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from fq_lint_qc import FQ_INSTALL_HINT, fq_available  # noqa: E402  (also puts qc_common/ on sys.path)
from qc_common.fq_lint import parse_fq_lint_output, run_fq_lint  # noqa: E402
from scrnaseq_common import DCQC_ERROR, DCQC_FAIL, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402


_MAX_ERRORS_SHOWN = 10


def lint_one(fastq_path: Path) -> Tuple[int, Optional[str], Optional[Dict[str, object]]]:
    """Lint one file. Returns (exit code, reason or None, metrics or None)."""
    name = fastq_path.name
    if not fastq_path.is_file() or not os.access(fastq_path, os.R_OK):
        return DCQC_ERROR, f"{name}: file not found or not readable", None

    text, returncode = run_fq_lint([fastq_path])
    parsed = parse_fq_lint_output(text)
    errors = parsed["errors"]
    tail = "; ".join(ln.strip() for ln in text.splitlines()[-3:] if ln.strip())

    if returncode != 0 and "fq-lint start" not in text:
        # fq never got as far as logging its own startup (e.g. the binary
        # can't execute) -- the check didn't run, so this isn't a verdict.
        return DCQC_ERROR, f"{name}: fq lint did not start (exit {returncode}): {tail}", None

    metrics = {
        "fname": name,
        "status": "pass" if returncode == 0 else "fail",
        "records": parsed["records"],
        "n_errors": len(errors),
        "n_warnings": len(parsed["warnings"]),
        "errors": errors[:_MAX_ERRORS_SHOWN],
    }
    if returncode == 0:
        return DCQC_PASS, None, metrics
    if not errors:
        # fq started on a file that exists and is readable (checked above)
        # but couldn't parse it at all -- e.g. a truncated or corrupt gzip.
        # That's a bad file, so it's a FAIL, not an operational ERROR.
        return DCQC_FAIL, f"{name}: fq lint could not parse the file: {tail}", metrics
    more = f" (+{len(errors) - _MAX_ERRORS_SHOWN} more)" if len(errors) > _MAX_ERRORS_SHOWN else ""
    return DCQC_FAIL, f"{name}: " + "; ".join(errors[:_MAX_ERRORS_SHOWN]) + more, metrics


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq_paths", nargs="+",
                   help="1 (SingleTarget) or 2 (PairedTarget) fastq filenames, relative to the current directory")
    args = p.parse_args()

    if len(args.fastq_paths) > 2:
        emit_dcqc_reason(f"Expected 1 or 2 fastq files (py-dcqc Single/PairedTarget), got {len(args.fastq_paths)}")
        return DCQC_ERROR

    if not fq_available():
        emit_dcqc_reason(FQ_INSTALL_HINT)
        return DCQC_ERROR

    results = [lint_one(Path(p)) for p in args.fastq_paths]
    codes = [code for code, _, _ in results]
    reasons = [reason for _, reason, _ in results if reason]

    if DCQC_ERROR in codes:
        # py-dcqc doesn't read metrics after an error, so only the reason
        # (every file's, so a FAIL on the other mate isn't lost) is written.
        emit_dcqc_reason("; ".join(reasons))
        return DCQC_ERROR

    emit_dcqc_metrics({"files": [metrics for _, _, metrics in results]})
    if DCQC_FAIL in codes:
        emit_dcqc_reason("; ".join(reasons))
        return DCQC_FAIL
    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
