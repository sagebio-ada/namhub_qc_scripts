"""Shared helpers for the scRNA-seq QC pipeline (findings format, report writer)."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Union

Finding = Dict[str, str]

# This pipeline is still at the preliminary/exploratory stage: we're using it
# to understand what these datasets actually contain, not to run an
# established, validated gate. So the report leads with `metric` (the actual
# reported number/fact) and pushes `status` later -- reading the numbers is
# the point right now, not triaging by verdict. Re-sort/reorder this once the
# pipeline moves past exploration and verdicts start carrying real weight.
REPORT_FIELDS = ["sample", "check", "metric", "status", "entity_id"]
# INFO is for findings that report a real number with no pass/fail judgment
# attached (no verified criterion, or a criterion we didn't trust enough to
# assert as a verdict) -- it's distinct from PASS, which means a check was
# actually applied and it succeeded. Don't reuse PASS for "here's some data."
_VALID_STATUSES = {"FAIL", "WARN", "PASS", "INFO"}


def first_value(raw: Optional[Union[list, tuple, str]]) -> str:
    """Synapse annotations are always lists; take the first value or ''."""
    if isinstance(raw, (list, tuple)):
        return str(raw[0]) if raw else ""
    return str(raw) if raw else ""


def make_finding(entity_id: str, sample: str, check: str, status: str, detail: str) -> Finding:
    assert status in _VALID_STATUSES, f"unknown status {status!r}"
    return {"entity_id": entity_id, "sample": sample, "check": check, "status": status, "metric": detail}


def write_report_csv(findings: Iterable[Finding], path: Path) -> None:
    # Sorted for readable grouping (by sample, then check) -- NOT by
    # severity. See the note on REPORT_FIELDS above.
    rows = sorted(findings, key=lambda f: (f["sample"], f["check"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=REPORT_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in REPORT_FIELDS})


def write_findings_json(findings: Iterable[Finding], path: Path) -> None:
    """Write Findings as a plain JSON array to an explicit path -- an extra,
    local-tooling-only artifact (never one of py-dcqc's three reserved
    output files) that lets our own driver script read back the full
    sample/check/status/metric shape for its own report, alongside the
    genuine py-dcqc stdout/stderr/exit-code contract below."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(list(findings), fh, indent=2)


# --- py-dcqc ExternalBaseTest contract ------------------------------------
#
# Per Sage-Bionetworks-Workflows/py-dcqc#96 ("Expanded architecture document
# to include Docker image instructions"): an external test's container runs
# one command and is judged ONLY by its exit code, stdout, and stderr --
# nf-dcqc captures these into exit_code.txt/std_out.txt/std_err.txt itself.
# Three exit codes, not two:
#   DCQC_PASS  (0): the file passed the check.
#   DCQC_FAIL  (1): the file failed the check -- a genuine tool-native
#                   bad-file verdict (what this pipeline's own philosophy
#                   calls a real verdict, not a judgment call).
#   DCQC_ERROR (2): the tool could not complete the check at all (a crash,
#                   a missing binary, a network failure) -- this is what
#                   this pipeline's own code has always called FAIL
#                   elsewhere (a genuine operational failure, never a
#                   data-quality opinion). Do not confuse the two: py-dcqc's
#                   FAIL and this pipeline's FAIL are not the same thing.
# Most of this pipeline's checks have no validated tool-native verdict at
# all (the whole point of the INFO status) -- for those, DCQC_FAIL is simply
# never used: the check only ever resolves to DCQC_PASS (normal) or
# DCQC_ERROR (couldn't run), with the actual numbers carried in the metrics
# JSON, not implied by the exit code.
DCQC_PASS = 0
DCQC_FAIL = 1
DCQC_ERROR = 2


def emit_dcqc_metrics(metrics: dict) -> None:
    """Write metrics as the one JSON object py-dcqc's external-test contract
    expects on stdout. Call this at most once, and never print anything
    else to stdout -- any other text makes the JSON invalid and py-dcqc
    reports the test as errored regardless of the exit code."""
    print(json.dumps(metrics))


def emit_dcqc_reason(reason: str) -> None:
    """Write a short failure/error reason to stderr -- the stream py-dcqc
    copies into status_reason for a failed or errored test. Skip calling
    this for a passing result; status_reason stays empty on pass."""
    if reason:
        print(reason, file=sys.stderr)


def summarize(findings: List[Finding]) -> str:
    n_fail = sum(1 for f in findings if f["status"] == "FAIL")
    n_warn = sum(1 for f in findings if f["status"] == "WARN")
    n_pass = sum(1 for f in findings if f["status"] == "PASS")
    n_info = sum(1 for f in findings if f["status"] == "INFO")
    return f"{n_pass} PASS, {n_warn} WARN, {n_fail} FAIL, {n_info} INFO"
