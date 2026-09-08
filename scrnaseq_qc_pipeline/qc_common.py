"""Shared helpers for the scRNA-seq QC pipeline (findings format, report writer)."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Union

Finding = Dict[str, str]

REPORT_FIELDS = ["sample", "entity_id", "check", "status", "detail"]
_SEVERITY_RANK = {"FAIL": 0, "WARN": 1, "PASS": 2}


def first_value(raw: Optional[Union[list, tuple, str]]) -> str:
    """Synapse annotations are always lists; take the first value or ''."""
    if isinstance(raw, (list, tuple)):
        return str(raw[0]) if raw else ""
    return str(raw) if raw else ""


def make_finding(entity_id: str, sample: str, check: str, status: str, detail: str) -> Finding:
    assert status in _SEVERITY_RANK, f"unknown status {status!r}"
    return {"entity_id": entity_id, "sample": sample, "check": check, "status": status, "detail": detail}


def write_report_csv(findings: Iterable[Finding], path: Path) -> None:
    rows = sorted(findings, key=lambda f: (_SEVERITY_RANK.get(f["status"], 1), f["sample"], f["check"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=REPORT_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in REPORT_FIELDS})


def summarize(findings: List[Finding]) -> str:
    n_fail = sum(1 for f in findings if f["status"] == "FAIL")
    n_warn = sum(1 for f in findings if f["status"] == "WARN")
    n_pass = sum(1 for f in findings if f["status"] == "PASS")
    return f"{n_pass} PASS, {n_warn} WARN, {n_fail} FAIL"
