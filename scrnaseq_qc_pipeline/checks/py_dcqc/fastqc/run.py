#!/usr/bin/env python3
"""FastQC on each FASTQ file, the R1/R2 read-count comparison, and a
MultiQC report over this run's FastQC output -- all in one container.

MultiQC is here, not in its own check, because py-dcqc runs every
ExternalBaseTest in its own container with only its target files staged:
MultiQC can only summarize output written in the same container. Kraken2 is
a separate check (kraken2_qc) and is not in this MultiQC report; its numbers
are in that check's own metrics.

Target: 1 or 2 fastq filenames (a PairedTarget's two mates, or a
SingleTarget's one file), relative to the current directory.

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: FastQC ran on every file and no module reported FAIL.
  1: on at least one file, a FastQC module reported FAIL, or FastQC could
     not process the file (e.g. a truncated gzip or a corrupt FASTQ -- the
     same verdict fq_lint gives that file). The ENA-registry comparison and
     MultiQC only report information: if one of them fails, the error is
     put in the metrics and the exit code is unchanged.
  2: the check itself could not run: FastQC is not installed, a file is
     missing or unreadable, or FastQC was stopped from outside (a signal,
     the timeout, or Java running out of memory). Every file is still
     checked first, and the reason also lists any FAIL found on another.

Metrics (one JSON object on stdout, on exit 0 and 1):
  "files": per file, FastQC's total_sequences, sequence_length, pct_gc,
     each module's PASS/WARN/FAIL, and its overrepresented sequences (the
     first 20, each with count, pct and FastQC's possible source, plus how
     many there were) -- MultiQC's report has only their total percentage.
     A file FastQC couldn't process has "error" instead.
  "aggregate": the R1/R2 read-count parity, and the ENA-registry
     comparisons (optional -- skipped unless --ena-metadata-json is given,
     which py-dcqc can't supply).
  "multiqc": everything in the MultiQC report (see
     multiqc_qc.read_multiqc_data): every table, what each General
     Statistics column means, the data behind every plot, the data only
     the HTML holds (FastQC's per-base sequence content), and every section
     with the data it shows or the note it shows instead (e.g. "No samples
     found with any adapter contamination > 0.1%"). nf-dcqc keeps only the
     exit code, stdout and stderr, so the HTML report is not kept: these
     metrics are the record.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from fastq_qc import (  # noqa: E402
    FASTQC_INSTALL_HINT, FastqcInterrupted, aggregate_fastq_checks, fastqc_available, overrepresented_sequences,
    parse_fastqc_data, run_fastqc, summarize_fastqc_data,
)
from multiqc_qc import read_multiqc_data, run_multiqc  # noqa: E402
from scrnaseq_common import DCQC_ERROR, DCQC_FAIL, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402

# At most this many overrepresented sequences per file go into the metrics.
_MAX_OVERREPRESENTED = 20


def _number(value):
    """FastQC's "2,000" / "49" -> 2000 / 49; anything else (e.g. a length
    range like "35-151") is returned unchanged."""
    text = str(value).replace(",", "").strip()
    try:
        return int(text)
    except ValueError:
        try:
            return float(text)
        except ValueError:
            return value


def _short_error(exc: Exception) -> str:
    """FastQC's log is long; keep the line that says what went wrong (e.g.
    "...SequenceFormatException: Ran out of data in the middle of a fastq
    entry. Your file is probably truncated")."""
    lines = [ln.strip() for ln in str(exc).splitlines() if ln.strip()]
    for line in lines:
        if "Exception" in line and not line.startswith("at "):
            return line
    return lines[-1] if lines else type(exc).__name__


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq_paths", nargs="+",
                   help="1 (SingleTarget) or 2 (PairedTarget) fastq filenames, relative to the current directory")
    p.add_argument("--ena-metadata-json", default=None,
                   help="ena_metadata.json filename produced by fastq_resolution, for a run resolved from "
                        "SRA/ENA. Omit for a directly-uploaded fastq with no ENA registry to compare against "
                        "-- the comparison is skipped in that case.")
    p.add_argument("--any-capped", action="store_true",
                   help="Set if any file's download was byte-capped for spot-QC")
    p.add_argument("--report-name", default="multiqc_report")
    args = p.parse_args()

    if len(args.fastq_paths) > 2:
        emit_dcqc_reason(f"Expected 1 or 2 fastq files (py-dcqc Single/PairedTarget), got {len(args.fastq_paths)}")
        return DCQC_ERROR
    if not fastqc_available():
        emit_dcqc_reason(FASTQC_INSTALL_HINT)
        return DCQC_ERROR

    fastqc_dir = Path("fastqc_out")
    per_file_metrics = []
    per_file_summaries = []
    fail_reasons = []
    error_reasons = []

    # Every file is checked before deciding the exit code, so an error on
    # one file never hides a FAIL on the other.
    for fastq_path_str in args.fastq_paths:
        fastq_path = Path(fastq_path_str)
        fname = fastq_path.name

        if not fastq_path.is_file() or not os.access(fastq_path, os.R_OK):
            error_reasons.append(f"{fname}: file not found or not readable")
            continue

        try:
            data_txt = run_fastqc(fastq_path, fastqc_dir)
        except FastqcInterrupted as exc:
            error_reasons.append(f"{fname}: {exc}")
            continue
        except Exception as exc:
            # FastQC is installed and the file is readable, so FastQC giving
            # up on it means the file is bad -- the same FAIL fq_lint gives
            # a truncated or corrupt file.
            message = f"FastQC could not process the file: {_short_error(exc)}"
            fail_reasons.append(f"{fname}: {message}")
            per_file_metrics.append({"fname": fname, "error": message})
            continue

        findings, summary = summarize_fastqc_data(fname, fname, fastq_path, data_txt, label=fname)
        per_file_summaries.append(summary)
        modules = {
            f["check"].split(":", 1)[1]: f["status"]
            for f in findings if f["check"].startswith("fastqc_module:")
        }
        failed_modules = [name for name, status in modules.items() if status == "FAIL"]
        if failed_modules:
            fail_reasons.append(f"{fname}: FastQC FAIL module(s): {', '.join(failed_modules)}")

        overrepresented = overrepresented_sequences(parse_fastqc_data(data_txt)["module_lines"])
        per_file_metrics.append({
            "fname": fname,
            "total_sequences": _number(summary.get("total_sequences")),
            "sequence_length": _number(summary.get("sequence_length")),
            "pct_gc": _number(summary.get("pct_gc")),
            "modules": modules,
            "n_overrepresented_sequences": len(overrepresented),
            "overrepresented_sequences": overrepresented[:_MAX_OVERREPRESENTED],
        })

    if error_reasons:
        # py-dcqc doesn't read metrics after an error, so MultiQC is not
        # run. The reason keeps every file's FAIL too.
        emit_dcqc_reason("; ".join(error_reasons + fail_reasons))
        return DCQC_ERROR

    # The R1/R2 parity check needs only the FastQC summaries, so it always
    # runs. The two ENA-registry comparisons are optional: they need
    # --ena-metadata-json, which py-dcqc has no way to supply (it stages only
    # the target files), so under py-dcqc they're always skipped.
    ena_links = []
    ena_error = None
    if args.ena_metadata_json:
        try:
            with open(args.ena_metadata_json) as fh:
                ena_links = json.load(fh).get("links", [])
        except Exception as exc:
            ena_error = f"error: could not read {args.ena_metadata_json}: {exc}"
    aggregate_findings = aggregate_fastq_checks(
        "aggregate", "aggregate", per_file_summaries, ena_links, args.any_capped,
    )
    by_aggregate_check = {f["check"]: f["metric"] for f in aggregate_findings}
    aggregate_metrics = {
        check: by_aggregate_check[check]
        for check in ("read_count_vs_ena", "paired_fastq_parity", "mean_length_vs_registry")
        if check in by_aggregate_check
    }
    if ena_error:
        aggregate_metrics["ena_comparison"] = ena_error
    elif not args.ena_metadata_json:
        aggregate_metrics["ena_comparison"] = "skipped: no --ena-metadata-json given"

    metrics = {"files": per_file_metrics, "aggregate": aggregate_metrics}
    try:
        # nf-dcqc doesn't keep the HTML report, so everything in it goes
        # into the metrics.
        run_multiqc(fastqc_dir, Path("."), report_name=args.report_name)
        metrics["multiqc"] = read_multiqc_data(Path("."), args.report_name)
    except Exception as exc:
        stderr = getattr(exc, "stderr", None)
        metrics["multiqc"] = {"error": f"{exc}" + (f": {stderr.strip()[-500:]}" if stderr else "")}

    emit_dcqc_metrics(metrics)

    if fail_reasons:
        emit_dcqc_reason("; ".join(fail_reasons))
        return DCQC_FAIL
    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
