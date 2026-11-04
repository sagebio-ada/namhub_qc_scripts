#!/usr/bin/env python3
"""Full per-run raw-read QC: FastQC + Kraken2 on each mate file, the
cross-mate/ENA-registry aggregate comparison, and a MultiQC report over this
run's own FastQC/Kraken2 output -- all in one container.

Merges what were four separate checks (fastqc, kraken2_qc, fastq_aggregate,
multiqc) into one py-dcqc test, per andrewelamb's PR #1 review
(2026-11-03). py-dcqc runs every ExternalBaseTest in its own isolated
container with only its own target file(s) staged, so fastq_aggregate's
cross-mate comparison and multiqc's "scan the working directory for other
tools' output" could never see another test's results -- those two were
never actually runnable as separate py-dcqc tests. Folding all four into
the container that produces the files they need is the only way to keep
this functionality inside py-dcqc at all.

Target: 1 or 2 fastq filenames (a PairedTarget run's two mates, or a
SingleTarget's one file), relative to the current directory, same as every
other check here.

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: FastQC ran on every mate and no module reported FAIL.
  1: at least one FastQC module reported FAIL on any mate -- the only
     tool-native verdict among the four merged checks (Kraken2, the
     ENA-registry aggregate comparison, and MultiQC are all INFO-only, same
     as their original standalone checks).
  2: FastQC or Kraken2 itself couldn't run on a mate, the Kraken2 database
     wasn't staged, the aggregate comparison's inputs were unreadable, or
     MultiQC crashed -- operational failures, not verdicts about the data.

Metrics (one JSON object on stdout): per-mate FastQC summary + modules,
per-mate Kraken2 classification, the ENA-registry aggregate comparisons
(omitted if --ena-metadata-json isn't given), and the MultiQC report
filename.

KNOWN DEVIATION: MultiQC here only ever aggregates this one run's own 1-2
mate files' FastQC/Kraken2 output, not every sample in a whole work
directory the way `run_scrnaseq_qc.py`'s own local driver runs it. A DCQC
test is inherently per-target, so cross-sample aggregation isn't reachable
from inside one.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from fastq_qc import FASTQC_INSTALL_HINT, aggregate_fastq_checks, fastqc_available, run_fastqc_and_summarize  # noqa: E402
from kraken_qc import KRAKEN2_INSTALL_HINT, kraken2_available, kraken2_db_staged, qc_kraken2  # noqa: E402
from multiqc_qc import run_multiqc  # noqa: E402
from scrnaseq_common import DCQC_ERROR, DCQC_FAIL, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402

# Baked into this image's final stage -- see the Dockerfile.
DEFAULT_KRAKEN2_DB = "/opt/kraken2-db"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq_paths", nargs="+",
                    help="1 (SingleTarget) or 2 (PairedTarget) fastq filenames, relative to the current directory")
    p.add_argument("--ena-metadata-json", default=None,
                    help="ena_metadata.json filename produced by fastq_resolution, for a run resolved from "
                         "SRA/ENA. Omit for a directly-uploaded fastq with no ENA registry to compare against "
                         "-- the aggregate comparison is skipped in that case.")
    p.add_argument("--any-capped", action="store_true",
                    help="Set if any mate file's download was byte-capped for spot-QC")
    p.add_argument("--db-path", default=DEFAULT_KRAKEN2_DB,
                    help=f"Path to the Kraken2 database (default: {DEFAULT_KRAKEN2_DB}, baked into this image)")
    p.add_argument("--expected-species", default=None)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--report-name", default="multiqc_report")
    args = p.parse_args()

    if len(args.fastq_paths) > 2:
        emit_dcqc_reason(
            f"Expected 1 or 2 fastq files (py-dcqc Single/PairedTarget), got {len(args.fastq_paths)}"
        )
        return DCQC_ERROR

    if not fastqc_available():
        emit_dcqc_reason(FASTQC_INSTALL_HINT)
        return DCQC_ERROR
    if not kraken2_available():
        emit_dcqc_reason(KRAKEN2_INSTALL_HINT)
        return DCQC_ERROR
    if not kraken2_db_staged(Path(args.db_path)):
        emit_dcqc_reason(
            f"No Kraken2 database found at {args.db_path} (expected hash.k2d/opts.k2d/taxo.k2d)"
        )
        return DCQC_ERROR

    fastqc_dir = Path("fastqc_out")
    kraken2_dir = Path("kraken2_out")

    per_mate_metrics = []
    per_file_summaries = []
    module_fail_details = []

    for fastq_path_str in args.fastq_paths:
        fastq_path = Path(fastq_path_str)
        fname = fastq_path.name

        try:
            fastqc_findings, summary = run_fastqc_and_summarize(
                fname, fname, fastq_path, fastqc_dir, label=fname,
            )
        except Exception as exc:
            emit_dcqc_reason(f"{fname}: FastQC failed: {exc}")
            return DCQC_ERROR
        if summary is None:
            reason = next(
                (f["metric"] for f in fastqc_findings if f["check"] == "fastqc"), f"{fname}: FastQC failed"
            )
            emit_dcqc_reason(reason)
            return DCQC_ERROR
        per_file_summaries.append(summary)

        modules = {
            f["check"].split(":", 1)[1]: f["status"]
            for f in fastqc_findings if f["check"].startswith("fastqc_module:")
        }
        module_details = {
            f["check"].split(":", 1)[1]: f["metric"]
            for f in fastqc_findings if f["check"].startswith("fastqc_module:")
        }
        failed_modules = [name for name, status in modules.items() if status == "FAIL"]
        if failed_modules:
            details = "; ".join(
                f["metric"] for f in fastqc_findings
                if f["check"].startswith("fastqc_module:") and f["status"] == "FAIL"
            )
            module_fail_details.append(f"{fname}: FAIL module(s) {', '.join(failed_modules)} -- {details}")

        kraken2_findings = qc_kraken2(
            fname, fname, fastq_path, args.db_path, kraken2_dir,
            expected_species=args.expected_species, threads=args.threads,
        )
        by_kraken2_check = {f["check"]: f for f in kraken2_findings}
        if "kraken2" in by_kraken2_check:
            # qc_kraken2()'s own "kraken2" finding is its missing-binary/crash
            # guard -- an operational failure, never a species-match verdict.
            emit_dcqc_reason(f"{fname}: {by_kraken2_check['kraken2']['metric']}")
            return DCQC_ERROR

        mate_metrics = {
            "fname": fname,
            "fastqc": {
                "total_sequences": summary.get("total_sequences"),
                "sequence_length": summary.get("sequence_length"),
                "pct_gc": summary.get("pct_gc"),
                "modules": modules,
                "module_details": module_details,
            },
            "kraken2": {
                "unclassified": by_kraken2_check.get("kraken2_unclassified", {}).get("metric"),
                "top_species": by_kraken2_check.get("kraken2_top_species", {}).get("metric"),
            },
        }
        if "kraken2_species_match" in by_kraken2_check:
            mate_metrics["kraken2"]["species_match"] = by_kraken2_check["kraken2_species_match"]["metric"]
        per_mate_metrics.append(mate_metrics)

    aggregate_metrics = {}
    if args.ena_metadata_json:
        try:
            with open(args.ena_metadata_json) as fh:
                ena_metadata = json.load(fh)
        except Exception as exc:
            emit_dcqc_reason(f"Could not read {args.ena_metadata_json}: {exc}")
            return DCQC_ERROR
        aggregate_findings = aggregate_fastq_checks(
            "aggregate", "aggregate", per_file_summaries, ena_metadata.get("links", []), args.any_capped,
        )
        by_aggregate_check = {f["check"]: f["metric"] for f in aggregate_findings}
        for check in ("read_count_vs_ena", "paired_fastq_parity", "mean_length_vs_registry"):
            if check in by_aggregate_check:
                aggregate_metrics[check] = by_aggregate_check[check]

    try:
        report_path = run_multiqc(Path("."), Path("."), report_name=args.report_name)
    except Exception as exc:
        emit_dcqc_reason(f"MultiQC failed: {exc}")
        return DCQC_ERROR

    emit_dcqc_metrics({
        "mates": per_mate_metrics,
        "aggregate": aggregate_metrics,
        "multiqc_report_file": report_path.name,
    })

    if module_fail_details:
        emit_dcqc_reason("; ".join(module_fail_details))
        return DCQC_FAIL

    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
