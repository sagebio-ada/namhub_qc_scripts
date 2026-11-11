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
per-mate Kraken2 classification (species is reported, never pass/fail),
the R1/R2 read-count parity, the ENA-registry aggregate comparisons
(optional -- skipped unless --ena-metadata-json is given, which py-dcqc
can't supply), the MultiQC report filename, and MultiQC's General
Statistics table. The MultiQC HTML itself is only kept if nf-dcqc
publishes it as an output.

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
from kraken_qc import (  # noqa: E402
    KRAKEN2_INSTALL_HINT, kraken2_available, kraken2_db_staged, kraken2_report_path, parse_kraken2_report, qc_kraken2,
)
from multiqc_qc import read_general_stats, run_multiqc  # noqa: E402
from scrnaseq_common import DCQC_ERROR, DCQC_FAIL, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402

# Baked into this image's final stage -- see the Dockerfile.
DEFAULT_KRAKEN2_DB = "/opt/kraken2-db"

_TOP_N_SPECIES = 5


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


def _kraken2_metrics(fastq_path: Path, kraken2_dir: Path, expected_species):
    """Numeric Kraken2 metrics from the --report file qc_kraken2() wrote."""
    rows = parse_kraken2_report(kraken2_report_path(fastq_path, kraken2_dir))
    unclassified = next((r for r in rows if r["rank"] == "U"), None)
    species = sorted((r for r in rows if r["rank"] == "S"), key=lambda r: -r["pct"])
    metrics = {
        "unclassified_pct": unclassified["pct"] if unclassified else 0.0,
        "top_species": [{"name": r["name"], "pct": r["pct"]} for r in species[:_TOP_N_SPECIES]],
    }
    if expected_species:
        # Reported only -- the species never affects pass/fail.
        match = next((r for r in species if r["name"].lower() == expected_species.lower()), None)
        metrics["expected_species"] = expected_species
        metrics["expected_species_pct"] = match["pct"] if match else 0.0
    return metrics


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
        failed_modules = [name for name, status in modules.items() if status == "FAIL"]
        if failed_modules:
            module_fail_details.append(f"{fname}: FastQC FAIL module(s): {', '.join(failed_modules)}")

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

        per_mate_metrics.append({
            "fname": fname,
            "fastqc": {
                "total_sequences": _number(summary.get("total_sequences")),
                "sequence_length": _number(summary.get("sequence_length")),
                "pct_gc": _number(summary.get("pct_gc")),
                "modules": modules,
            },
            "kraken2": _kraken2_metrics(fastq_path, kraken2_dir, args.expected_species),
        })

    # The R1/R2 parity check needs only the FastQC summaries, so it always
    # runs. The two ENA-registry comparisons are optional: they need
    # --ena-metadata-json, which py-dcqc has no way to supply (it stages only
    # the target files), so under py-dcqc they're always skipped.
    ena_links = []
    if args.ena_metadata_json:
        try:
            with open(args.ena_metadata_json) as fh:
                ena_links = json.load(fh).get("links", [])
        except Exception as exc:
            emit_dcqc_reason(f"Could not read {args.ena_metadata_json}: {exc}")
            return DCQC_ERROR
    aggregate_findings = aggregate_fastq_checks(
        "aggregate", "aggregate", per_file_summaries, ena_links, args.any_capped,
    )
    by_aggregate_check = {f["check"]: f["metric"] for f in aggregate_findings}
    aggregate_metrics = {
        check: by_aggregate_check[check]
        for check in ("read_count_vs_ena", "paired_fastq_parity", "mean_length_vs_registry")
        if check in by_aggregate_check
    }
    if not args.ena_metadata_json:
        aggregate_metrics["ena_comparison"] = "skipped: no --ena-metadata-json given"

    try:
        report_path = run_multiqc(Path("."), Path("."), report_name=args.report_name)
    except Exception as exc:
        emit_dcqc_reason(f"MultiQC failed: {exc}")
        return DCQC_ERROR

    emit_dcqc_metrics({
        "mates": per_mate_metrics,
        "aggregate": aggregate_metrics,
        "multiqc_report_file": report_path.name,
        # The HTML report itself is a file in the container's working
        # directory; nf-dcqc has to publish it as an output to keep it. Its
        # General Statistics table is carried here so the numbers survive
        # even if the HTML doesn't.
        "multiqc_general_stats": read_general_stats(Path("."), args.report_name),
    })

    if module_fail_details:
        emit_dcqc_reason("; ".join(module_fail_details))
        return DCQC_FAIL

    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
