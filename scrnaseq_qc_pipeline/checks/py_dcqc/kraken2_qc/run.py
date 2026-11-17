#!/usr/bin/env python3
"""Kraken2 taxonomic classification of each FASTQ file, against the
Standard-8 database baked into this image. Reports only: the species found
never makes the check fail.

Target: 1 or 2 fastq filenames (a PairedTarget's two mates, or a
SingleTarget's one file), relative to the current directory.

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: Kraken2 classified every file. The results are metrics only.
  1: Kraken2 could not process a readable file (e.g. a truncated gzip or a
     corrupt FASTQ) -- the same verdict fq_lint and fastqc give that file.
  2: the check itself could not run: Kraken2 or its database is missing, a
     file is missing or unreadable, or Kraken2 was stopped from outside (a
     signal, or running out of memory -- the database is loaded into
     memory, about 8 GB). Every file is still checked first, and the
     reason also lists any FAIL found on another.

Metrics (one JSON object on stdout, on exit 0 and 1): {"files": [...]},
per file: fname, n_reads, classified_pct, unclassified_pct, top_species
[{name, pct}] (the 5 largest), and expected_species /
expected_species_pct when --expected-species is given. A file Kraken2
couldn't process has "error" instead.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from kraken_qc import (  # noqa: E402
    KRAKEN2_INSTALL_HINT, CorruptGzipError, kraken2_available, kraken2_db_staged, parse_kraken2_report,
    run_kraken2,
)
from scrnaseq_common import DCQC_ERROR, DCQC_FAIL, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402

# Baked into this image -- see the Dockerfile.
DEFAULT_KRAKEN2_DB = "/opt/kraken2-db"

_TOP_N_SPECIES = 5

# Kraken2 messages that mean it ran out of memory (a resource limit, not a
# bad file).
_OUT_OF_MEMORY_MARKERS = ("bad_alloc", "Cannot allocate memory", "out of memory")


def kraken2_metrics(report_path: Path, expected_species=None) -> dict:
    """Numeric metrics from a Kraken2 --report file."""
    rows = parse_kraken2_report(report_path)
    unclassified = next((r for r in rows if r["rank"] == "U"), None)
    root = next((r for r in rows if r["rank"] == "R"), None)
    species = sorted((r for r in rows if r["rank"] == "S"), key=lambda r: -r["pct"])
    metrics = {
        "n_reads": (unclassified["clade_count"] if unclassified else 0) + (root["clade_count"] if root else 0),
        "classified_pct": root["pct"] if root else 0.0,
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
    p.add_argument("--db-path", default=DEFAULT_KRAKEN2_DB,
                   help=f"Path to the Kraken2 database (default: {DEFAULT_KRAKEN2_DB}, baked into this image)")
    p.add_argument("--expected-species", default=None,
                   help="Report the percentage of reads classified as this species (never pass/fail)")
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args()

    if len(args.fastq_paths) > 2:
        emit_dcqc_reason(f"Expected 1 or 2 fastq files (py-dcqc Single/PairedTarget), got {len(args.fastq_paths)}")
        return DCQC_ERROR
    if not kraken2_available():
        emit_dcqc_reason(KRAKEN2_INSTALL_HINT)
        return DCQC_ERROR
    if not kraken2_db_staged(Path(args.db_path)):
        emit_dcqc_reason(f"No Kraken2 database found at {args.db_path} (expected hash.k2d/opts.k2d/taxo.k2d)")
        return DCQC_ERROR

    outdir = Path("kraken2_out")
    per_file_metrics = []
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
            report_path = run_kraken2(fastq_path, args.db_path, outdir, threads=args.threads)
        except CorruptGzipError as exc:
            # A truncated or corrupt .gz -- the same FAIL fq_lint and fastqc
            # give that file. (Kraken2 alone would classify the reads before
            # the break and pass.)
            fail_reasons.append(f"{fname}: {exc}")
            per_file_metrics.append({"fname": fname, "error": str(exc)})
            continue
        except subprocess.CalledProcessError as exc:
            stderr = (exc.stderr or "").strip()
            if exc.returncode < 0:
                error_reasons.append(f"{fname}: Kraken2 was stopped by signal {-exc.returncode}")
            elif any(marker in stderr for marker in _OUT_OF_MEMORY_MARKERS):
                error_reasons.append(f"{fname}: Kraken2 ran out of memory: {stderr[-300:]}")
            else:
                # Kraken2 and its database are present and the file is
                # readable, so Kraken2 giving up on it means the file is bad.
                message = f"Kraken2 could not process the file (exit {exc.returncode}): {stderr[-300:]}"
                fail_reasons.append(f"{fname}: {message}")
                per_file_metrics.append({"fname": fname, "error": message})
            continue
        except Exception as exc:
            error_reasons.append(f"{fname}: Kraken2 did not produce a report: {exc}")
            continue

        per_file_metrics.append({"fname": fname, **kraken2_metrics(report_path, args.expected_species)})

    if error_reasons:
        # py-dcqc doesn't read metrics after an error; the reason keeps
        # every file's FAIL too.
        emit_dcqc_reason("; ".join(error_reasons + fail_reasons))
        return DCQC_ERROR

    emit_dcqc_metrics({"files": per_file_metrics})
    if fail_reasons:
        emit_dcqc_reason("; ".join(fail_reasons))
        return DCQC_FAIL
    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
