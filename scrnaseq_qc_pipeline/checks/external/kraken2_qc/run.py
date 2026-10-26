#!/usr/bin/env python3
"""Species/contamination screening via Kraken2 taxonomic classification.
Wraps kraken_qc.qc_kraken2(). Has no validated tool-native verdict of its
own -- Kraken2 doesn't assert what counts as "too much" unclassified or
contamination for RNA-seq, so this check never resolves to a real FAIL.

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: classification completed. The numbers are in the metrics; there is no
     verdict attached to them (consistent with this pipeline's standing
     policy against inventing an unvalidated threshold).
  2: the check could not run at all -- the kraken2 binary is missing, the
     database could not be staged, or classification crashed.
(Exit code 1 is defined by the contract but genuinely unused by this
check -- there is no condition under which it reports FAIL.)

KNOWN DEVIATION from py-dcqc's "What Makes a Good Test" rule ("the tool must
not download data when it runs" / "put all the reference data into the
image"): this auto-downloads the database on demand (see ensure_kraken2_db
below) instead of shipping it baked into the image. See the Dockerfile for
why, and the fully-compliant alternative.

Metrics (one JSON object on stdout): unclassified %, top species, and the
expected-species match % if --expected-species was given. Working directory
stays clean -- Kraken2's own report file is written to a throwaway temp
directory, not the task's own working directory.

Runs on demand: if --db-path isn't already a staged Kraken2 database, it's
downloaded and extracted there automatically before classification (Kraken2's
own Standard-8 build by default, or --db-url). Mount a *persistent* volume
at --db-path (a host directory or named Docker volume, not the container's
own ephemeral filesystem) so that download only ever happens once, not on
every run:

    docker run --rm -v kraken2-db-cache:/db -v $(pwd):$(pwd) -w $(pwd) \\
        qc-kraken2 --db-path /db sample.fastq.gz

Invoked with the fastq's filename, relative to the current directory
(nf-dcqc stages the file there as a symlink) -- do not pass a directory.
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from kraken_qc import ensure_kraken2_db, kraken2_db_staged, qc_kraken2  # noqa: E402
from qc_common import DCQC_ERROR, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq_path", help="Fastq filename, relative to the current directory")
    p.add_argument("--db-path", required=True,
                    help="Path to the Kraken2 database (ideally a mounted persistent volume). "
                         "Auto-downloaded here on first use if not already staged.")
    p.add_argument("--db-url", default=None,
                    help="Override the database tarball URL to auto-download (default: Kraken2's own "
                         "Standard-8 build). Ignored if --db-path is already staged.")
    p.add_argument("--expected-species", default=None)
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args()

    fastq_path = Path(args.fastq_path)
    db_path = Path(args.db_path)

    db_findings = ensure_kraken2_db(db_path, url=args.db_url)
    if not kraken2_db_staged(db_path):
        emit_dcqc_reason(db_findings[0]["metric"])
        return DCQC_ERROR

    with tempfile.TemporaryDirectory() as tmp:
        findings = qc_kraken2(
            fastq_path.name, fastq_path.name, fastq_path, str(db_path), Path(tmp),
            expected_species=args.expected_species, threads=args.threads,
        )

    by_check = {f["check"]: f for f in findings}
    if "kraken2" in by_check:
        # qc_kraken2()'s own "kraken2" finding is its missing-binary/crash
        # guard -- an operational failure, never a species-match verdict.
        emit_dcqc_reason(by_check["kraken2"]["metric"])
        return DCQC_ERROR

    metrics = {
        "fname": fastq_path.name,
        "unclassified": by_check.get("kraken2_unclassified", {}).get("metric"),
        "top_species": by_check.get("kraken2_top_species", {}).get("metric"),
    }
    if "kraken2_species_match" in by_check:
        metrics["species_match"] = by_check["kraken2_species_match"]["metric"]
    emit_dcqc_metrics(metrics)
    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
