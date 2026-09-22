#!/usr/bin/env python3
"""Dockerizable check: species/contamination screening via Kraken2
taxonomic classification. Wraps kraken_qc.qc_kraken2().

Runs on demand: if --db-path isn't already a staged Kraken2 database, it's
downloaded and extracted there automatically before classification (Kraken2's
own Standard-8 build by default, or --db-url). Mount a *persistent* volume
at --db-path (a host directory or named Docker volume, not the container's
own ephemeral filesystem) so that download only ever happens once, not on
every run:

    docker run --rm -v kraken2-db-cache:/db -v $(pwd)/work:/data qc-kraken2 \\
        --db-path /db --fastq-path /data/sample.fastq.gz ...
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from kraken_qc import ensure_kraken2_db, kraken2_db_staged, qc_kraken2  # noqa: E402
from qc_common import write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--entity-id", required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--fastq-path", required=True)
    p.add_argument("--db-path", required=True,
                    help="Path (inside the container, ideally a mounted persistent volume) to the Kraken2 "
                         "database. Auto-downloaded here on first use if not already staged.")
    p.add_argument("--db-url", default=None,
                    help="Override the database tarball URL to auto-download (default: Kraken2's own "
                         "Standard-8 build). Ignored if --db-path is already staged.")
    p.add_argument("--outdir", required=True, help="Where to write the *.kraken2report.txt")
    p.add_argument("--expected-species", default=None)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    db_path = Path(args.db_path)
    findings = list(ensure_kraken2_db(db_path, url=args.db_url))
    if kraken2_db_staged(db_path):
        findings.extend(qc_kraken2(
            args.entity_id, args.sample_name, Path(args.fastq_path), str(db_path), Path(args.outdir),
            expected_species=args.expected_species, threads=args.threads,
        ))
    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
