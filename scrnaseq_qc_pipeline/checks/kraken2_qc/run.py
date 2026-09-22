#!/usr/bin/env python3
"""Dockerizable check: species/contamination screening via Kraken2
taxonomic classification. Wraps kraken_qc.qc_kraken2().

Needs a pre-built Kraken2 database directory mounted into the container (see
kraken_qc.py's module docstring for where to get one) and passed via
--db-path.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from kraken_qc import qc_kraken2  # noqa: E402
from qc_common import write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--entity-id", required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--fastq-path", required=True)
    p.add_argument("--db-path", required=True, help="Path (inside the container) to the extracted Kraken2 database")
    p.add_argument("--outdir", required=True, help="Where to write the *.kraken2report.txt")
    p.add_argument("--expected-species", default=None)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    findings = qc_kraken2(
        args.entity_id, args.sample_name, Path(args.fastq_path), args.db_path, Path(args.outdir),
        expected_species=args.expected_species, threads=args.threads,
    )
    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
