#!/usr/bin/env python3
"""Dockerizable check: resolve an SRA run accession to ENA's directly-
downloadable fastq(.gz) link(s) plus their registered size/md5/read_count/
base_count. Wraps fastq_qc.resolve_ena_links(). Network-only; writes
ena_metadata.json for fastq_download and fastq_aggregate to consume.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from fastq_qc import build_session, resolve_ena_links  # noqa: E402
from qc_common import write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--entity-id", required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--srr-accession", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--ena-metadata-out", required=True)
    args = p.parse_args()

    session = build_session()
    findings, links = resolve_ena_links(args.entity_id, args.sample_name, args.srr_accession, session)
    write_findings_json(findings, Path(args.output))

    out_path = Path(args.ena_metadata_out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump({"srr_accession": args.srr_accession, "links": links}, fh, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
