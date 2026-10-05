#!/usr/bin/env python3
"""Dockerizable check: fetch a GEO series' per-sample SOFT metadata.

Wraps geo_metadata_qc.fetch_geo_metadata_by_gsm(). Writes the resulting
{GSM: row} map as geo_metadata.json (consumed by geo_annotation_check, and by
the driver for the sample_organism_ch1 value it threads into fastqc/kraken2_qc
as --expected-species), plus a Findings JSON reporting whether the fetch
itself succeeded.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from geo_metadata_qc import fetch_geo_metadata_by_gsm  # noqa: E402
from qc_common import make_finding, write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geo-accession", required=True, help="GEO series accession, e.g. GSE293390")
    p.add_argument("--output", required=True, help="Path to write Findings JSON")
    p.add_argument("--geo-metadata-out", required=True, help="Path to write geo_metadata.json")
    args = p.parse_args()

    try:
        samples = fetch_geo_metadata_by_gsm(args.geo_accession)
    except Exception as exc:
        findings = [make_finding(
            args.geo_accession, args.geo_accession, "geo_metadata_fetch", "FAIL",
            f"Could not fetch GEO metadata for {args.geo_accession}: {exc}",
        )]
        write_findings_json(findings, Path(args.output))
        return 0

    findings = [make_finding(
        args.geo_accession, args.geo_accession, "geo_metadata_fetch", "INFO",
        f"{len(samples)} sample(s) found in GEO series {args.geo_accession}",
    )]
    write_findings_json(findings, Path(args.output))

    out_path = Path(args.geo_metadata_out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump({"geo_accession": args.geo_accession, "samples": samples}, fh, indent=2, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main())
