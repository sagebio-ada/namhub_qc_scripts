#!/usr/bin/env python3
"""Dockerizable check: cross-check curated Synapse annotations against a
sample's own GEO record.

Wraps geo_metadata_qc.check_geo_annotations(), plus the geo_lookup
presence/absence check that gates it. That gating logic used to live inline
in run_scrnaseq_qc.py's main() -- it's moved here because it's part of this
check's own verdict (found vs. not found in GEO), not orchestrator glue.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from geo_metadata_qc import check_geo_annotations  # noqa: E402
from qc_common import make_finding, write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--entity-id", required=True)
    p.add_argument("--sample-name", required=True, help="GSM id / Synapse entity name to look up in GEO")
    p.add_argument("--annotations-json", required=True, help="JSON file of the Synapse entity's annotations dict")
    p.add_argument("--geo-metadata-json", required=True, help="geo_metadata.json produced by geo_metadata_fetch")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    with open(args.annotations_json) as fh:
        annotations = json.load(fh)
    with open(args.geo_metadata_json) as fh:
        geo_metadata = json.load(fh)

    geo_accession = geo_metadata.get("geo_accession", "")
    geo_row = geo_metadata.get("samples", {}).get(args.sample_name)

    findings = []
    if geo_row:
        findings.append(make_finding(
            args.entity_id, args.sample_name, "geo_lookup", "INFO",
            f"{args.sample_name} found in GEO series {geo_accession}",
        ))
        findings.extend(check_geo_annotations(args.entity_id, args.sample_name, annotations, geo_row))
    else:
        # A genuine blocker, not a judgment call -- there's no GEO row at
        # all to compare against for this sample.
        findings.append(make_finding(
            args.entity_id, args.sample_name, "geo_lookup", "FAIL",
            f"{args.sample_name} not found in GEO series {geo_accession}",
        ))

    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
