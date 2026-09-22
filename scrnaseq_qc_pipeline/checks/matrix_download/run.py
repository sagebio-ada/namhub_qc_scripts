#!/usr/bin/env python3
"""Dockerizable check: verify a processed CellRanger sample has all three
required roles (barcodes/features/matrix) and download them from Synapse.
Wraps matrix_qc.check_processed_completeness() + matrix_qc.download_matrix_files()
-- fused into one script since they're sequential and cheap; the second is a
no-op if the first fails.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from matrix_qc import check_processed_completeness, download_matrix_files  # noqa: E402
from qc_common import write_findings_json  # noqa: E402

try:
    import synapseclient
except ImportError:
    print("synapseclient is required: pip install synapseclient", file=sys.stderr)
    raise


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--barcodes-entity-id", required=True)
    p.add_argument("--features-entity-id", required=True)
    p.add_argument("--matrix-entity-id", required=True)
    p.add_argument("--work-dir", required=True, help="Shared work directory; files land under <work-dir>/<sample-id>/")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    files = {
        "barcodes": args.barcodes_entity_id,
        "features": args.features_entity_id,
        "matrix": args.matrix_entity_id,
    }

    findings, complete = check_processed_completeness(args.sample_id, files)
    if complete:
        syn = synapseclient.Synapse()
        syn.login(authToken=os.environ.get("SYNAPSE_AUTH_TOKEN"), silent=True)
        dl_findings, _paths = download_matrix_files(syn, args.sample_id, files, Path(args.work_dir))
        findings.extend(dl_findings)

    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
