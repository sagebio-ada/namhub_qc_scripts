#!/usr/bin/env python3
"""Dockerizable check: download a directly-uploaded (no SRA indirection)
fastq from Synapse. Wraps fastq_qc.download_from_synapse().

This is the qc_direct_fastq() analog of fastq_download -- no ENA registry to
check integrity against here, so it's download-only.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from fastq_qc import download_from_synapse  # noqa: E402
from qc_common import write_findings_json  # noqa: E402

try:
    import synapseclient
except ImportError:
    print("synapseclient is required: pip install synapseclient", file=sys.stderr)
    raise


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--entity-id", required=True, help="Synapse file entity id of the uploaded fastq")
    p.add_argument("--sample-name", required=True)
    p.add_argument("--work-dir", required=True, help="Shared work directory; the file lands under <work-dir>/<sample-name>/")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    syn = synapseclient.Synapse()
    syn.login(authToken=os.environ.get("SYNAPSE_AUTH_TOKEN"), silent=True)

    findings, _fastq_path = download_from_synapse(syn, args.entity_id, args.sample_name, Path(args.work_dir))
    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
