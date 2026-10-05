#!/usr/bin/env python3
"""Dockerizable check: download one ENA fastq link and, for a full (uncapped)
download, verify its size/md5 against the ENA registry row. Wraps
fastq_qc.download_and_check_fastq() -- download and its own-registry
integrity check stay fused (sequential, same file, no benefit to splitting).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from fastq_qc import build_session, download_and_check_fastq  # noqa: E402
from qc_common import write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--entity-id", required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--url", required=True)
    p.add_argument("--dest", required=True, help="Path (on the shared work volume) to write the downloaded fastq to")
    p.add_argument("--registered-size", default="", help="ENA's registered file_size for this file, if known")
    p.add_argument("--registered-md5", default="", help="ENA's registered md5 for this file, if known")
    p.add_argument("--max-download-mb", type=float, default=None)
    p.add_argument("--keep-files", action="store_true",
                    help="Don't delete the destination if a capped download can't be sanitized")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    session = build_session()
    dest = Path(args.dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    max_bytes = int(args.max_download_mb * 1024 * 1024) if args.max_download_mb else None
    registered = {"file_size": args.registered_size, "md5": args.registered_md5}

    findings, _dest_path, _was_capped = download_and_check_fastq(
        args.entity_id, args.sample_name, session, args.url, dest, registered,
        max_bytes=max_bytes, keep_files=args.keep_files,
    )
    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
