#!/usr/bin/env python3
"""Dockerizable check: CellRanger matrix content QC (file_integrity through
mito_fraction). Wraps matrix_qc.qc_matrix_content() -- kept fused since these
all operate on the same small in-memory scipy sparse matrix, which is cheap
to re-parse and has no clean cross-container serialization.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from matrix_qc import qc_matrix_content  # noqa: E402
from qc_common import write_findings_json  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--barcodes-path", required=True)
    p.add_argument("--features-path", required=True)
    p.add_argument("--matrix-path", required=True)
    p.add_argument("--cellranger-output-class", default="")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    findings = qc_matrix_content(
        args.sample_id, Path(args.barcodes_path), Path(args.features_path), Path(args.matrix_path),
        cellranger_output_class=args.cellranger_output_class,
    )
    write_findings_json(findings, Path(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
