#!/usr/bin/env python3
"""CellRanger matrix content QC (file_integrity through mito_fraction).
Wraps matrix_qc.qc_matrix_content() -- kept fused since these all operate on
the same small in-memory scipy sparse matrix, which is cheap to re-parse and
has no clean cross-container serialization.

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: the matrix parsed and has at least one usable (non-empty) barcode.
  1: the matrix parsed but is completely degenerate (zero non-empty
     barcodes) -- a genuine tool-native bad-data verdict.
  2: one of the three files couldn't be read/parsed at all -- an
     operational failure, not a verdict about the data. A failure
     computing the mitochondrial fraction specifically is treated as
     non-fatal (noted in the metrics, not an error) since the core
     integrity/dimensions/count checks already succeeded by that point.

Metrics (one JSON object on stdout): parsed shape, barcode/feature counts,
median UMI/genes per cell, empty-barcode fraction, median mito fraction
(if computable).

Note: this needs three input files (barcodes/features/matrix), beyond
py-dcqc's currently-tested Target cardinality (SingleTarget=1 file,
PairedTarget=2) -- "No test does this today" per the architecture doc, so
a 3-file external test is untested territory in the actual nf-dcqc runtime.

Invoked with three filenames, relative to the current directory (nf-dcqc
stages files there as symlinks) -- do not pass directories.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from matrix_qc import qc_matrix_content  # noqa: E402
from qc_common import DCQC_ERROR, DCQC_FAIL, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("barcodes_path", help="barcodes.tsv.gz filename, relative to the current directory")
    p.add_argument("features_path", help="features.tsv.gz filename, relative to the current directory")
    p.add_argument("matrix_path", help="matrix.mtx.gz filename, relative to the current directory")
    p.add_argument("--cellranger-output-class", default="")
    args = p.parse_args()

    sample_id = Path(args.matrix_path).name
    findings = qc_matrix_content(
        sample_id, Path(args.barcodes_path), Path(args.features_path), Path(args.matrix_path),
        cellranger_output_class=args.cellranger_output_class,
    )
    by_check = {f["check"]: f for f in findings}

    for check in ("file_integrity", "matrix_integrity"):
        if check in by_check and by_check[check]["status"] == "FAIL":
            emit_dcqc_reason(by_check[check]["metric"])
            return DCQC_ERROR

    metrics = {
        "matrix_integrity": by_check.get("matrix_integrity", {}).get("metric"),
        "matrix_dimensions": by_check.get("matrix_dimensions", {}).get("metric"),
        "matrix_metrics": by_check.get("matrix_metrics", {}).get("metric"),
        "empty_barcode_check": by_check.get("empty_barcode_check", {}).get("metric"),
    }
    if "mito_fraction" in by_check:
        metrics["mito_fraction"] = by_check["mito_fraction"]["metric"]
    emit_dcqc_metrics(metrics)

    if by_check.get("matrix_metrics", {}).get("status") == "FAIL":
        # Completely degenerate (zero non-empty barcodes) -- a real
        # tool-native verdict, not a crash.
        emit_dcqc_reason(by_check["matrix_metrics"]["metric"])
        return DCQC_FAIL

    return DCQC_PASS


if __name__ == "__main__":
    sys.exit(main())
