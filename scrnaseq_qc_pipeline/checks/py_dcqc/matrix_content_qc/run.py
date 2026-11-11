#!/usr/bin/env python3
"""CellRanger matrix content QC (file_integrity through mito_fraction).
Wraps matrix_qc.qc_matrix_content() -- kept fused since these all operate on
the same small in-memory scipy sparse matrix, which is cheap to re-parse and
has no clean cross-container serialization.

py-dcqc ExternalBaseTest contract (see
https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96):
  0: the matrix parsed and has at least one usable (non-empty) barcode.
  1: the input is bad: an archive that can't be read or doesn't hold
     exactly one barcodes/features/matrix file, a file that can't be
     read/parsed, or a matrix that is completely degenerate (zero
     non-empty barcodes).
  2: the check couldn't run (a wrong number of arguments, or an input file
     that doesn't exist) -- not a verdict about the data. A failure
     computing the mitochondrial fraction specifically is treated as
     non-fatal (noted in the metrics, not an error) since the core
     integrity/dimensions/count checks already succeeded by that point.

Metrics (one JSON object on stdout): parsed shape, barcode/feature counts,
median UMI/genes per cell, empty-barcode fraction, median mito fraction
(if computable).

Input: either
  - ONE tar archive (.tar.gz/.tgz/.tar) holding barcodes.tsv[.gz],
    features.tsv[.gz] (or CellRanger v2's genes.tsv[.gz]) and
    matrix.mtx[.gz], at any depth -- e.g. GEO's
    <sample>_filtered_feature_bc_matrix.tar.gz, or one made with
    package_matrix.py. This is the py-dcqc form: one file, a SingleTarget.
  - or the three files themselves (barcodes, features, matrix, in that
    order), for running outside py-dcqc.
For an archive, a missing/duplicated member or an unreadable archive is
exit 1 -- the file itself is bad, not an operational failure.

Filenames are relative to the current directory (nf-dcqc stages files
there as symlinks) -- do not pass directories.
"""
from __future__ import annotations

import argparse
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from matrix_qc import extract_matrix_archive, qc_matrix_content  # noqa: E402
from scrnaseq_common import DCQC_ERROR, DCQC_FAIL, DCQC_PASS, emit_dcqc_metrics, emit_dcqc_reason  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paths", nargs="+",
                   help="One matrix tar archive, or the barcodes, features and matrix files (in that order)")
    p.add_argument("--cellranger-output-class", default="")
    args = p.parse_args()

    paths = [Path(x) for x in args.paths]
    if len(paths) not in (1, 3):
        emit_dcqc_reason(f"Expected 1 archive or 3 files (barcodes, features, matrix), got {len(paths)}")
        return DCQC_ERROR
    for path in paths:
        if not path.is_file():
            emit_dcqc_reason(f"{path.name}: file not found")
            return DCQC_ERROR

    with tempfile.TemporaryDirectory() as tmp:
        if len(paths) == 1:
            sample_id = paths[0].name
            try:
                extracted = extract_matrix_archive(paths[0], Path(tmp))
            except (ValueError, tarfile.TarError, EOFError, OSError) as exc:
                emit_dcqc_reason(f"{sample_id}: not a usable matrix archive: {exc}")
                return DCQC_FAIL
            barcodes, features, matrix = extracted["barcodes"], extracted["features"], extracted["matrix"]
        else:
            barcodes, features, matrix = paths
            sample_id = matrix.name
        findings = qc_matrix_content(
            sample_id, barcodes, features, matrix,
            cellranger_output_class=args.cellranger_output_class,
        )
    by_check = {f["check"]: f for f in findings}

    for check in ("file_integrity", "matrix_integrity"):
        if check in by_check and by_check[check]["status"] == "FAIL":
            # The files exist (checked above), so unreadable contents mean
            # a bad file, not a check that couldn't run. py-dcqc doesn't
            # read metrics after an error, but does after a fail, so write
            # an empty metrics object to keep stdout valid JSON.
            emit_dcqc_metrics({})
            emit_dcqc_reason(by_check[check]["metric"])
            return DCQC_FAIL

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
