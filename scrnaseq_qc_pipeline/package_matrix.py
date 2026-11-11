#!/usr/bin/env python3
"""Package one sample's barcodes/features/matrix files into a single
.tar.gz, so the matrix content check can run on it as one py-dcqc file.

The files go into the archive under their standard CellRanger names
(barcodes.tsv[.gz], features.tsv[.gz], matrix.mtx[.gz]), whatever they are
called on disk -- GEO usually prefixes them, e.g. GSM123_barcodes.tsv.gz.
Each file is stored as-is; a .gz file stays gzipped inside the archive.

Example:
  python3 package_matrix.py GSM123_barcodes.tsv.gz GSM123_features.tsv.gz \\
      GSM123_matrix.mtx.gz -o GSM123_matrix.tar.gz

An archive GEO already provides (e.g. <sample>_filtered_feature_bc_matrix.tar.gz)
doesn't need repackaging: the check finds the three files in it by name.
"""
from __future__ import annotations

import argparse
import sys
import tarfile
from pathlib import Path


def package_matrix(barcodes: Path, features: Path, matrix: Path, out_path: Path) -> Path:
    """Write barcodes/features/matrix into out_path (.tar.gz) under their
    standard names. Returns out_path."""
    members = (("barcodes.tsv", barcodes), ("features.tsv", features), ("matrix.mtx", matrix))
    for _, path in members:
        if not path.is_file():
            raise FileNotFoundError(f"{path}: file not found")
    with tarfile.open(out_path, "w:gz") as tar:
        for standard_name, path in members:
            suffix = ".gz" if path.name.endswith(".gz") else ""
            tar.add(path, arcname=standard_name + suffix, recursive=False)
    return out_path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("barcodes", type=Path, help="barcodes.tsv[.gz]")
    p.add_argument("features", type=Path, help="features.tsv[.gz] (or CellRanger v2 genes.tsv[.gz])")
    p.add_argument("matrix", type=Path, help="matrix.mtx[.gz]")
    p.add_argument("-o", "--output", type=Path, required=True, help="Output archive, e.g. <sample>_matrix.tar.gz")
    args = p.parse_args()

    if not args.output.name.endswith((".tar.gz", ".tgz")):
        print(f"Output name must end in .tar.gz or .tgz: {args.output}", file=sys.stderr)
        return 2
    try:
        package_matrix(args.barcodes, args.features, args.matrix, args.output)
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 2
    print(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
