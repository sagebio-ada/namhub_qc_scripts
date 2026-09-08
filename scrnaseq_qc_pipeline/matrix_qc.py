"""Processed-data QC for CellRanger-style output (barcodes/features/matrix.mtx)."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from qc_common import Finding, make_finding

REQUIRED_ROLES = {"barcodes", "features", "matrix"}


def _count_gz_lines(path: Path) -> int:
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh)


def _mito_fraction(features_path: Path, matrix) -> Optional[np.ndarray]:
    """Fraction of counts per cell coming from MT- genes, or None if none found."""
    mito_rows = []
    with gzip.open(features_path, "rt") as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            symbol = parts[1] if len(parts) > 1 else parts[0]
            if symbol.upper().startswith("MT-"):
                mito_rows.append(i)
    if not mito_rows:
        return None
    total = np.asarray(matrix.sum(axis=0)).ravel()
    mito = np.asarray(matrix[mito_rows, :].sum(axis=0)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(total > 0, mito / total, 0.0)


def qc_processed_sample(
    syn,
    sample_id: str,
    files: Dict[str, str],
    work_dir: Path,
    cellranger_output_class: str = "",
) -> List[Finding]:
    findings: List[Finding] = []

    def add(check: str, status: str, detail: str) -> None:
        findings.append(make_finding(sample_id, sample_id, check, status, detail))

    missing = REQUIRED_ROLES - files.keys()
    if missing:
        add("processed_completeness", "FAIL", f"Missing role(s) {sorted(missing)} for sample {sample_id}")
        return findings
    present = ", ".join(f"{role}={entity_id}" for role, entity_id in sorted(files.items()))
    add("processed_completeness", "PASS", f"barcodes/features/matrix all present: {present}")

    sample_dir = work_dir / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, Path] = {}
    for role, entity_id in files.items():
        try:
            ent = syn.get(entity_id, downloadLocation=str(sample_dir), ifcollision="overwrite.local")
            paths[role] = Path(ent.path)
        except Exception as exc:
            add(f"download:{role}", "FAIL", f"{entity_id}: {exc}")
    if len(paths) < 3:
        return findings

    try:
        n_barcodes = _count_gz_lines(paths["barcodes"])
        n_features = _count_gz_lines(paths["features"])
    except Exception as exc:
        add("file_integrity", "FAIL", f"Could not read barcodes/features files: {exc}")
        return findings
    add("file_integrity", "PASS", f"{n_barcodes} barcodes, {n_features} features readable")

    try:
        import scipy.io
        with gzip.open(paths["matrix"], "rb") as fh:
            mat = scipy.io.mmread(fh).tocsc()
    except Exception as exc:
        add("matrix_integrity", "FAIL", f"Could not parse matrix.mtx.gz: {exc}")
        return findings

    n_genes_mtx, n_cells_mtx = mat.shape
    if n_genes_mtx != n_features or n_cells_mtx != n_barcodes:
        add("matrix_dimensions", "FAIL",
            f"matrix header {mat.shape} != (features={n_features}, barcodes={n_barcodes})")
    else:
        add("matrix_dimensions", "PASS", f"matrix shape {mat.shape} matches barcodes/features counts")

    counts_per_cell = np.asarray(mat.sum(axis=0)).ravel()
    genes_per_cell = np.asarray((mat > 0).sum(axis=0)).ravel()
    nonzero_counts = counts_per_cell[counts_per_cell > 0]
    nonzero_genes = genes_per_cell[genes_per_cell > 0]
    # INFO for the normal case (just reporting numbers); WARN is a real flag --
    # zero non-empty barcodes in the whole matrix is a genuine integrity problem,
    # not a threshold call.
    add("matrix_metrics", "INFO" if len(nonzero_counts) else "WARN",
        f"median UMI/cell={np.median(nonzero_counts) if len(nonzero_counts) else 0:.0f}, "
        f"median genes/cell={np.median(nonzero_genes) if len(nonzero_genes) else 0:.0f}, "
        f"barcodes_with_zero_counts={(counts_per_cell == 0).sum()}/{len(counts_per_cell)}, nnz={mat.nnz}")

    # No pass/fail judgment -- what fraction of empty barcodes is "expected"
    # for a raw vs. filtered matrix depends on the dataset and isn't something
    # we've validated a cutoff for. Report the fraction and the declared class
    # side by side; let whoever reads the report judge whether they agree.
    frac_empty = (counts_per_cell == 0).sum() / max(len(counts_per_cell), 1)
    add("empty_barcode_check", "INFO",
        f"{frac_empty:.1%} of {len(counts_per_cell)} barcodes have zero counts "
        f"(declared CellrangerOutputClass={cellranger_output_class!r})")

    try:
        mito_frac = _mito_fraction(paths["features"], mat)
        if mito_frac is not None:
            nz_mito = mito_frac[counts_per_cell > 0]
            if len(nz_mito):
                add("mito_fraction", "INFO", f"median mito fraction (non-empty barcodes)={np.median(nz_mito):.2%}")
    except Exception as exc:
        add("mito_fraction", "WARN", f"Could not compute mitochondrial fraction: {exc}")

    return findings
