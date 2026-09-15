"""Optional aggregated HTML report via MultiQC.

MultiQC doesn't compute any new QC metrics -- it scans a directory for
already-existing output from other tools it recognizes and renders one
interactive HTML report comparing all samples. In this pipeline that means:

  - FastQC: MultiQC's `fastqc` module reads the same fastqc_data.txt files
    fastq_qc.py already writes into <work_dir>/<sample>/fastqc/.
  - Kraken2: MultiQC's `kraken` module reads the same *.kraken2report.txt
    files kraken_qc.py already writes into <work_dir>/<sample>/kraken2/.

It has no module for `fq lint`, and its `cellranger` module needs
CellRanger's own web_summary.html run report, which isn't part of this
dataset -- only the final matrix files were deposited, so that module has
nothing to parse here. It also knows nothing about this pipeline's own
checks (the GEO/ENA metadata cross-checks, matrix sanity checks) --
MultiQC is a visualization layer on top of the FastQC/Kraken2 portion only,
never a replacement for scrnaseq_qc_report.csv.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

MULTIQC_INSTALL_HINT = (
    "multiqc was not found on PATH. Install it with 'pip install multiqc'. "
    "pip sometimes installs its console script to a user directory that "
    "isn't on PATH (e.g. ~/Library/Python/3.x/bin on macOS) -- if so, add "
    "that directory to PATH or invoke multiqc there directly."
)


def multiqc_available() -> bool:
    return shutil.which("multiqc") is not None


def run_multiqc(work_dir: Path, output_dir: Path, report_name: str = "multiqc_report") -> Path:
    """Scan work_dir for FastQC/Kraken2 output and produce one aggregated HTML report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["multiqc", str(work_dir), "-o", str(output_dir), "-n", report_name, "-f"],
        check=True, capture_output=True, text=True,
    )
    report_path = output_dir / f"{report_name}.html"
    if not report_path.exists():
        raise FileNotFoundError(f"Expected MultiQC report not found: {report_path}")
    return report_path
