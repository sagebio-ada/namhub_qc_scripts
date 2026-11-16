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

import json
import shutil
import subprocess
from pathlib import Path
from typing import Dict

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


# Plot config keys worth keeping: they say what a plot and each of its
# series are. Everything else in a plot's config (colors, axis limits,
# tooltips) only affects how the HTML draws it.
_PLOT_CONFIG_KEYS = ("id", "title", "xlab", "ylab", "data_labels")


def _short_general_stats_key(column: str) -> str:
    """"FastQC_mqc-generalstats-fastqc-percent_gc" -> "fastqc.percent_gc"."""
    module, _, rest = column.partition("_mqc-generalstats-")
    return f"{module.lower()}.{rest.rsplit('-', 1)[-1]}" if rest else column


def _strip_colors(value):
    """Drop "color"/"colors" keys at any depth of a plot's data."""
    if isinstance(value, dict):
        return {k: _strip_colors(v) for k, v in value.items() if k not in ("color", "colors")}
    if isinstance(value, list):
        return [_strip_colors(v) for v in value]
    return value


def read_multiqc_data(output_dir: Path, report_name: str = "multiqc_report") -> Dict[str, object]:
    """Every number in the MultiQC report, from the multiqc_data.json that
    MultiQC writes next to the HTML (<output_dir>/<report_name>_data/).

    Returns {"version", "tables", "plots"}:
      - "tables": every table MultiQC saved, as {table: {sample: {column:
        value}}} -- e.g. the full FastQC table (multiqc_fastqc) and General
        Statistics (multiqc_general_stats, with keys shortened to
        "<module>.<metric>").
      - "plots": the data behind every plot, as {plot_id: {"plot_type",
        "datasets", "samples" (bar graphs only), "config"}}, with only the
        config keys that label the data (_PLOT_CONFIG_KEYS) and no colors.
    Returns {} if MultiQC wrote no data file (it found no supported tool
    output).
    """
    path = output_dir / f"{report_name}_data" / "multiqc_data.json"
    if not path.exists():
        return {}
    with open(path) as fh:
        data = json.load(fh)

    tables = dict(data.get("report_saved_raw_data", {}))
    if "multiqc_general_stats" in tables:
        tables["multiqc_general_stats"] = {
            sample: {_short_general_stats_key(col): val for col, val in row.items()}
            for sample, row in tables["multiqc_general_stats"].items()
        }

    plots = {}
    for plot_id, plot in data.get("report_plot_data", {}).items():
        config = plot.get("config", {})
        plots[plot_id] = _strip_colors({
            **{k: v for k, v in plot.items() if k != "config"},
            "config": {k: config[k] for k in _PLOT_CONFIG_KEYS if k in config},
        })

    return {"version": data.get("config_version"), "tables": tables, "plots": plots}
