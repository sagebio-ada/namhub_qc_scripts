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

import html
import json
import re
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


_SECTION_RE = re.compile(r'<div class="mqc-section mqc-section-[^"]*">(.*?)(?=<div class="mqc-section mqc-section-|</body>)', re.S)
# Data blocks the HTML embeds outside the section that draws them, and the
# section each one belongs to.
_HTML_BLOCK_SECTIONS = {
    "fastqc_seq_content": "fastqc_per_base_sequence_content",
    "fastqc_passfails": "fastqc_status_checks",
}
_HTML_JSON_RE = re.compile(r'<script type="application/json" class="([^"]+)">(.*?)</script>', re.S)


def _html_text(fragment: str) -> str:
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", fragment))).strip()


def read_html_data(html_path: Path) -> Dict[str, object]:
    """Data that only the HTML report holds, not multiqc_data.json:
      - every JSON block a module embeds in the page, by its class -- e.g.
        FastQC's per-base sequence content heatmap (fastqc_seq_content) and
        its per-module status (fastqc_passfails). A [module, data] block is
        stored as its data.
      - every report section, in order: {"id", "title", "data" (the ids of
        the plots/blocks it draws), "notes" (the messages it shows instead
        of, or next to, a plot -- e.g. "All samples have sequences of a
        single length (150bp)." or "No samples found with any adapter
        contamination > 0.1%")}.
    """
    text = html_path.read_text(errors="replace")
    blocks: Dict[str, object] = {}
    for cls, body in _HTML_JSON_RE.findall(text):
        try:
            value = json.loads(body)
        except ValueError:
            continue
        if isinstance(value, list) and len(value) == 2 and isinstance(value[0], str):
            value = value[1]
        blocks[cls] = value

    sections = []
    for body in _SECTION_RE.findall(text):
        heading = re.search(r'<h3 id="([^"]+)">(.*?)</h3>', body, re.S)
        if not heading:
            continue
        section_id, title_html = heading.groups()
        title = _html_text(re.sub(r"<button.*?</button>", "", title_html, flags=re.S))
        ids = re.findall(r'\bid="([^"]+)"', body)
        notes = [_html_text(n) for n in re.findall(r'<div class="alert[^"]*">(.*?)</div>', body, re.S)]
        # Drop instructions for the interactive page ("Click a sample row ...").
        notes = [n for n in notes if not n.startswith("Click ")]
        sections.append({"id": section_id, "title": title, "plot_ids": ids, "notes": notes})
    return {"blocks": blocks, "sections": sections}


def read_multiqc_data(output_dir: Path, report_name: str = "multiqc_report") -> Dict[str, object]:
    """Everything in the MultiQC report, as data: what MultiQC writes to
    <output_dir>/<report_name>_data/multiqc_data.json, plus what only the
    HTML (<output_dir>/<report_name>.html) holds.

    Returns {"version", "tables", "general_stats_headers", "plots",
    "html_data", "sections"}:
      - "tables": every table MultiQC saved, as {table: {sample: {column:
        value}}} -- e.g. the full FastQC table (multiqc_fastqc) and General
        Statistics (multiqc_general_stats, with keys shortened to
        "<module>.<metric>").
      - "general_stats_headers": what each General Statistics column is
        ({key: {"title", "description", "suffix"}}).
      - "plots": the data behind every plot, as {plot_id: {"plot_type",
        "datasets", "samples" (bar graphs only), "config"}}, with only the
        config keys that label the data (_PLOT_CONFIG_KEYS) and no colors.
      - "html_data": data embedded only in the HTML, e.g. FastQC's per-base
        sequence content (see read_html_data).
      - "sections": every section of the report, in order, with the plot
        or data ids it shows and any note it shows instead (see
        read_html_data), so a section with nothing to plot (e.g. no adapter
        contamination found) is still recorded.
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

    general_stats_headers = {}
    for module_headers in data.get("report_general_stats_headers", []):
        for column, header in module_headers.items():
            # Same "<module>.<metric>" keys as the General Statistics table.
            key = f"{str(header.get('namespace', '')).lower()}.{column}" if header.get("namespace") else column
            general_stats_headers[key] = {
                k: header[k] for k in ("title", "description", "suffix") if header.get(k) not in (None, "")
            }

    plots = {}
    for plot_id, plot in data.get("report_plot_data", {}).items():
        config = plot.get("config", {})
        plots[plot_id] = _strip_colors({
            **{k: v for k, v in plot.items() if k != "config"},
            "config": {k: config[k] for k in _PLOT_CONFIG_KEYS if k in config},
        })

    html_data: Dict[str, object] = {}
    sections = []
    html_path = output_dir / f"{report_name}.html"
    if html_path.exists():
        from_html = read_html_data(html_path)
        html_data = from_html["blocks"]
        known = set(plots) | set(html_data)
        for section in from_html["sections"]:
            # Keep only ids that name data in these metrics (plots, html
            # blocks), not the page's layout ids.
            section["data"] = [i for i in section.pop("plot_ids") if i in known] + [
                block for block, owner in _HTML_BLOCK_SECTIONS.items()
                if owner == section["id"] and block in html_data
            ]
            sections.append(section)

    return {
        "version": data.get("config_version"),
        "tables": tables,
        "general_stats_headers": general_stats_headers,
        "plots": plots,
        "html_data": html_data,
        "sections": sections,
    }
