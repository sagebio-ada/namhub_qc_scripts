"""Sections for FASTQ-based data (RNA-seq, and reusable for other assays).

Expected ctx keys:
  root        pipeline output directory (fq_lint/ lives here)
  fastqc      {"raw": {sample: row}, "trimmed": {...}} from the MultiQC FastQC tables
  contam_dir  directory holding kraken2/kraken_reports and bracken/
  species     expected species for the contamination check
  fastqc_id   optional, sample -> the id the FastQC tables use (default: the sample name itself);
              for libraries whose FastQC rows are named differently from the sample
"""
import re
import sys

from .util import pct, read_table, to_num

ANSI = re.compile(r"\x1b\[[0-9;]*m")
FASTQC_MODULES = [
    "basic_statistics", "per_base_sequence_quality", "per_sequence_quality_scores",
    "per_base_sequence_content", "per_sequence_gc_content", "per_base_n_content",
    "sequence_length_distribution", "sequence_duplication_levels",
    "overrepresented_sequences", "adapter_content",
]
STAGES = ["raw", "trimmed"]

FQ_LINT_KEYS = [f"fq_lint_{s}_{c}" for s in STAGES for c in ("status", "records", "issues")]


def fq_lint(ctx, sample):
    row = {}
    for stage in STAGES:
        path = ctx["root"] / "fq_lint" / stage / f"{sample}.fq_lint.txt"
        if not path.exists():  # e.g. no trimmed FASTQs saved: leave this stage blank
            print(f"WARNING [{sample}] fq_lint: missing {path}", file=sys.stderr)
            continue
        text = ANSI.sub("", path.read_text())
        issues = [ln for ln in text.splitlines() if re.search(r"\b(ERROR|WARN)\b", ln)]
        recs = re.findall(r"read (\d+)(?: \* \d+)? records", text)
        row[f"fq_lint_{stage}_status"] = (
            "fail" if issues else "pass" if "fq-lint end" in text else "incomplete")
        row[f"fq_lint_{stage}_records"] = int(recs[-1]) if recs else ""
        row[f"fq_lint_{stage}_issues"] = len(issues)
    return row


def fastqc_id(ctx, sample):
    return ctx.get("fastqc_id", str)(sample)


def fastqc_rows(table, sample):
    """FastQC rows for read 1 (and read 2). Paired rows are '<sample>_1'/'_2'; single-end '<sample>'."""
    rows = [table[f"{sample}_{r}"] for r in (1, 2) if f"{sample}_{r}" in table]
    if not rows and sample in table:
        rows = [table[sample]]
    return rows


READ_LENGTH_KEYS = ["read_length_r1", "read_length_r2",
                    "read_length_trimmed_r1", "read_length_trimmed_r2"]


def read_length(ctx, sample):
    row = {}
    for stage, prefix in (("raw", "read_length"), ("trimmed", "read_length_trimmed")):
        rows = fastqc_rows(ctx["fastqc"][stage], fastqc_id(ctx, sample))
        lens = [to_num(r["median_sequence_length"]) for r in rows] + ["", ""]
        row[f"{prefix}_r1"], row[f"{prefix}_r2"] = lens[0], lens[1]
    return row


CONTAM_KEYS = ["contam_unclassified_pct", "contam_top5_species", "contam_expected_species",
               "contam_expected_pct_of_classified", "contam_expected_pct_of_all_reads"]


def contaminants(ctx, sample):
    cdir = ctx["contam_dir"]
    unclassified = classified = None
    with open(cdir / "kraken2" / "kraken_reports" / f"{sample}.kraken2.report.txt") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[3] == "U":
                unclassified = int(f[1])
            elif f[3] == "R" and f[4] == "1":
                classified = int(f[1])
    if unclassified is None or classified is None:
        raise ValueError("no unclassified/root line in the kraken2 report")
    total = unclassified + classified

    species = [(r["name"], int(r["new_est_reads"])) for r in
               read_table(cdir / "bracken" / f"{sample}.tsv") if r["taxonomy_lvl"] == "S"]
    species.sort(key=lambda x: -x[1])
    est_total = sum(n for _, n in species)
    top5 = "; ".join(f"{name}:{n}({pct(n, est_total):.2f}%)" for name, n in species[:5])
    expected = dict((n.lower(), c) for n, c in species).get(ctx["species"].lower(), 0)
    return {"contam_unclassified_pct": pct(unclassified, total), "contam_top5_species": top5,
            "contam_expected_species": ctx["species"],
            "contam_expected_pct_of_classified": pct(expected, est_total),
            "contam_expected_pct_of_all_reads": pct(expected, total)}


FASTQC_KEYS = (["raw_reads", "trimmed_reads", "trimmed_pct_retained"]
               + [f"fastqc_{s}_{m}_{r}" for s in STAGES
                  for m in ("gc_pct", "dedup_pct", "poor_quality_seqs") for r in ("r1", "r2")]
               + [f"fastqc_{s}_{m}" for s in STAGES for m in FASTQC_MODULES])


def fastqc(ctx, sample):
    row, totals = {}, {}
    for stage in STAGES:
        rows = fastqc_rows(ctx["fastqc"][stage], fastqc_id(ctx, sample))
        if not rows:
            raise KeyError(f"no {stage} FastQC rows for {sample}")
        totals[stage] = int(to_num(rows[0]["Total Sequences"]))
        for name, col in (("gc_pct", "%GC"), ("dedup_pct", "total_deduplicated_percentage"),
                          ("poor_quality_seqs", "Sequences flagged as poor quality")):
            vals = [to_num(r[col]) for r in rows] + ["", ""]
            row[f"fastqc_{stage}_{name}_r1"], row[f"fastqc_{stage}_{name}_r2"] = vals[0], vals[1]
        for m in FASTQC_MODULES:
            row[f"fastqc_{stage}_{m}"] = "/".join(r[m] for r in rows)
    row.update({"raw_reads": totals["raw"], "trimmed_reads": totals["trimmed"],
                "trimmed_pct_retained": pct(totals["trimmed"], totals["raw"])})
    return row


