#!/usr/bin/env python
"""Summarize QC metrics from an nf-core/atacseq run (bwa, narrow peaks) into one TSV row per library.

Input: the nf-core/atacseq output directory (--atac-dir), the samplesheet used to run it
(--samplesheet), and the table written by atacseq_encode_qc_metrics.py (--encode-qc), and
optionally the table(s) written by atacseq_idr_concordance.py (--idr-metrics). A library is
<sample>_REP<replicate> (the samplesheet's sample and replicate columns), as in nf-core's file
names. Samplesheet rows of one library are technical replicates (sequencing runs): nf-core names
them <library>_T1, <library>_T2, ... in sheet order and merges them after alignment.
A missing or unparsable file leaves that section's columns blank and prints a warning.

Output: one table with two kinds of rows, joined on the `library` column.
  run rows      sample = <library>_T<n>, row_type = run. Pre-merge sections only (fq lint, read
                length, contamination, FastQC); the post-merge columns are blank.
  library rows  sample = library = <library>, row_type = library. Post-merge sections only
                (ENCODE table, FRiP, TSS enrichment); the per-FASTQ columns are blank.
Nothing is aggregated across runs. Each library has its run rows followed by its library row.

fq lint and kraken2/Bracken are not part of nf-core/atacseq. This script does not run them: run
atacseq_tools_qc.py first with the same --qc-dir (default <atac-dir>/qc_metrics). Where their
output is missing, the fq_lint_* / contam_* columns are blank and a warning is printed.

Columns (paths are relative to --atac-dir)
  sample, library, row_type
      Row id, the library it belongs to (join key), and "run" or "library".
  fq_lint_{raw,trimmed}_{status,records,issues}
      As in rnaseq_qc_metrics.py. Raw = the samplesheet FASTQs. nf-core/atacseq does not keep
      trimmed FASTQs unless asked to, so the trimmed columns are blank unless --trimmed-fastq-dir
      was given to atacseq_tools_qc.py (Trim Galore names: <library>_T<n>_1_val_1.fq.gz, ..._2_val_2.fq.gz).
  read_length_r1, read_length_r2, read_length_trimmed_r1, read_length_trimmed_r2
      Median read length. This MultiQC has no median column, so it is computed from the
      FastQC zips (fastqc/zips, trimgalore/fastqc/zips), as the midpoint-weighted median of FastQC's
      Sequence Length Distribution.
  contam_*  As in rnaseq_qc_metrics.py (kraken2 + Bracken at species level). Which FASTQs were
      screened, and how many reads, is set when running atacseq_tools_qc.py.
  raw_reads, trimmed_reads, trimmed_pct_retained, fastqc_*
      As in rnaseq_qc_metrics.py, from multiqc/<peak-type>/multiqc_data/multiqc_fastqc*.txt.
  alignment_rate_raw, alignment_rate_nonmt, nonmt_nondup_aligned_reads, NRF, PBC1, PBC2,
  fragment_size_median, nfr_fraction, mono_fraction, nfr_over_mono, nfr_peak_found, nfr_peak_bp,
  mono_peak_found, mono_peak_bp
      Copied from the atacseq_encode_qc_metrics.py table (--encode-qc).
  frip
      bwa/merged_library/macs2/<peak-type>/qc/<library>.mLb.clN_peaks.FRiP_mqc.tsv
  ataqv_hqaa_overlapping_peaks_pct, tss_enrichment
      hqaa_overlapping_peaks_percent and tss_enrichment from
      bwa/merged_library/ataqv/<peak-type>/<library>.ataqv.json.
  idr_Nt, idr_Np, idr_N_optimal, idr_rescue_ratio, idr_self_consistency_ratio,
  idr_reproducibility, overlap_N_optimal
      Copied from the atacseq_idr_concordance.py table(s) (--idr-metrics), columns Nt, Np, N_optimal,
      rescue_ratio, self_consistency_ratio, reproducibility and overlap_N_optimal (IDR peaks after
      the blacklist/chromosome filter). That table has one row per group (the sample, i.e. the
      library name without _REP<n>), so every replicate library of a group gets the same values.
      Blank if --idr-metrics is not given.
"""
import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from atacseq_common import fastqc_tables, read_samplesheet, select_runs  # noqa: E402
from qc_common.fastq_sections import (  # noqa: E402
    CONTAM_KEYS, FASTQC_KEYS, FQ_LINT_KEYS, contaminants, fastqc, fq_lint)
from qc_common.fastqc_zip import median_read_length  # noqa: E402
from qc_common.runner import collect, write_tsv  # noqa: E402
from qc_common.util import read_table, to_num  # noqa: E402

ENCODE_KEYS = ["alignment_rate_raw", "alignment_rate_nonmt", "nonmt_nondup_aligned_reads",
               "NRF", "PBC1", "PBC2", "fragment_size_median", "nfr_fraction", "mono_fraction",
               "nfr_over_mono", "nfr_peak_found", "nfr_peak_bp", "mono_peak_found", "mono_peak_bp"]
READ_LENGTH_KEYS = ["read_length_r1", "read_length_r2",
                    "read_length_trimmed_r1", "read_length_trimmed_r2"]
FRIP_KEYS = ["frip", "ataqv_hqaa_overlapping_peaks_pct"]
ATAQV_KEYS = ["tss_enrichment"]
# column in the IDR table -> column here
IDR_COLUMNS = {"Nt": "idr_Nt", "Np": "idr_Np", "N_optimal": "idr_N_optimal",
               "rescue_ratio": "idr_rescue_ratio",
               "self_consistency_ratio": "idr_self_consistency_ratio",
               "reproducibility": "idr_reproducibility", "overlap_N_optimal": "overlap_N_optimal"}
IDR_KEYS = list(IDR_COLUMNS.values())


def read_length(ctx, sample):
    fid, row = sample, {}
    paired = ctx["paired"][sample]
    for stage, prefix in (("raw", "read_length"), ("trimmed", "read_length_trimmed")):
        zips = ctx["zips"][stage]
        lens = []
        for r in ((1, 2) if paired else (None,)):
            z = zips(fid, r)
            lens.append(median_read_length(z) if z.exists() else "")
            if not z.exists():
                print(f"WARNING [{sample}] read length: missing {z}", file=sys.stderr)
        lens += ["", ""]
        row[f"{prefix}_r1"], row[f"{prefix}_r2"] = lens[0], lens[1]
    return row


def encode(ctx, sample):
    return {k: ctx["encode"][sample][k] for k in ENCODE_KEYS}


def ataqv(ctx, sample):
    path = ctx["merged_lib"] / "ataqv" / ctx["peak_type"] / f"{sample}.ataqv.json"
    with open(path) as fh:
        return json.load(fh)[0]["metrics"]


def frip(ctx, sample):
    path = (ctx["merged_lib"] / "macs2" / ctx["peak_type"] / "qc" /
            f"{sample}.mLb.clN_peaks.FRiP_mqc.tsv")
    lines = [ln.split("\t") for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    return {"frip": to_num(next(f[1] for f in lines if f[0].strip("'\"") == sample)),
            "ataqv_hqaa_overlapping_peaks_pct": ataqv(ctx, sample)["hqaa_overlapping_peaks_percent"]}


def tss(ctx, sample):
    return {"tss_enrichment": ataqv(ctx, sample)["tss_enrichment"]}


def idr(ctx, sample):
    if not ctx["idr"]:  # --idr-metrics not given
        return {}
    group = re.sub(r"_REP\d+$", "", sample)
    if group not in ctx["idr"]:
        raise KeyError(f"no row for group {group} in --idr-metrics (IDR needs >= 2 replicates)")
    return {new: ctx["idr"][group][old] for old, new in IDR_COLUMNS.items()}


# Pre-merge sections run on each sequencing run (<library>_T<n>), post-merge ones on the library.
RUN_SECTIONS = [
    ("fq_lint", FQ_LINT_KEYS, fq_lint),
    ("read length", READ_LENGTH_KEYS, read_length),
    ("contaminants", CONTAM_KEYS, contaminants),
    ("FastQC", FASTQC_KEYS, fastqc),
]
LIBRARY_SECTIONS = [
    ("ENCODE QC", ENCODE_KEYS, encode),
    ("FRiP", FRIP_KEYS, frip),
    ("TSS enrichment", ATAQV_KEYS, tss),
    ("IDR", IDR_KEYS, idr),
]
ID_KEYS = ["library", "row_type"]  # written after "sample"; only used for the header


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--atac-dir", required=True, type=Path, help="nf-core/atacseq output directory")
    p.add_argument("--samplesheet", required=True, type=Path,
                   help="samplesheet CSV used to run nf-core/atacseq (sample, fastq_1, fastq_2, replicate)")
    p.add_argument("--encode-qc", required=True, type=Path,
                   help="TSV written by atacseq_encode_qc_metrics.py")
    p.add_argument("--idr-metrics", nargs="+", type=Path,
                   help="concordance_metrics.tsv file(s) from atacseq_idr_concordance.py (one per "
                        "run of that script); adds the IDR / overlap columns. Optional")
    p.add_argument("--out", required=True, type=Path, help="Output TSV")
    p.add_argument("--qc-dir", type=Path,
                   help="Where atacseq_tools_qc.py wrote the fq lint and kraken2/Bracken output "
                        "(default: <atac-dir>/qc_metrics)")
    p.add_argument("--peak-type", default="narrow_peak", help="narrow_peak (default) or broad_peak")
    p.add_argument("--species", default="Homo sapiens",
                   help="Expected species for the contamination check (default: Homo sapiens)")
    p.add_argument("--samples", nargs="+", help="Only these libraries (default: all in the samplesheet)")
    args = p.parse_args()

    root = args.atac_dir
    qc_dir = args.qc_dir or root / "qc_metrics"
    libs = read_samplesheet(args.samplesheet)
    samples, runs = select_runs(libs, args.samples)
    fq_tables = fastqc_tables(root / "multiqc" / args.peak_type / "multiqc_data")

    encode_rows = {r["sample"]: r for r in read_table(args.encode_qc)}
    ctx = {"root": qc_dir, "contam_dir": qc_dir / "contaminants", "species": args.species,
           "fastqc": fq_tables,
           "merged_lib": root / "bwa" / "merged_library", "peak_type": args.peak_type,
           "encode": encode_rows,
           "idr": {r["group"]: r for f in args.idr_metrics or [] for r in read_table(f)},
           "paired": {s: r2 is not None for s, (_, _, r2) in runs.items()},
           "zips": {"raw": lambda i, r: root / "fastqc" / "zips" / (f"{i}_{r}_fastqc.zip" if r else f"{i}_fastqc.zip"),
                    "trimmed": lambda i, r: root / "trimgalore" / "fastqc" / "zips" /
                    (f"{i}_{r}_val_{r}_fastqc.zip" if r else f"{i}_trimmed_fastqc.zip")}}

    rows = []
    for lib in samples:
        for run in (r for r, v in runs.items() if v[0] == lib):
            rows.append({**collect(run, ctx, RUN_SECTIONS), "library": lib, "row_type": "run"})
        rows.append({**collect(lib, ctx, LIBRARY_SECTIONS), "library": lib, "row_type": "library"})
    fieldnames = write_tsv(rows, [("ids", ID_KEYS, None)] + RUN_SECTIONS + LIBRARY_SECTIONS, args.out)
    print(f"Wrote {len(rows)} rows ({len(runs)} runs, {len(samples)} libraries) x "
          f"{len(fieldnames) - 1} columns to {args.out}")


if __name__ == "__main__":
    main()
