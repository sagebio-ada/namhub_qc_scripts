#!/usr/bin/env python
"""Run fq lint and kraken2/Bracken on the FASTQs of an nf-core/atacseq run.

Neither tool is part of nf-core/atacseq. Their output goes to --qc-dir, where atacseq_qc_metrics.py
reads it to fill the fq_lint_* and contam_* columns:
  <qc-dir>/fq_lint/{raw,trimmed}/<library>_T<n>.fq_lint.txt
  <qc-dir>/contaminants/kraken2/kraken_reports/<library>_T<n>.kraken2.report.txt
  <qc-dir>/contaminants/bracken/<library>_T<n>.tsv
Anything already there is skipped, so a rerun resumes. Kraken2 only runs if --kraken-db is given.

The sequencing runs (<library>_T<n>, one per samplesheet row) come from the samplesheet, as in
atacseq_qc_metrics.py. --samples limits the work to some libraries, so the libraries can be
processed in parallel (e.g. a Slurm array); --list-libraries prints the library names.
Needs `fq`, `kraken2` and `bracken` on PATH (the namhub_atac_qc env).
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from atacseq_common import fastqc_tables, read_samplesheet, select_runs  # noqa: E402
from qc_common import tools  # noqa: E402
from qc_common.util import to_num  # noqa: E402


def fastq_pair(run, r1, r2, trimmed_dir):
    """(raw FASTQs, trimmed FASTQs or None) of one sequencing run, run = <library>_T<n>."""
    raw = [r1] + ([r2] if r2 else [])
    trimmed = None
    if trimmed_dir:
        names = [f"{run}_1_val_1.fq.gz", f"{run}_2_val_2.fq.gz"] if r2 else [f"{run}_trimmed.fq.gz"]
        trimmed = [trimmed_dir / n for n in names]
        if not all(p.exists() for p in trimmed):
            print(f"WARNING [{run}] trimmed FASTQs not found in {trimmed_dir}", file=sys.stderr)
            trimmed = None
    return raw, trimmed


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--samplesheet", required=True, type=Path,
                   help="samplesheet CSV used to run nf-core/atacseq (sample, fastq_1, fastq_2, replicate)")
    p.add_argument("--list-libraries", action="store_true",
                   help="Print the library names, one per line, and exit")
    p.add_argument("--atac-dir", type=Path, help="nf-core/atacseq output directory (for the FastQC tables)")
    p.add_argument("--qc-dir", type=Path,
                   help="Where the output goes (default: <atac-dir>/qc_metrics)")
    p.add_argument("--peak-type", default="narrow_peak", help="narrow_peak (default) or broad_peak")
    p.add_argument("--trimmed-fastq-dir", type=Path,
                   help="Directory with the Trim Galore FASTQs (not saved by nf-core/atacseq by default)")
    p.add_argument("--kraken-db", help="Kraken2 database directory (with Bracken kmer_distrib files); "
                                       "without it kraken2/Bracken is not run")
    p.add_argument("--kraken-max-reads", type=int, default=1_000_000,
                   help="Screen only the first N reads (pairs) of each file; 0 = all (default: 1000000)")
    p.add_argument("--threads", type=int, default=8, help="kraken2 threads (default: 8)")
    p.add_argument("--samples", nargs="+", help="Only these libraries (default: all in the samplesheet)")
    args = p.parse_args()

    libs = read_samplesheet(args.samplesheet)
    if args.list_libraries:
        print("\n".join(libs))
        return
    if not args.atac_dir:
        sys.exit("--atac-dir is required")
    qc_dir = args.qc_dir or args.atac_dir / "qc_metrics"
    _, runs = select_runs(libs, args.samples)

    tools.need("fq")
    if args.kraken_db:
        tools.need("kraken2", "bracken")
    fq_tables = fastqc_tables(args.atac_dir / "multiqc" / args.peak_type / "multiqc_data")
    for s, (_, r1, r2) in runs.items():
        raw, trimmed = fastq_pair(s, r1, r2, args.trimmed_fastq_dir)
        tools.run_fq_lint(raw, qc_dir / "fq_lint" / "raw" / f"{s}.fq_lint.txt")
        if trimmed:
            tools.run_fq_lint(trimmed, qc_dir / "fq_lint" / "trimmed" / f"{s}.fq_lint.txt")
        if args.kraken_db:
            stage = "trimmed" if trimmed else "raw"
            row = fq_tables[stage].get(f"{s}_1") or fq_tables[stage].get(s)
            length = int(to_num(row["avg_sequence_length"])) if row else 100  # picks the Bracken k-mer file
            tools.run_kraken_bracken(trimmed or raw, s, qc_dir / "contaminants", args.kraken_db, length,
                                     threads=args.threads, max_reads=args.kraken_max_reads)
        print(f"[{s}] done", flush=True)


if __name__ == "__main__":
    main()
