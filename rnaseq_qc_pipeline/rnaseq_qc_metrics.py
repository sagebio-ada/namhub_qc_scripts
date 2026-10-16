#!/usr/bin/env python
"""Summarize QC metrics from an nf-core/rnaseq run (aligner star_salmon) into one TSV row per sample.

Input: the nf-core/rnaseq output directory (--rnaseq-dir). Samples are the unique values of
the "sample" column of the nf-core samplesheet (--samplesheet, CSV, in order of first
appearance; replicate rows share a sample). Without --samplesheet they are the ones with a
<outdir>/star_salmon/log/<sample>.Log.final.out (STAR only). Needs numpy and scipy (replicate
correlation); everything else uses the standard library.

NOTE: only tested with aligner star_salmon. Sample discovery no longer depends on STAR when
--samplesheet is given, but the STAR, Qualimap, featureCounts, Picard and Salmon sections still
read star_salmon-style paths under <outdir>/<aligner>/ and may be blank (with warnings) or need
fixes for other aligners. Not yet tested.
A missing or unparsable file leaves that section's columns blank and prints a warning.

Columns (paths are relative to --rnaseq-dir)
  fq_lint_{raw,trimmed}_{status,records,issues}
      fq_lint/{raw,trimmed}/<sample>.fq_lint.txt. status = pass if the log has no ERROR/WARN
      lines and ends with "fq-lint end", fail if it has any, incomplete if it never finished.
      records = read pairs (paired-end) or reads (single-end) as fq-lint counted them.
  read_length_r1, read_length_r2, read_length_trimmed_r1, read_length_trimmed_r2
      Median read length from FastQC, raw (multiqc_fastqc_fastqc_raw.txt) and trimmed
      (multiqc_fastqc_fastqc_trimmed.txt). Each pair is blank if its table has no row for the
      sample; r2 is blank for single-end data.
  strand_provided, strand_salmon, strand_rseqc, strand_inferred_agree, strand_status
      multiqc_strand_check_summary_table.txt. agree = Salmon and RSeQC inferred the same
      strandedness and it is not "undetermined". status is the pipeline's own value.
  contam_unclassified_pct, contam_top5_species, contam_expected_species,
  contam_expected_pct_of_classified, contam_expected_pct_of_all_reads
      kraken2 report (unclassified) and Bracken TSV (species). Top 5 are species-level by
      Bracken-estimated reads, as name:reads(% of classified reads). The expected species
      (--species, default Homo sapiens) is given as % of classified reads and % of all
      reads (classified + unclassified).
  raw_reads, trimmed_reads, trimmed_pct_retained
      FastQC "Total Sequences" of read 1 (= read pairs for paired-end).
  fastqc_{raw,trimmed}_{gc_pct,dedup_pct,poor_quality_seqs}_{r1,r2}
      FastQC %GC, % remaining after deduplication, and sequences flagged as poor quality.
  fastqc_{raw,trimmed}_<module>
      FastQC pass/warn/fail for each module, "r1/r2" for paired-end.
  salmon_fraglen_mean, salmon_fraglen_sd
      star_salmon/<sample>/logs/salmon_quant.log, "fragment length mean (sd)".
  star_input_reads, star_total_aligned, star_total_aligned_pct
      STAR Log.final.out. Total aligned (ENCODE definition) = uniquely mapped + mapped to
      multiple loci (STAR --outFilterMultimapNmax 20 in nf-core/rnaseq).
  star_unique_{n,pct}, star_multi_{n,pct}, star_too_many_loci_{n,pct}, star_unmapped_{n,pct}
      multi = "mapped to multiple loci" only; too many loci is reported separately; unmapped =
      too many mismatches + too short + other. Percentages are summed from STAR's rounded values.
  qualimap_exonic_pct, qualimap_intronic_pct, qualimap_intergenic_pct,
  qualimap_read_pairs_aligned
      star_salmon/qualimap/<sample>/rnaseq_qc_results.txt.
  featurecounts_assigned, featurecounts_no_features, featurecounts_ambiguous,
  featurecounts_multimapping, featurecounts_assigned_pct
      star_salmon/featurecounts/<sample>.featureCounts.tsv.summary. nf-core runs featureCounts
      with -p but not --countReadPairs, so for paired-end data each mate is counted (about 2x
      the read pairs); multimapping alignments are not counted as assigned. assigned_pct =
      assigned / all counted alignments.
  rrna_pct_nfcore, rrna_reads, mt_rrna_reads, rrna_plus_mt_pct
      nf-core's % rRNA (biotype_counts_rrna_mqc.tsv; this is a percentage, and only the exact
      "rRNA" biotype), and read counts of the rRNA and Mt_rRNA biotypes from
      biotype_counts_mqc.tsv. rrna_plus_mt_pct = (rRNA + Mt_rRNA) / all biotype-assigned reads.
  picard_duplication_rate, picard_estimated_library_size
      star_salmon/picard_metrics/<sample>.markdup.sorted.metrics.txt: PERCENT_DUPLICATION
      (a fraction, 0-1) and ESTIMATED_LIBRARY_SIZE (blank for single-end data).
  spearman_isogenic, spearman_anisogenic
      Only with --replicate-sheet. Spearman correlation of gene-level TPM
      (star_salmon/salmon.merged.gene_tpm.tsv) between this sample and the other samples of its
      group, as the median over those partners, separately for isogenic and anisogenic pairs
      (blank if the sample has no such partner). Genes are chosen as in ENCODE's MAD.R: both TPM
      > 0 and mean(log2 TPM) of the pair > 0, per pair. ENCODE uses RSEM FPKM; this uses Salmon
      TPM, so the cutoff is close to ENCODE's but not identical. Reported only, no pass/fail.
      Every pair is also written to <out>_replicate_pairs.tsv (sample_a, sample_b, group,
      relationship, n_genes, spearman).
      The replicate sheet is a CSV with columns sample, group, donor (and optionally replicate, a
      label that is not used). Samples in the same group are compared pairwise: same donor
      (genetic source: donor, individual or cell line) = isogenic, different donors = anisogenic.
"""
import argparse
import csv
import re
import statistics
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qc_common.fastq_sections import (  # noqa: E402
    CONTAM_KEYS, FASTQC_KEYS, FQ_LINT_KEYS, READ_LENGTH_KEYS, contaminants, fastqc, fq_lint,
    read_length)
from qc_common.runner import collect, write_tsv  # noqa: E402
from qc_common.samplesheet import read_samples  # noqa: E402
from qc_common.util import pct, read_table, to_num  # noqa: E402


def star_fields(path):
    """STAR Log.final.out: 'label | value' lines -> {label: value}."""
    out = {}
    with open(path) as fh:
        for ln in fh:
            if "|" in ln:
                k, v = ln.split("|", 1)
                out[k.strip()] = v.strip()
    return out


# ---------------------------------------------------------------- sections
# Each section function takes (ctx, sample) and returns a dict with exactly its KEYS.
# fq_lint, read length, FastQC and contaminants come from qc_common.fastq_sections.

STRAND_KEYS = ["strand_provided", "strand_salmon", "strand_rseqc", "strand_inferred_agree",
               "strand_status"]


def strandedness(ctx, sample):
    r = ctx["strand"][sample]
    agree = r["salmon_inferred"] == r["rseqc_inferred"] and r["salmon_inferred"] != "undetermined"
    return {"strand_provided": r["provided"], "strand_salmon": r["salmon_inferred"],
            "strand_rseqc": r["rseqc_inferred"], "strand_inferred_agree": agree,
            "strand_status": r["status"]}


FRAGLEN_KEYS = ["salmon_fraglen_mean", "salmon_fraglen_sd"]


def fragment_length(ctx, sample):
    text = (ctx["aligner_dir"] / sample / "logs" / "salmon_quant.log").read_text()
    m = re.search(r"fragment length mean \(sd\):\s*([\d.]+)\s*\(([\d.]+)\)", text)
    if not m:
        raise ValueError("no 'fragment length mean (sd)' line")
    return {"salmon_fraglen_mean": float(m[1]), "salmon_fraglen_sd": float(m[2])}


STAR_KEYS = (["star_input_reads", "star_total_aligned", "star_total_aligned_pct"]
             + [f"star_{c}_{u}" for c in ("unique", "multi", "too_many_loci", "unmapped")
                for u in ("n", "pct")])


def star(ctx, sample):
    f = star_fields(ctx["aligner_dir"] / "log" / f"{sample}.Log.final.out")
    n = lambda k: to_num(f[k])  # noqa: E731
    unmapped = [("Number of reads unmapped: " + w, "% of reads unmapped: " + w)
                for w in ("too many mismatches", "too short", "other")]
    uniq, multi = n("Uniquely mapped reads number"), n("Number of reads mapped to multiple loci")
    uniq_p, multi_p = n("Uniquely mapped reads %"), n("% of reads mapped to multiple loci")
    return {"star_input_reads": n("Number of input reads"),
            "star_total_aligned": uniq + multi, "star_total_aligned_pct": round(uniq_p + multi_p, 2),
            "star_unique_n": uniq, "star_unique_pct": uniq_p,
            "star_multi_n": multi, "star_multi_pct": multi_p,
            "star_too_many_loci_n": n("Number of reads mapped to too many loci"),
            "star_too_many_loci_pct": n("% of reads mapped to too many loci"),
            "star_unmapped_n": sum(n(a) for a, _ in unmapped),
            "star_unmapped_pct": round(sum(n(b) for _, b in unmapped), 2)}


QUALIMAP_KEYS = ["qualimap_exonic_pct", "qualimap_intronic_pct", "qualimap_intergenic_pct",
                 "qualimap_read_pairs_aligned"]


def qualimap(ctx, sample):
    text = (ctx["aligner_dir"] / "qualimap" / sample / "rnaseq_qc_results.txt").read_text()
    row = {}
    for key, label in (("exonic", "exonic"), ("intronic", "intronic"), ("intergenic", "intergenic")):
        m = re.search(rf"^\s*{label}\s*=\s*[\d,]+\s*\(([\d.]+)%\)", text, re.M)
        row[f"qualimap_{key}_pct"] = float(m[1]) if m else ""
    m = re.search(r"^\s*read pairs aligned\s*=\s*([\d,]+)", text, re.M)
    row["qualimap_read_pairs_aligned"] = to_num(m[1]) if m else ""
    return row


FEATURECOUNTS_KEYS = ["featurecounts_assigned", "featurecounts_no_features",
                      "featurecounts_ambiguous", "featurecounts_multimapping",
                      "featurecounts_assigned_pct"]


def featurecounts(ctx, sample):
    path = ctx["aligner_dir"] / "featurecounts" / f"{sample}.featureCounts.tsv.summary"
    counts = {r[0]: int(r[1]) for r in (ln.rstrip("\n").split("\t") for ln in open(path).readlines()[1:])}
    return {"featurecounts_assigned": counts["Assigned"],
            "featurecounts_no_features": counts["Unassigned_NoFeatures"],
            "featurecounts_ambiguous": counts["Unassigned_Ambiguity"],
            "featurecounts_multimapping": counts["Unassigned_MultiMapping"],
            "featurecounts_assigned_pct": pct(counts["Assigned"], sum(counts.values()))}


RRNA_KEYS = ["rrna_pct_nfcore", "rrna_reads", "mt_rrna_reads", "rrna_plus_mt_pct"]


def rrna(ctx, sample):
    fdir = ctx["aligner_dir"] / "featurecounts"
    pct_nf = [float(r["percent_rRNA"]) for r in read_table(fdir / f"{sample}.biotype_counts_rrna_mqc.tsv")
              if r["Sample"].strip("'\"") == sample]
    with open(fdir / f"{sample}.biotype_counts_mqc.tsv") as fh:
        biotypes = {f[0]: int(float(f[1])) for f in
                    (ln.rstrip("\n").split("\t") for ln in fh if ln.strip() and not ln.startswith("#"))}
    r, mt = biotypes.get("rRNA", 0), biotypes.get("Mt_rRNA", 0)
    return {"rrna_pct_nfcore": pct_nf[0] if pct_nf else "", "rrna_reads": r, "mt_rrna_reads": mt,
            "rrna_plus_mt_pct": pct(r + mt, sum(biotypes.values()))}


PICARD_KEYS = ["picard_duplication_rate", "picard_estimated_library_size"]


def picard(ctx, sample):
    path = ctx["aligner_dir"] / "picard_metrics" / f"{sample}.markdup.sorted.metrics.txt"
    lines = path.read_text().splitlines()
    i = next(i for i, ln in enumerate(lines) if ln.startswith("LIBRARY\t"))
    r = dict(zip(lines[i].split("\t"), lines[i + 1].split("\t")))
    return {"picard_duplication_rate": to_num(r["PERCENT_DUPLICATION"]),
            "picard_estimated_library_size": to_num(r.get("ESTIMATED_LIBRARY_SIZE", ""))}


REPCORR_KEYS = ["spearman_isogenic", "spearman_anisogenic"]
RELATIONSHIPS = ("isogenic", "anisogenic")


def read_replicate_sheet(path):
    """CSV (sample, group, donor) -> {sample: (group, donor)}."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = {"sample", "group", "donor"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"{path}: missing column(s) {', '.join(sorted(missing))}")
        sheet = {}
        for r in reader:
            group, donor = r["group"].strip(), r["donor"].strip()
            if not donor:
                sys.exit(f"{path}: no donor for {r['sample']}")
            sheet[r["sample"].strip()] = (group, donor)
    return sheet


def read_gene_tpm(path):
    """salmon.merged.gene_tpm.tsv -> {sample: TPM array}, in file (gene) order."""
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        cols = list(zip(*reader))
    return {name: np.array(col, dtype=float) for name, col in zip(header[2:], cols[2:])}


def spearman_pair(x, y):
    """Spearman on the genes ENCODE's MAD.R keeps: expressed in both and mean log2 > 0."""
    keep = (x > 0) & (y > 0)
    keep[keep] = (np.log2(x[keep]) + np.log2(y[keep])) / 2 > 0
    if keep.sum() < 3:
        return int(keep.sum()), ""
    return int(keep.sum()), float(spearmanr(x[keep], y[keep])[0])


def replicate_correlations(tpm, sheet):
    """Pairs of samples in the same group -> list of row dicts.

    relationship: isogenic if both samples have the same donor, otherwise anisogenic.
    """
    rows = []
    for a, b in combinations([s for s in sheet if s in tpm], 2):
        (ga, da), (gb, db) = sheet[a], sheet[b]
        if ga != gb:
            continue
        n, rho = spearman_pair(tpm[a], tpm[b])
        rows.append({"sample_a": a, "sample_b": b, "group": ga,
                     "relationship": "isogenic" if da == db else "anisogenic",
                     "n_genes": n, "spearman": rho})
    return rows


def replicate_correlation(ctx, sample):
    row = {}
    for kind in RELATIONSHIPS:
        vals = [r["spearman"] for r in ctx["rep_pairs"]
                if sample in (r["sample_a"], r["sample_b"]) and r["relationship"] == kind
                and r["spearman"] != ""]
        row[f"spearman_{kind}"] = statistics.median(vals) if vals else ""
    return row


SECTIONS = [
    ("fq_lint", FQ_LINT_KEYS, fq_lint),
    ("read length", READ_LENGTH_KEYS, read_length),
    ("strandedness", STRAND_KEYS, strandedness),
    ("contaminants", CONTAM_KEYS, contaminants),
    ("FastQC", FASTQC_KEYS, fastqc),
    ("salmon fragment length", FRAGLEN_KEYS, fragment_length),
    ("STAR", STAR_KEYS, star),
    ("Qualimap", QUALIMAP_KEYS, qualimap),
    ("featureCounts", FEATURECOUNTS_KEYS, featurecounts),
    ("rRNA", RRNA_KEYS, rrna),
    ("Picard", PICARD_KEYS, picard),
    ("replicate correlation", REPCORR_KEYS, replicate_correlation),
]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rnaseq-dir", required=True, type=Path,
                   help="nf-core/rnaseq output directory (contains fq_lint/, star_salmon/, multiqc/)")
    p.add_argument("--out", required=True, type=Path, help="Output TSV")
    p.add_argument("--samplesheet", type=Path,
                   help="nf-core/rnaseq samplesheet CSV (with a 'sample' column) used to run the "
                        "pipeline; sample names are taken from it (default: STAR Log.final.out files)")
    p.add_argument("--aligner", default="star_salmon",
                   help="Aligner/quantifier subdirectory (default: star_salmon)")
    p.add_argument("--species", default="Homo sapiens",
                   help="Expected species for the contamination check (default: Homo sapiens)")
    p.add_argument("--replicate-sheet", type=Path,
                   help="CSV with columns sample, group, donor (replicate optional); "
                        "enables spearman_isogenic / spearman_anisogenic")
    p.add_argument("--samples", nargs="+",
                   help="Only these samples (default: all in the samplesheet, or all with a STAR Log.final.out)")
    args = p.parse_args()

    root = args.rnaseq_dir
    aligner_dir = root / args.aligner
    mqc = root / "multiqc" / args.aligner / "multiqc_report_data"

    if args.samplesheet:
        samples = read_samples(args.samplesheet)
    else:
        samples = sorted(f.name[: -len(".Log.final.out")]
                         for f in (aligner_dir / "log").glob("*.Log.final.out"))
        if not samples:
            sys.exit(f"No *.Log.final.out found in {aligner_dir / 'log'}; "
                     "use --samplesheet for a non-STAR run")
    if args.samples:
        unknown = set(args.samples) - set(samples)
        if unknown:
            sys.exit(f"Samples not found: {', '.join(sorted(unknown))}")
        samples = args.samples

    def table(name, key):
        path = mqc / name
        if not path.exists():
            print(f"WARNING: {path} not found", file=sys.stderr)
            return {}
        return {r[key]: r for r in read_table(path)}

    ctx = {"root": root, "aligner_dir": aligner_dir, "contam_dir": aligner_dir / "contaminants",
           "species": args.species,
           "fastqc": {"raw": table("multiqc_fastqc_fastqc_raw.txt", "Sample"),
                      "trimmed": table("multiqc_fastqc_fastqc_trimmed.txt", "Sample")},
           "strand": table("multiqc_strand_check_summary_table.txt", "Sample")}

    ctx["rep_pairs"] = []
    if args.replicate_sheet:
        tpm_path = aligner_dir / "salmon.merged.gene_tpm.tsv"
        if tpm_path.exists():
            ctx["rep_pairs"] = replicate_correlations(read_gene_tpm(tpm_path),
                                                      read_replicate_sheet(args.replicate_sheet))
            pairs_path = args.out.with_name(args.out.stem + "_replicate_pairs.tsv")
            with open(pairs_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
                    "sample_a", "sample_b", "group", "relationship", "n_genes", "spearman"])
                w.writeheader()
                w.writerows(ctx["rep_pairs"])
        else:
            print(f"WARNING: {tpm_path} not found; replicate correlations left blank", file=sys.stderr)

    rows = [collect(s, ctx, SECTIONS) for s in samples]
    fieldnames = write_tsv(rows, SECTIONS, args.out)
    print(f"Wrote {len(rows)} samples x {len(fieldnames) - 1} metrics to {args.out}")


if __name__ == "__main__":
    main()
