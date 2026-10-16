#!/usr/bin/env python
"""Compute ATAC-seq QC metrics from nf-core/atacseq duplicate-marked BAMs.

Input: <outdir>/bwa/merged_library/*.mLb.mkD.sorted.bam (coordinate-sorted,
duplicates flagged but not removed, unmapped reads retained). These are only
published when the pipeline is run with `save_align_intermeds: true`.

Requires `samtools` and `awk` on PATH (e.g. `module load samtools`).

Metrics per sample
  Alignment rate (raw BAM, all contigs incl. MT, before dedup/filtering; approximates the
    mapped % ENCODE takes from SAMstats on its raw bowtie2 BAM)
  Non-MT: total / unique / multi-mapped / unmapped primary reads (+ proportions)
  Mito fraction (MT mapped reads / all mapped reads, raw BAM). Not ENCODE's method: ENCODE
    counts MT reads (Rm) by aligning all reads again to a mito-only index and reports
    Rm / (Rn + Rm), Rn = non-MT mapped reads of the main alignment
  ENCODE-style duplication rate (flag 0x400 among reads passing
    -F 1804(minus dup bit) -f 2 -q 30, all contigs incl. MT). Approximate: nf-core
    marked duplicates on the unfiltered BAM, ENCODE marks them on the filtered one.
  Non-MT duplication rate (flag 0x400 among mapped primary reads)
  Non-MT non-duplicate ENCODE-"aligned" reads (-F 1804 -f 2 -q 30)
  Non-MT NRF, PBC1, PBC2 (ENCODE definitions, on fragment positions)
  Fragment-size metrics following the ENCODE ATAC-seq pipeline (encode_task_fraglen_stat_pe.py:
    same cutoffs, windows and peak finder), from the Picard insert-size histogram of nf-core's final
    filtered BAM (<bam-dir>/picard_metrics/<sample>.mLb.clN.CollectMultipleMetrics.insert_size_metrics):
    fraction of read pairs in the NFR (< 150 bp), NFR / mono-nucleosome read ratio (mono =
    > 150 and <= 300 bp), mono fraction, and whether a NFR (20-90 bp), mono-nucleosome (120-250 bp)
    and di-nucleosome (300-500 bp) peak was found with scipy find_peaks_cwt (width 25),
    with its fragment length (the tallest peak if several fall in a window).
    Reported only; no pass/fail. Counts are divided by READ_PAIRS
    from the Picard summary because the histogram itself is trimmed (median + 10 x MAD).
    Differences from ENCODE: ENCODE runs Picard on the first 5,000,000 reads only
    (STOP_AFTER, histogram width 1000) of its own filtered BAM (bowtie2, MAPQ >= 30), and divides
    by the histogram sum; here the whole nf-core clN BAM is used (the histogram sum and READ_PAIRS
    differed by 0.01% in the GSE195688 test).
    The windows and peak-finder width are the ENCODE_* constants below, or command-line options.

NRF/PBC differ slightly from ENCODE's `bedtools bamtobed -bedpe` method: each pair
is counted once via its leftmost mate (TLEN > 0) after the per-read -q filter, with
no fixmate/orphan removal. ENCODE's fragment key (R1 start, R2 end, both strands) also
separates otherwise identical fragments by which mate is read 1; the key here does not (the
leftmost mate of a proper pair is almost always on the forward strand, so its strand adds
little). A pair whose rightmost mate fails MAPQ is still counted
(ENCODE drops it); one whose leftmost mate fails is dropped by both. Across six
test samples this inflated the fragment total by 0.15-0.21% vs ENCODE; NRF and
PBC1 agreed within 0.0002 and PBC2 within 0.004.
"""
import argparse
import csv
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks_cwt

SUFFIX = ".mLb.mkD.sorted.bam"
PICARD_SUFFIX = ".mLb.clN.CollectMultipleMetrics.insert_size_metrics"

# Fragment-size defaults, as in ENCODE's encode_task_fraglen_stat_pe.py (overridable on the
# command line). Windows are in bp of fragment length.
ENCODE_NFR_UPPER = 150          # NFR reads: fragment < this
ENCODE_MONO_LOWER = 150         # mono-nucleosome reads: lower < fragment <= upper
ENCODE_MONO_UPPER = 300
ENCODE_PEAK_WINDOWS = {"nfr": (20, 90), "mono": (120, 250), "di": (300, 500)}  # bp, inclusive
ENCODE_CWT_WIDTH = 25           # find_peaks_cwt width (in bins; histogram bins are 1 bp)

# 0x100 secondary + 0x800 supplementary: pass as -F to keep primary reads only.
NONPRIMARY = 0x900

# 0x4 unmapped: pass as -f to keep unmapped reads only.
UNMAPPED = 0x4

# 0x400 duplicate: pass as -f to keep duplicates only. (set by Picard MarkDuplicates)
DUPLICATES = 0x400

# 0x2 properly paired (both mates mapped, correct orientation/distance apart)
PROPERLY_PAIRED = 0x2

# 0x40 first mate in pair: used to count fragments by leftmost mate only.
FIRST_MATE = 0x40

# ENCODE filter: unmapped(4) + mate unmapped(8) + secondary(256) + QC fail(512)
# + duplicate(1024) = 1804
ENCODE_FILTER = 1804

# Same as above without the duplicate bit: used for NRF/PBC (computed before dedup)
ENCODE_FILTER_KEEP_DUPS = ENCODE_FILTER - 1024

# Counts fragments by (chr, pos, tlen, strand) for the leftmost mate of each pair.
# Input is coordinate-sorted so a position's records are contiguous; the array is
# flushed whenever (chr, pos) changes, keeping memory flat.
# Per position, cnt[key] = number of read-pairs sharing (tlen, strand); chr/pos is
# implicit because cnt is cleared at each new (chr, pos). A key is one "fragment";
# cnt[key] > 1 means identical read-pairs (duplicates of that fragment).
# flush() adds each position's cnt into running totals over the whole file:
#   total    = all read-pairs (duplicates included)  -> frag_total
#   distinct = number of distinct fragments (keys)   -> frag_distinct
#   m1       = fragments seen in exactly one read-pair  (cnt == 1)
#   m2       = fragments seen in exactly two read-pairs (cnt == 2, i.e. one duplicate)
AWK_COMPLEXITY = r"""
function flush(   k) {
    for (k in cnt) {
        distinct++; total += cnt[k]
        if (cnt[k] == 1) m1++
        else if (cnt[k] == 2) m2++
    }
    delete cnt
}
BEGIN { FS = "\t" }
$9 > 0 {
    if ($3 != chr || $4 != pos) { flush(); chr = $3; pos = $4 }
    cnt[$9 SUBSEP (int($2 / 16) % 2)]++
}
END { flush(); printf "%d\t%d\t%d\t%d\n", total, distinct, m1, m2 }
"""


def run(cmd):
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"Command failed: {' '.join(cmd)}\n{res.stderr}")
    return res.stdout


def count(bam, regions, threads, include=0, exclude=0, min_mapq=0):
    cmd = ["samtools", "view", "-c", "-@", str(threads)]
    if include:
        cmd += ["-f", str(include)]
    if exclude:
        cmd += ["-F", str(exclude)]
    if min_mapq:
        cmd += ["-q", str(min_mapq)]
    cmd += [str(bam)] + regions
    return int(run(cmd).strip())


def frac(num, den):
    return num / den if den else float("nan")


def contigs(bam, mito):
    """Return (all contig names, non-mito contig names) from the BAM index."""
    names = [
        line.split("\t")[0]
        for line in run(["samtools", "idxstats", str(bam)]).splitlines()
        if line and not line.startswith("*")
    ]
    if mito not in names:
        sys.exit(f"Mito contig '{mito}' not found in {bam}; use --mito-name")
    return names, [n for n in names if n != mito]


def complexity(bam, regions, threads, min_mapq):
    pipeline = (
        f"samtools view -@ {threads} -F {ENCODE_FILTER_KEEP_DUPS} -f 2 "
        f"-q {min_mapq} {shlex.quote(str(bam))} "
        f"{' '.join(shlex.quote(r) for r in regions)} "
        f"| awk {shlex.quote(AWK_COMPLEXITY)}"
    )
    out = run(["bash", "-o", "pipefail", "-c", pipeline]).strip().split("\t")
    total, distinct, m1, m2 = (int(x) for x in out)
    return {
        "frag_total": total,
        "frag_distinct": distinct,
        "frag_m1": m1,
        "frag_m2": m2,
        "NRF": frac(distinct, total),
        "PBC1": frac(m1, distinct),
        "PBC2": frac(m1, m2),
    }


def read_insert_histogram(path):
    """Return (dense counts array indexed by insert size in bp, Picard summary dict)."""
    sizes, counts, summary, in_hist, header = [], [], {}, False, None
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("MEDIAN_INSERT_SIZE\t"):
                header = f
            elif header and not summary and line.strip():
                summary = dict(zip(header, f))
            elif line.startswith("insert_size\t"):
                in_hist = True
            elif in_hist:
                if not line.strip():
                    break
                sizes.append(int(f[0]))
                counts.append(sum(float(x) for x in f[1:] if x))
    if not sizes or "READ_PAIRS" not in summary:
        sys.exit(f"No insert size histogram or summary in {path}")
    hist = np.zeros(max(sizes) + 1)
    hist[sizes] = counts
    return hist, summary


def fragment_size_metrics(picard_file, nfr_upper, mono_lower, mono_upper, peak_windows, cwt_width):
    """ENCODE fraglen_stat_pe numbers: NFR fraction, NFR / mono ratio, fragment-length peaks."""
    keys = ["fragment_size_median", "nfr_fraction", "mono_fraction", "nfr_over_mono"] + [
        f"{k}_peak_{c}" for k in peak_windows for c in ("found", "bp")]
    if not picard_file.exists():
        print(f"WARNING: {picard_file} not found; fragment-size metrics left blank", file=sys.stderr)
        return dict.fromkeys(keys, "")
    hist, summary = read_insert_histogram(picard_file)
    total = float(summary["READ_PAIRS"])
    nfr = hist[:nfr_upper].sum()
    mono = hist[mono_lower + 1:mono_upper + 1].sum()
    # ENCODE runs the peak finder on the histogram starting at its first bin
    first = int(np.flatnonzero(hist)[0])
    peaks = [int(i) + first for i in find_peaks_cwt(hist[first:], np.array([cwt_width]))]
    row = {"fragment_size_median": summary["MEDIAN_INSERT_SIZE"],
           "nfr_fraction": frac(nfr, total), "mono_fraction": frac(mono, total),
           "nfr_over_mono": frac(nfr, mono)}
    for k, (lo, hi) in peak_windows.items():
        in_win = [x for x in peaks if lo <= x <= hi]
        best = max(in_win, key=lambda x: hist[x]) if in_win else ""
        row[f"{k}_peak_found"], row[f"{k}_peak_bp"] = bool(in_win), best
    return row


def process(bam, mito, threads, min_mapq, unique_mapq, frag):
    sample = bam.name[: -len(SUFFIX)]
    all_contigs, nonmt = contigs(bam, mito)

    # Non-MT regions plus '*' (reads with no coordinates, i.e. fully unmapped pairs)
    nonmt_regions = nonmt + ["*"]

    # --- ENCODE alignment rate: raw BAM, all contigs incl. MT, incl. duplicates ---
    raw_total = count(bam, [], threads, exclude=NONPRIMARY)
    raw_unmapped = count(bam, [], threads, include=UNMAPPED)
    raw_mapped = raw_total - raw_unmapped

    # --- Non-MT mapping categories (primary reads only) ---
    total = count(bam, nonmt_regions, threads, exclude=NONPRIMARY)
    unmapped = count(bam, nonmt_regions, threads, include=UNMAPPED)
    unique = count(bam, nonmt_regions, threads, exclude=NONPRIMARY | UNMAPPED, min_mapq=unique_mapq)
    mapped = total - unmapped
    multi = mapped - unique

    # --- Non-MT duplication rate ---
    dups = count(bam, nonmt, threads, include=DUPLICATES, exclude=NONPRIMARY | UNMAPPED)
    mapped_nonmt_regions = count(bam, nonmt, threads, exclude=NONPRIMARY | UNMAPPED)

    # --- Mito fraction (ENCODE: mito reads / all mapped reads, raw BAM) ---
    mito_reads = count(bam, [mito], threads, exclude=NONPRIMARY | UNMAPPED)

    # --- ENCODE-style duplication rate: after the ENCODE filter, MT included ---
    encode_filt = count(bam, [], threads, include=PROPERLY_PAIRED,
                        exclude=ENCODE_FILTER_KEEP_DUPS, min_mapq=min_mapq)
    encode_dups = count(bam, [], threads, include=PROPERLY_PAIRED | DUPLICATES,
                        exclude=ENCODE_FILTER_KEEP_DUPS, min_mapq=min_mapq)

    # --- Non-MT, non-duplicate ENCODE-aligned reads ---
    aligned = count(bam, nonmt, threads, include=PROPERLY_PAIRED, exclude=ENCODE_FILTER,
                    min_mapq=min_mapq)
    aligned_r1 = count(bam, nonmt, threads, include=PROPERLY_PAIRED | FIRST_MATE,
                       exclude=ENCODE_FILTER, min_mapq=min_mapq)

    row = {
        "sample": sample,
        "raw_primary_reads": raw_total,
        "raw_mapped_reads": raw_mapped,
        "alignment_rate_raw": frac(raw_mapped, raw_total),
        "nonmt_primary_reads": total,
        "nonmt_mapped_reads": mapped,
        "nonmt_unique": unique,
        "nonmt_multimapped": multi,
        "nonmt_unmapped": unmapped,
        "nonmt_unique_prop": frac(unique, total),
        "nonmt_multimapped_prop": frac(multi, total),
        "nonmt_unmapped_prop": frac(unmapped, total),
        "alignment_rate_nonmt": frac(mapped, total),
        "mito_reads": mito_reads,
        "mito_fraction": frac(mito_reads, raw_mapped),
        "encode_filtered_reads": encode_filt,
        "encode_duplicate_reads": encode_dups,
        "encode_duplication_rate": frac(encode_dups, encode_filt),
        "nonmt_duplicate_reads": dups,
        "nonmt_duplication_rate": frac(dups, mapped_nonmt_regions),
        "nonmt_nondup_aligned_reads": aligned,
        "nonmt_nondup_aligned_pairs": aligned_r1,
    }
    row.update(complexity(bam, nonmt, threads, min_mapq))
    row.update(fragment_size_metrics(frag["picard_dir"] / f"{sample}{PICARD_SUFFIX}", frag["nfr_upper"],
                                     frag["mono_lower"], frag["mono_upper"], frag["peak_windows"],
                                     frag["cwt_width"]))
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam-dir", required=True, type=Path,
                   help="Directory with *.mLb.mkD.sorted.bam "
                        "(<outdir>/bwa/merged_library)")
    p.add_argument("--out", required=True, type=Path, help="Output TSV")
    p.add_argument("--mito-name", default="MT",
                   help="Mitochondrial contig name (default: MT; UCSC: chrM)")
    p.add_argument("--min-mapq", type=int, default=30,
                   help="MAPQ cutoff for ENCODE-aligned reads and NRF/PBC "
                        "(default: 30, as ENCODE)")
    p.add_argument("--unique-mapq", type=int, default=1,
                   help="Mapped reads with MAPQ >= this are 'unique', otherwise "
                        "'multi-mapped' (default: 1; bwa gives MAPQ 0 to "
                        "multi-mappers)")
    p.add_argument("--picard-dir", type=Path,
                   help="Directory with the Picard *.insert_size_metrics files "
                        "(default: <bam-dir>/picard_metrics)")
    p.add_argument("--nfr-upper", type=int, default=ENCODE_NFR_UPPER,
                   help=f"NFR reads are fragments < this many bp (default: {ENCODE_NFR_UPPER})")
    p.add_argument("--mono-lower", type=int, default=ENCODE_MONO_LOWER,
                   help=f"Mono-nucleosome reads: lower < size <= upper (default: {ENCODE_MONO_LOWER})")
    p.add_argument("--mono-upper", type=int, default=ENCODE_MONO_UPPER,
                   help=f"Mono-nucleosome upper size, bp (default: {ENCODE_MONO_UPPER})")
    for k, (lo, hi) in ENCODE_PEAK_WINDOWS.items():
        p.add_argument(f"--{k}-peak-window", type=int, nargs=2, metavar=("LO", "HI"), default=[lo, hi],
                       help=f"bp range in which the {k} peak is looked for (default: {lo} {hi})")
    p.add_argument("--cwt-width", type=int, default=ENCODE_CWT_WIDTH,
                   help=f"find_peaks_cwt width (default: {ENCODE_CWT_WIDTH})")
    p.add_argument("--threads", type=int, default=4,
                   help="samtools threads per sample (default: 4)")
    p.add_argument("--jobs", type=int, default=1,
                   help="Samples processed in parallel (default: 1)")
    args = p.parse_args()

    bams = sorted(args.bam_dir.glob(f"*{SUFFIX}"))
    if not bams:
        sys.exit(f"No *{SUFFIX} found in {args.bam_dir}. Re-run nf-core/atacseq "
                 "with save_align_intermeds: true to publish them.")

    frag = {"picard_dir": args.picard_dir or args.bam_dir / "picard_metrics",
            "nfr_upper": args.nfr_upper, "mono_lower": args.mono_lower, "mono_upper": args.mono_upper,
            "peak_windows": {"nfr": tuple(args.nfr_peak_window), "mono": tuple(args.mono_peak_window),
                             "di": tuple(args.di_peak_window)},
            "cwt_width": args.cwt_width}
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        rows = list(ex.map(
            lambda b: process(b, args.mito_name, args.threads, args.min_mapq,
                              args.unique_mapq, frag),
            bams))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} samples to {args.out}")


if __name__ == "__main__":
    main()
