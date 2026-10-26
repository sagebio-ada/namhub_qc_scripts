#!/usr/bin/env python
"""Replicate concordance (IDR) metrics for ATAC-seq, following the ENCODE pipeline.

Input: <outdir>/bwa/merged_library/<group>_REP<n>.mLb.clN.sorted.bam (+ .bai), the
final filtered, de-duplicated library BAMs of nf-core/atacseq. Replicates of a group
are found from the file names; groups with fewer than 2 replicates are skipped.

Requires on PATH: samtools, bedtools, bedClip, idr, bash, shuf, openssl, split (conda env
`namhub_atac_qc`), plus MACS2 via --macs2 (macs2_cluster.sh; the conda macs2 does not load here).

Follows ENCODE-DCC/atac-seq-pipeline (encode_task_bam2ta / spr / pool_ta / macs2_atac /
idr / reproducibility, and atac.wdl). Per group:
  1. BAM -> tagAlign, Tn5-shifted (+ strand start +4, - strand end -5). Paired-end: one
     line per mate; single-end: one line per read. Layout is detected from the BAMs.
     First, ENCODE's read filter (--min-mapq, default 30): -F 1804 -q 30, and for paired-end
     -f 2 then fixmate -r, so a pair is kept only if both mates have MAPQ >= 30.
  2. Self-pseudoreplicates (pr1/pr2) of each replicate: pairs (SE: reads) shuffled, split in two
  3. Pooled replicates, pooled pr1, pooled pr2
  4. MACS2 peaks on every set (-p 0.01, --shift -75 --extsize 150, top 300k by p-value)
  5. IDR (rank p.value, threshold 0.05) on
       - every pair of true replicates                        -> Nt (max over pairs)
       - each replicate's pr1 vs pr2 (self-pseudoreplicates)  -> N1, N2, ...
       - pooled pr1 vs pooled pr2                             -> Np
     The reference peak list (idr --peak-list) is the pooled peaks for the replicate pairs and
     for pooled pr1 vs pr2, and the replicate's own peaks for its pr1 vs pr2.
     IDR peaks are blacklist/chromosome-filtered before they are counted (the main columns).
     The same metrics on the unfiltered IDR peaks (the protocol's NPEAKS_IDR) are written to
     the *_unfiltered columns.
  6. rescue_ratio = max(Np, Nt) / min(Np, Nt)
     self_consistency_ratio = max(N1, N2, ...) / min(N1, N2, ...)
     reproducibility: pass; borderline if either ratio > 2; fail if both > 2
     (not ENCODE: a NaN ratio, i.e. 0 peaks on both sides, counts as > 2, and Nt = Np = 0
     is a fail; ENCODE's pipeline stops with an error on an empty IDR result instead)
  7. Naive overlap (ENCODE encode_task_overlap, the "replicated peak" set): the same comparisons
     and reference lists as in 5, but a reference peak is kept if it overlaps a peak of the first
     list and a peak of the second list by >= 50% of either peak's width (no ranks, no IDR). The
     same metrics (6) are written with an overlap_ prefix.

Deviations from ENCODE
  - Input: nf-core's final clN BAMs, not ENCODE's. With --min-mapq 30 the MAPQ filter is
    ENCODE's (nf-core itself only drops MAPQ 0: about 5% of clN reads had MAPQ 1-29 in the
    GSE195688 test), but nf-core also removes reads ENCODE keeps (soft-clipped, > 4
    mismatches, insert > 2 kb, blacklisted regions), so peak counts are close to, not
    identical to, what ENCODE would get. ENCODE also maps with bowtie2, nf-core with bwa.
  - MACS2 runs without -B --SPMR (they only add bedGraph output; peak calls are unchanged).
  - Zero peaks are reported (see 6), not an error.
Every step is skipped if its output already exists, so an interrupted run can be restarted
with the same --out-dir. The settings that change results are kept in
<out-dir>/<group>/settings.json, and a rerun with different settings stops.
"""
import argparse
import csv
import json
import math
import re
import shlex
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
from pathlib import Path

BAM_RE = re.compile(r"^(?P<group>.+)_REP(?P<rep>\d+)\.mLb\.clN\.sorted\.bam$")
# Options that change the results; recorded per group in settings.json (see check_settings)
RESULT_SETTINGS = ["min_mapq", "chroms", "seed", "gensz", "pval_thresh", "smooth_win",
                   "cap_num_peak", "idr_thresh", "blacklist", "keep_chr_regex"]
TOOLS = ["samtools", "bedtools", "bedClip", "idr", "bash", "shuf", "openssl", "split"]

# bedtools bamtobed -bedpe -mate1 columns: chr1 start1 end1 chr2 start2 end2 name score strand1 strand2
# -> two tagAlign lines per pair (mate 1 then mate 2)
AWK_BEDPE_TO_TA = r"""BEGIN{OFS="\t"}{
    print $1,$2,$3,"N",1000,$9
    print $4,$5,$6,"N",1000,$10 }"""

# Single-end: bamtobed BED6 -> tagAlign (name "N", score 1000)
AWK_BED_TO_TA = r"""BEGIN{OFS="\t"}{$4="N"; $5=1000; print $0}"""

# Tn5 shift: + strand start +4, - strand end -5 (keep start < end)
AWK_TN5 = r"""BEGIN{OFS="\t"}{
    if ($6=="+") $2=$2+4; else if ($6=="-") $3=$3-5
    if ($2>=$3) { if ($6=="+") $2=$3-1; else $3=$2+1 }
    print $0 }"""

# Pair file (12 columns per line) back to two tagAlign lines
AWK_SPLIT_PAIR = r"""BEGIN{OFS="\t"}{
    print $1,$2,$3,$4,$5,$6
    print $7,$8,$9,$10,$11,$12 }"""

# Rank by p-value (already sorted), rename, clamp coordinates, keep the top `cap` peaks
AWK_OVERLAP = r"""BEGIN{FS="\t"; OFS="\t"}{
    s1=$3-$2; s2=$13-$12
    if ((s1>0 && $21/s1 >= 0.5) || (s2>0 && $21/s2 >= 0.5)) print $0 }"""

AWK_RENAME_CAP =r"""BEGIN{OFS="\t"} NR<=cap {
    $4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0
    if ($10==-1) $10=$2+int(($3-$2+1)/2.0)
    print $0 }"""


class StepTimer:
    """Logs wall time of each stage and appends it to <path> (group, step, seconds, h:m:s).

    Steps skipped because their output exists are fast, so on a restarted run the times
    cover only the work that was actually redone.
    """

    def __init__(self, group, path):
        self.group, self.path, self.name, self.t0 = group, Path(path), None, None
        if not self.path.exists():
            self.path.write_text("group\tstep\tseconds\thms\n")

    def start(self, name):
        self.stop()
        self.name, self.t0 = name, time.monotonic()
        print(f"[{self.group}] {name}", flush=True)

    def stop(self):
        if self.name is None:
            return
        sec = time.monotonic() - self.t0
        hms = time.strftime("%H:%M:%S", time.gmtime(sec))
        print(f"[{self.group}] done: {self.name} ({hms})", flush=True)
        with open(self.path, "a") as fh:
            fh.write(f"{self.group}\t{self.name}\t{sec:.1f}\t{hms}\n")
        self.name = None


def sh(cmd):
    res = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"Command failed: {cmd}\n{res.stderr}")
    return res.stdout


def q(x):
    return shlex.quote(str(x))


def produce(out, cmd):
    """Run `cmd` (which writes to @OUT@) unless `out` exists; the write is atomic."""
    out = Path(out)
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    sh(cmd.replace("@OUT@", q(tmp)))
    tmp.rename(out)
    return out


def parallel(fns, jobs):
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futures = [ex.submit(f) for f in fns]
        return [f.result() for f in futures]


def n_lines(path):
    return int(sh(f"zcat -f {q(path)} | wc -l").strip())


def frac(num, den):
    """max/min style ratio: inf if the denominator is 0, nan if both are 0."""
    if den:
        return num / den
    return float("inf") if num else float("nan")


def is_paired(bam):
    """True if the first alignment of the BAM has the paired flag (0x1) set."""
    proc = subprocess.Popen(["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True)
    line = proc.stdout.readline()
    proc.stdout.close()
    proc.terminate()
    proc.wait()
    if not line:
        sys.exit(f"No alignments in {bam}")
    return bool(int(line.split("\t")[1]) & 1)


def make_ta(bam, out, args, paired):
    """BAM -> Tn5-shifted tagAlign. With --min-mapq > 0, first apply ENCODE's filter: -F 1804
    (-f 2 if paired) -q <min-mapq>; paired-end then fixmate -r and -F 1804 -f 2 again, so a pair
    is kept only if both mates pass (ENCODE encode_task_filter rm_unmapped_lowq_reads_pe)."""
    regions = " ".join(q(c) for c in args.chroms)
    t, mq = args.threads, args.min_mapq
    sort_n = f"samtools sort -n -u -@ {t} -m {args.sort_mem} -T {q(str(out) + '.sorttmp')} -"
    if not mq:
        view = f"samtools view -b -@ {t} {q(bam)} {regions}"
    elif paired:
        view = (f"samtools view -u -@ {t} -F 1804 -f 2 -q {mq} {q(bam)} {regions} | {sort_n} "
                f"| samtools fixmate -r -u - - | samtools view -u -F 1804 -f 2 -")
    else:
        view = f"samtools view -u -@ {t} -F 1804 -q {mq} {q(bam)} {regions}"
    if paired:
        sort = "" if mq else f"| {sort_n} "  # the MAPQ filter already name-sorted
        to_ta = f"{sort}| bedtools bamtobed -bedpe -mate1 -i stdin | awk {q(AWK_BEDPE_TO_TA)}"
    else:
        to_ta = f"| bedtools bamtobed -i stdin | awk {q(AWK_BED_TO_TA)}"
    return produce(out, f"{view} {to_ta} | awk {q(AWK_TN5)} | gzip -nc > @OUT@")


def make_pseudoreps(ta, pr1, pr2, seed, paired):
    """ENCODE spr / spr_pe: shuffle reads (SE) or read pairs (PE) with a seed, split into two halves."""
    if pr1.exists() and pr2.exists():
        return
    if paired:
        nlines = int((n_lines(ta) / 2 + 1) / 2)
        join = "| sed 'N;s/\\n/\\t/' "
    else:
        nlines = (n_lines(ta) + 1) // 2
        join = ""
    seed = seed or int(sh(f"zcat -f {q(ta)} | wc -c").strip())
    prefix = str(pr1)[: -len(".pr1.tagAlign.gz")] + ".spr"
    sh(f"zcat -f {q(ta)} {join}"
       f"| shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:{seed} -nosalt </dev/zero 2>/dev/null) "
       f"| split -d -l {nlines} - {q(prefix)}.")
    for tmp, out in ((f"{prefix}.00", pr1), (f"{prefix}.01", pr2)):
        to_ta = f"awk {q(AWK_SPLIT_PAIR)} {q(tmp)}" if paired else f"cat {q(tmp)}"
        produce(out, f"{to_ta} | gzip -nc > @OUT@")
        Path(tmp).unlink()


def pool(tas, out):
    return produce(out, f"zcat -f {' '.join(q(t) for t in tas)} | gzip -nc > @OUT@")


def call_peaks(ta, out, args):
    """ENCODE macs2_atac: MACS2, sort by p-value, rename, cap, clip to chromosome sizes."""
    out = Path(out)
    if out.exists():
        return out
    work = out.parent / (out.name + ".macs2tmp")
    shiftsize = -int(round(args.smooth_win / 2.0))
    sh(f"mkdir -p {q(work)} && {q(args.macs2)} callpeak -t {q(ta)} -f BED -n x --outdir {q(work)} "
       f"-g {args.gensz} -p {args.pval_thresh} --shift {shiftsize} --extsize {args.smooth_win} "
       f"--nomodel --keep-dup all --call-summits")
    produce(out, (
        f"LC_COLLATE=C sort -k 8gr,8gr {q(work / 'x_peaks.narrowPeak')} "
        f"| awk -v cap={args.cap_num_peak} {q(AWK_RENAME_CAP)} > {q(work / 'capped.bed')} "
        f"&& bedClip {q(work / 'capped.bed')} {q(args.chrom_sizes)} {q(work / 'clipped.bed')} -truncate "
        f"&& gzip -nc {q(work / 'clipped.bed')} > @OUT@"))
    shutil.rmtree(work)
    return out


def blacklist_filter(peak, out, args):
    """ENCODE blacklist_filter: drop blacklisted and non-standard-chromosome peaks."""
    keep = (f"awk -v re={q(args.keep_chr_regex)} 'BEGIN{{OFS=\"\\t\"}} re==\"\" || $1 ~ re "
            f"{{if ($5>1000) $5=1000; print $0}}'")
    if args.blacklist:
        cmd = (f"bedtools intersect -nonamecheck -v -a <(zcat -f {q(peak)}) -b {q(args.blacklist)} "
               f"| {keep} | gzip -nc > @OUT@")
    else:
        cmd = f"zcat -f {q(peak)} | {keep} | gzip -nc > @OUT@"
    return produce(out, cmd)


def run_idr(label, peak1, peak2, pooled, idr_dir, args):
    """ENCODE idr(): returns (IDR peaks before blacklist/chromosome filter, filtered IDR peaks)."""
    prefix = idr_dir / f"{label}.idr{args.idr_thresh}"
    bfilt = Path(f"{prefix}.bfilt.narrowPeak.gz")
    if bfilt.exists():
        return Path(f"{prefix}.narrowPeak.gz"), bfilt
    idr_dir.mkdir(parents=True, exist_ok=True)
    work = idr_dir / f"{label}.idrtmp"
    work.mkdir(exist_ok=True)
    thresh_log = -math.log10(args.idr_thresh)
    sh(f"idr --samples {q(peak1)} {q(peak2)} --peak-list {q(pooled)} --input-file-type narrowPeak "
       f"--output-file {q(work / 'unthresholded.txt')} --rank p.value "
       f"--soft-idr-threshold {args.idr_thresh} --plot --use-best-multisummit-IDR "
       f"--log-output-file {q(str(prefix) + '.log')}")
    sh(f"bedClip {q(work / 'unthresholded.txt')} {q(args.chrom_sizes)} {q(work / 'clipped.txt')} -truncate")
    peaks = produce(Path(f"{prefix}.narrowPeak.gz"), (
        f"awk 'BEGIN{{OFS=\"\\t\"}} $12>={thresh_log} {{if ($2<0) $2=0; "
        f"print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' {q(work / 'clipped.txt')} "
        f"| LC_ALL=C sort | uniq | LC_ALL=C sort -grk8,8 "
        f"| cut -f 1-10 | gzip -nc > @OUT@"))
    sh(f"gzip -nc {q(work / 'clipped.txt')} > {q(str(prefix) + '.unthresholded-peaks.txt.gz')}")
    plot = work / "unthresholded.txt.png"
    if plot.exists():
        shutil.copy(plot, f"{prefix}.unthresholded-peaks.txt.png")
    shutil.rmtree(work)
    return peaks, blacklist_filter(peaks, bfilt, args)


def run_overlap(label, peak1, peak2, pooled, ov_dir, args):
    """ENCODE overlap(): returns (overlap peaks before blacklist/chromosome filter, filtered peaks).

    Reference (pooled) peaks are intersected with peak1, then the survivors with peak2; a pair
    passes if the overlap is >= 50% of either peak. Output keeps the reference peak's columns.
    """
    prefix = ov_dir / f"{label}.overlap"
    bfilt = Path(f"{prefix}.bfilt.narrowPeak.gz")
    if bfilt.exists():
        return Path(f"{prefix}.narrowPeak.gz"), bfilt
    keep = f"awk {q(AWK_OVERLAP)} | cut -f 1-10 | LC_ALL=C sort | uniq"
    peaks = produce(Path(f"{prefix}.narrowPeak.gz"), (
        f"bedtools intersect -wo -nonamecheck -a <(zcat -f {q(pooled)}) -b <(zcat -f {q(peak1)}) "
        f"| {keep} "
        f"| bedtools intersect -wo -nonamecheck -a stdin -b <(zcat -f {q(peak2)}) "
        f"| {keep} | gzip -nc > @OUT@"))
    return peaks, blacklist_filter(peaks, bfilt, args)


def reproducibility(nt_pairs, n_pr, np_):
    """ENCODE encode_task_reproducibility (multiple-replicate case)."""
    nt_label = max(nt_pairs, key=nt_pairs.get)
    nt = nt_pairs[nt_label]
    rescue = frac(max(np_, nt), min(np_, nt))
    selfc = frac(max(n_pr.values()), min(n_pr.values()))
    # A ratio is NaN when both counts are 0 (no IDR peaks); that must not pass
    rescue_bad = math.isnan(rescue) or rescue > 2.0
    selfc_bad = math.isnan(selfc) or selfc > 2.0
    status = "pass"
    if rescue_bad or selfc_bad:
        status = "borderline"
    if (rescue_bad and selfc_bad) or (nt == 0 and np_ == 0):
        status = "fail"
    if nt > np_:
        opt_set, n_opt = nt_label, nt
    else:
        opt_set, n_opt = "pooled-pr1_vs_pooled-pr2", np_
    return {"Nt": nt, "Np": np_, "N_optimal": n_opt, "N_conservative": nt,
            "optimal_set": opt_set, "conservative_set": nt_label,
            "rescue_ratio": rescue, "self_consistency_ratio": selfc,
            "reproducibility": status}


def check_settings(path, args):
    """Steps are skipped when their output exists, so a rerun in the same --out-dir must use the
    same settings. Record them in <group>/settings.json and stop on a mismatch."""
    settings = {k: getattr(args, k) for k in RESULT_SETTINGS}
    settings = {k: str(v) if isinstance(v, Path) else v for k, v in settings.items()}
    if path.exists():
        old = json.loads(path.read_text())
        diff = {k: (old.get(k), v) for k, v in settings.items() if old.get(k) != v}
        if diff:
            sys.exit(f"{path.parent}: earlier run used other settings {diff} (old, new); "
                     "use a new --out-dir or delete the group directory")
    elif any(path.parent.glob("tagAlign/*.tagAlign.gz")):
        sys.exit(f"{path.parent}: output from a run without settings.json (before --min-mapq "
                 "existed); use a new --out-dir or delete the group directory")
    else:
        path.write_text(json.dumps(settings, indent=1) + "\n")


def process_group(group, reps, args):
    """reps: list of (rep number, BAM). Returns one output row."""
    gdir = args.out_dir / group
    ta_dir, peak_dir, idr_dir = gdir / "tagAlign", gdir / "peaks", gdir / "idr"
    ov_dir = gdir / "overlap"
    ta_dir.mkdir(parents=True, exist_ok=True)
    nums = [n for n, _ in reps]
    check_settings(gdir / "settings.json", args)
    timer = StepTimer(group, gdir / "step_timings.tsv")
    t_group = time.monotonic()
    print(f"[{group}] replicates: {', '.join(f'rep{n}' for n in nums)}", flush=True)

    timer.start("layout detection")
    layouts = {is_paired(b) for _, b in reps}
    if len(layouts) > 1:
        sys.exit(f"[{group}] replicates mix single-end and paired-end BAMs")
    paired = layouts.pop()
    print(f"[{group}] {'paired' if paired else 'single'}-end", flush=True)

    timer.start("BAM -> tagAlign")
    tas = dict(zip(nums, parallel(
        [lambda n=n, b=b: make_ta(b, ta_dir / f"{group}_rep{n}.tagAlign.gz", args, paired)
         for n, b in reps],
        args.jobs)))

    timer.start("pseudoreplicates and pooled sets")
    pr = {n: (ta_dir / f"{group}_rep{n}.pr1.tagAlign.gz", ta_dir / f"{group}_rep{n}.pr2.tagAlign.gz")
          for n in nums}
    parallel([lambda n=n: make_pseudoreps(tas[n], *pr[n], args.seed, paired) for n in nums], args.jobs)
    sets = {f"rep{n}": tas[n] for n in nums}
    for n in nums:
        sets[f"rep{n}-pr1"], sets[f"rep{n}-pr2"] = pr[n]
    sets["pooled"] = pool([tas[n] for n in nums], ta_dir / f"{group}.pooled.tagAlign.gz")
    sets["pooled-pr1"] = pool([pr[n][0] for n in nums], ta_dir / f"{group}.pooled-pr1.tagAlign.gz")
    sets["pooled-pr2"] = pool([pr[n][1] for n in nums], ta_dir / f"{group}.pooled-pr2.tagAlign.gz")

    timer.start(f"MACS2 peaks on {len(sets)} sets")
    names = list(sets)
    peaks = dict(zip(names, parallel(
        [lambda s=s: call_peaks(sets[s], peak_dir / f"{group}_{s}.narrowPeak.gz", args) for s in names],
        args.jobs)))

    jobs = {}
    for i, j in combinations(nums, 2):
        jobs[f"rep{i}_vs_rep{j}"] = (peaks[f"rep{i}"], peaks[f"rep{j}"], peaks["pooled"])
    for n in nums:
        jobs[f"rep{n}-pr1_vs_rep{n}-pr2"] = (peaks[f"rep{n}-pr1"], peaks[f"rep{n}-pr2"], peaks[f"rep{n}"])
    jobs["pooled-pr1_vs_pooled-pr2"] = (peaks["pooled-pr1"], peaks["pooled-pr2"], peaks["pooled"])

    timer.start(f"IDR: {len(jobs)} comparisons")
    labels = list(jobs)
    idr_peaks = dict(zip(labels, parallel(
        [lambda l=l: run_idr(f"{group}_{l}", *jobs[l], idr_dir, args) for l in labels], args.jobs)))
    # counts: after blacklist/chromosome filter (ENCODE repo); raw: IDR-threshold peaks before
    # that filter (the protocol's NPEAKS_IDR)
    counts = {l: n_lines(p[1]) for l, p in idr_peaks.items()}
    raw = {l: n_lines(p[0]) for l, p in idr_peaks.items()}

    def split(c):
        nt = {f"rep{i}_vs_rep{j}": c[f"rep{i}_vs_rep{j}"] for i, j in combinations(nums, 2)}
        n_self = {n: c[f"rep{n}-pr1_vs_rep{n}-pr2"] for n in nums}
        return nt, n_self, c["pooled-pr1_vs_pooled-pr2"]

    nt_pairs, n_pr, np_ = split(counts)
    raw_nt_pairs, raw_n_pr, raw_np = split(raw)

    timer.start(f"overlap: {len(jobs)} comparisons")
    ov_peaks = dict(zip(labels, parallel(
        [lambda l=l: run_overlap(f"{group}_{l}", *jobs[l], ov_dir, args) for l in labels], args.jobs)))
    ov_nt_pairs, ov_n_pr, ov_np = split({l: n_lines(p[1]) for l, p in ov_peaks.items()})
    timer.start("summarize")
    row = {"group": group, "replicates": ",".join(f"rep{n}" for n in nums)}
    row.update(reproducibility(nt_pairs, n_pr, np_))
    row["Nt_per_pair"] = ";".join(f"{k}={v}" for k, v in nt_pairs.items())
    row["N_self_per_rep"] = ";".join(f"rep{n}={v}" for n, v in n_pr.items())
    r = reproducibility(raw_nt_pairs, raw_n_pr, raw_np)
    row.update({"Nt_unfiltered": r["Nt"], "Np_unfiltered": r["Np"],
                "rescue_ratio_unfiltered": r["rescue_ratio"],
                "self_consistency_ratio_unfiltered": r["self_consistency_ratio"],
                "reproducibility_unfiltered": r["reproducibility"],
                "Nt_per_pair_unfiltered": ";".join(f"{k}={v}" for k, v in raw_nt_pairs.items()),
                "N_self_per_rep_unfiltered": ";".join(f"rep{n}={v}" for n, v in raw_n_pr.items())})
    row["macs2_peaks_per_rep"] = ";".join(f"rep{n}={n_lines(peaks[f'rep{n}'])}" for n in nums)
    row["macs2_peaks_pooled"] = n_lines(peaks["pooled"])
    row.update({f"overlap_{k}": v for k, v in reproducibility(ov_nt_pairs, ov_n_pr, ov_np).items()})
    row["overlap_Nt_per_pair"] = ";".join(f"{k}={v}" for k, v in ov_nt_pairs.items())
    row["overlap_N_self_per_rep"] = ";".join(f"rep{n}={v}" for n, v in ov_n_pr.items())
    timer.stop()
    total = time.monotonic() - t_group
    print(f"[{group}] total {time.strftime('%H:%M:%S', time.gmtime(total))}", flush=True)
    with open(timer.path, "a") as fh:
        fh.write(f"{group}\tTOTAL\t{total:.1f}\t{time.strftime('%H:%M:%S', time.gmtime(total))}\n")
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam-dir", required=True, type=Path,
                   help="Directory with *_REP<n>.mLb.clN.sorted.bam "
                        "(<outdir>/bwa/merged_library)")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Working/output directory; metrics go to <out-dir>/concordance_metrics.tsv")
    p.add_argument("--chrom-sizes", required=True, type=Path,
                   help="2-column chromosome sizes file (e.g. the nf-core saved genome/*.fa.sizes)")
    p.add_argument("--min-mapq", type=int, default=30,
                   help="Use only reads with MAPQ >= this (paired-end: pairs whose mates both pass), "
                        "as ENCODE's filtered BAM (default: 30). nf-core's clN BAMs only drop MAPQ 0. "
                        "0 = no extra filter")
    p.add_argument("--blacklist", type=Path,
                   help="Optional blacklist BED applied to IDR and overlap peaks, as ENCODE does (same "
                        "chromosome naming as the BAMs): peaks overlapping it are dropped. nf-core/atacseq "
                        "already removes the reads in blacklisted regions from the clN BAMs, so only peaks "
                        "reaching into a blacklisted region from outside remain (about 0.02%% of IDR peaks "
                        "in the GSE195688 test)")
    p.add_argument("--keep-chr-regex", default=r"^([0-9]+|X|Y)$",
                   help="Keep IDR peaks on chromosomes matching this regex "
                        r"(default: Ensembl names; UCSC: ^chr([0-9]+|X|Y)$). '' keeps all")
    p.add_argument("--groups", nargs="+", help="Only these groups (default: all with >= 2 replicates)")
    p.add_argument("--macs2", default="macs2",
                   help="MACS2 executable (default: macs2 on PATH). On this cluster use "
                        "macs2_cluster.sh; the conda macs2 build does not load here")
    p.add_argument("--gensz", default="hs", help="MACS2 genome size (default: hs)")
    p.add_argument("--pval-thresh", type=float, default=0.01, help="MACS2 p-value (ENCODE: 0.01)")
    p.add_argument("--smooth-win", type=int, default=150, help="MACS2 --extsize (ENCODE: 150)")
    p.add_argument("--cap-num-peak", type=int, default=300000, help="Top N peaks kept (ENCODE: 300000)")
    p.add_argument("--idr-thresh", type=float, default=0.05, help="IDR threshold (ENCODE ATAC: 0.05)")
    p.add_argument("--seed", type=int, default=0,
                   help="Pseudoreplication seed; 0 = size of the tagAlign in bytes, as ENCODE")
    p.add_argument("--chroms", nargs="+", default=[],
                   help="Restrict BAMs to these contigs (for quick tests only)")
    p.add_argument("--threads", type=int, default=4, help="samtools threads (default: 4)")
    p.add_argument("--sort-mem", default="1G", help="samtools sort memory per thread (default: 1G)")
    p.add_argument("--jobs", type=int, default=2, help="Steps run in parallel (default: 2)")
    args = p.parse_args()

    missing = [t for t in TOOLS + [args.macs2] if shutil.which(t) is None]
    if missing:
        sys.exit(f"Missing from PATH: {', '.join(missing)} (activate the namhub_atac_qc env)")

    groups = {}
    for bam in sorted(args.bam_dir.glob("*.mLb.clN.sorted.bam")):
        m = BAM_RE.match(bam.name)
        if m:
            groups.setdefault(m["group"], []).append((int(m["rep"]), bam))
    if not groups:
        sys.exit(f"No <group>_REP<n>.mLb.clN.sorted.bam found in {args.bam_dir}")
    if args.groups:
        unknown = set(args.groups) - set(groups)
        if unknown:
            sys.exit(f"Groups not found: {', '.join(sorted(unknown))}")
        groups = {g: groups[g] for g in args.groups}

    if not args.blacklist:
        print("WARNING: no --blacklist: IDR/overlap peaks overlapping blacklisted regions are kept "
              "(ENCODE drops them); pass the blacklist given to nf-core", file=sys.stderr, flush=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for group, reps in groups.items():
        if len(reps) < 2:
            print(f"[{group}] skipped: {len(reps)} replicate (IDR needs >= 2)", flush=True)
            continue
        rows.append(process_group(group, sorted(reps), args))
    if not rows:
        print("No group with >= 2 replicates; nothing to do (no metrics file written)")
        return

    out = args.out_dir / "concordance_metrics.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} groups to {out}")


if __name__ == "__main__":
    main()
