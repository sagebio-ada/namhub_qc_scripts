# Data dictionary: `atacseq_idr_concordance.py` output

`<out-dir>/concordance_metrics.tsv`: one row per group, tab-separated. A group is the part of the library name before `_REP<n>` (BAM `<group>_REP<n>.mLb.clN.sorted.bam`); only groups with at least 2 replicates are included. The method follows the ENCODE ATAC-seq pipeline (ENCODE-DCC/atac-seq-pipeline); the steps are in the script docstring and the [README](README.md). The combined table (`atacseq_qc_metrics.py --idr-metrics`) copies `Nt`, `Np`, `N_optimal`, `rescue_ratio`, `self_consistency_ratio`, `reproducibility` and `overlap_N_optimal` from it.

Terms
- **Replicates** are `rep1`, `rep2`, ... (the `<n>` of `_REP<n>`, sorted). **pr1 / pr2** are the two pseudoreplicates made by shuffling and splitting a replicate's Tn5-shifted tagAlign in half (**self-pseudoreplicates**). **pooled** is all replicates merged; **pooled-pr1 / pooled-pr2** merge all pr1 and all pr2.
- **Reads used**: those of nf-core's final `clN` BAM with MAPQ >= `--min-mapq` (default 30, as ENCODE); for paired-end data both mates must pass. See the README for the remaining differences from ENCODE's input.
- MACS2 peaks are called on every set (p 0.01, extsize 150, top 300,000 by p-value kept). These are re-called ENCODE-style and differ from nf-core's own MACS2 peaks used for FRiP.
- **IDR peaks** are peaks passing IDR (rank p.value, threshold `--idr-thresh`, default 0.05) between two peak lists, using a reference peak list. Counts are the number of peaks in the file.
- **Comparisons** (labels used in the per-pair columns): every pair of true replicates `rep<i>_vs_rep<j>` (reference: pooled peaks); each replicate's `rep<n>-pr1_vs_rep<n>-pr2` (reference: that replicate's peaks); and `pooled-pr1_vs_pooled-pr2` (reference: pooled peaks).
- **Filtered vs unfiltered.** The main columns count IDR peaks after the blacklist / chromosome filter (`--blacklist`: peaks overlapping a blacklisted region by even 1 bp are dropped, as ENCODE does; `submit_atacseq_qc.sh` passes nf-core's blacklist by default. `--keep-chr-regex`: by default only autosomes and X/Y are kept). The `_unfiltered` columns count the IDR-threshold peaks before that filter (the protocol's `NPEAKS_IDR`).
- Ratios: `max / min` of two counts; `inf` if the smaller is 0 and the larger is not, `nan` if both are 0 (a `nan` ratio counts as above 2). ENCODE's pipeline stops with an error in these cases instead, so this handling is this script's own.
- "Per pair / per replicate" columns hold `;`-separated `label=count` text.

## Identifier

| Column | Definition |
|---|---|
| `group` | Group name (library name without `_REP<n>`). |
| `replicates` | Replicates used, comma-separated (e.g. `rep1,rep2,rep3`). |

## IDR reproducibility (filtered peaks)

| Column | Definition |
|---|---|
| `Nt` | True-replicate IDR peaks: the largest count over all pairs of true replicates. |
| `Np` | Pooled-pseudoreplicate IDR peaks: `pooled-pr1_vs_pooled-pr2`. |
| `N_optimal` | `Nt` if `Nt` > `Np`, otherwise `Np`. |
| `N_conservative` | Equal to `Nt`. |
| `optimal_set` | Comparison that gives `N_optimal`: the replicate pair with the most peaks if `Nt` > `Np`, otherwise `pooled-pr1_vs_pooled-pr2`. |
| `conservative_set` | The replicate-pair comparison that gives `Nt`. |
| `rescue_ratio` | max(`Np`, `Nt`) / min(`Np`, `Nt`). ENCODE flags values above 2. |
| `self_consistency_ratio` | max / min of the self-pseudoreplicate counts (`N_self_per_rep`) over all replicates. ENCODE flags values above 2. |
| `reproducibility` | `pass`: both ratios <= 2. `borderline`: either ratio > 2. `fail`: both ratios > 2 (ENCODE's rule). Not ENCODE: a `nan` ratio counts as > 2, and `Nt` and `Np` both 0 is `fail`. |
| `Nt_per_pair` | IDR peak count of every true-replicate pair, `rep1_vs_rep2=N;rep1_vs_rep3=N;...`. |
| `N_self_per_rep` | IDR peak count of each replicate's pr1 vs pr2, `rep1=N;rep2=N;...`. |

## IDR reproducibility (unfiltered peaks)

Same definitions as above, on the IDR peaks before the blacklist / chromosome filter. No `N_optimal`, `N_conservative` or set columns are written for these.

| Column | Definition |
|---|---|
| `Nt_unfiltered`, `Np_unfiltered` | As `Nt` and `Np`. |
| `rescue_ratio_unfiltered`, `self_consistency_ratio_unfiltered` | As `rescue_ratio` and `self_consistency_ratio`. |
| `reproducibility_unfiltered` | As `reproducibility`. |
| `Nt_per_pair_unfiltered`, `N_self_per_rep_unfiltered` | As `Nt_per_pair` and `N_self_per_rep`. |

## MACS2 peaks

| Column | Definition |
|---|---|
| `macs2_peaks_per_rep` | Peaks called on each replicate (after the 300,000 cap and chromosome-size clipping), `rep1=N;rep2=N;...`. |
| `macs2_peaks_pooled` | Peaks called on the pooled replicates. |

## Naive overlap ("replicated peaks")

The same comparisons and reference peak lists as the IDR step, but a reference peak is kept if it overlaps a peak of the first list and a peak of the second list by at least 50% of either peak's width. No ranks and no IDR. Counts are after the blacklist / chromosome filter; there are no unfiltered overlap columns. All columns are defined as their IDR namesakes, on the overlap peaks.

| Column | Definition |
|---|---|
| `overlap_Nt`, `overlap_Np`, `overlap_N_optimal`, `overlap_N_conservative` | As `Nt`, `Np`, `N_optimal`, `N_conservative`. |
| `overlap_optimal_set`, `overlap_conservative_set` | As `optimal_set`, `conservative_set`. |
| `overlap_rescue_ratio`, `overlap_self_consistency_ratio` | As `rescue_ratio`, `self_consistency_ratio`. |
| `overlap_reproducibility` | As `reproducibility`. |
| `overlap_Nt_per_pair`, `overlap_N_self_per_rep` | As `Nt_per_pair`, `N_self_per_rep`. |

## Other output files (under `<out-dir>/<group>/`)

| File | Content |
|---|---|
| `tagAlign/` | Tn5-shifted tagAlign of each replicate, its pseudoreplicates, and the pooled sets. |
| `peaks/<group>_<set>.narrowPeak.gz` | MACS2 peaks of each set. |
| `idr/<group>_<comparison>.idr<thresh>.narrowPeak.gz`, `.bfilt.narrowPeak.gz` | IDR peaks before and after the blacklist / chromosome filter, with the IDR log, unthresholded peaks and plot. |
| `overlap/<group>_<comparison>.overlap.narrowPeak.gz`, `.bfilt.narrowPeak.gz` | Naive-overlap peaks before and after the filter. |
| `settings.json` | The settings that change results (`min_mapq`, MACS2/IDR settings, seed, filters). A rerun in the same folder with other settings stops. |
| `step_timings.tsv` | Wall time of each stage (`group, step, seconds, hms`) plus a `TOTAL` row. On a restarted run only the work that was redone is timed. |
