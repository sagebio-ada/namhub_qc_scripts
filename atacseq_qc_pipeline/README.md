# ATAC-seq QC pipeline

Scripts for QC of an nf-core/atacseq run (bwa, MACS2 narrow peaks). Some metrics nf-core does not
produce, or computes differently from ENCODE, so there are two metric scripts plus the combined
table.

| Script | What it does |
|---|---|
| `atacseq_encode_qc_metrics.py` | ENCODE-style metrics that nf-core does not report (alignment rates, NRF/PBC, fragment-size metrics, duplication) from the duplicate-marked BAMs and Picard insert-size files. Needs `save_align_intermeds: true`. Columns: [atacseq_encode_qc_metrics_data_dictionary.md](atacseq_encode_qc_metrics_data_dictionary.md). |
| `atacseq_qc_metrics.py` | **The combined table**: one row per sequencing run (fq lint, read length, kraken2/Bracken contamination, FastQC) plus one row per library (selected columns of the ENCODE table, FRiP, TSS enrichment, and optionally the IDR results), joined on `library`. Every column is defined in [atacseq_qc_metrics_data_dictionary.md](atacseq_qc_metrics_data_dictionary.md). |
| `atacseq_idr_concordance.py` | IDR / peak concordance between replicates, re-calling peaks the ENCODE way (its main results can be added to the combined table with `--idr-metrics`; columns: [atacseq_idr_concordance_data_dictionary.md](atacseq_idr_concordance_data_dictionary.md)). `macs2_cluster.sh` is the MACS2 wrapper it uses. |
| `atacseq_tools_qc.py` | fq lint and kraken2/Bracken, which are not part of nf-core/atacseq (through [`../qc_common/tools.py`](../qc_common/tools.py)). Writes the files the combined table reads the `fq_lint_*` and `contam_*` columns from. |
| `submit_atacseq_qc.sh` | Submits all of the above as separate Slurm jobs that run at the same time, with the combined table last. |

`atacseq_common.py` holds the samplesheet and FastQC helpers used by `atacseq_qc_metrics.py` and
`atacseq_tools_qc.py`.

## Running

Create the environment once, then submit the jobs (edit the paths at the top of the script first):

```
conda env create -f environment.yml          # namhub_atac_qc: fq, kraken2, bracken, samtools, bedtools, bedClip, idr, numpy, scipy
bash submit_atacseq_qc.sh                    # output in <atac-dir>/qc_metrics
QC_DIR=/path/to/out bash submit_atacseq_qc.sh   # output elsewhere
```

`submit_atacseq_qc.sh` only calls `sbatch`, so it can be run from any node. It submits: the ENCODE
metrics (one job), IDR (one job per group with at least 2 replicates), fq lint and kraken2/Bracken
(a job array, one task per library), and the combined table, which starts when all of those have
finished successfully (`--dependency=afterok`). Resources for each are set at the top of the script.
Logs are in `<QC_DIR>/logs/`. If one job fails, the table job is cancelled; fix the cause and run the
script again, since finished work is skipped.

The scripts are heavy (BAM processing, kraken2 loads a ~16 GB database): run them through Slurm,
not on the login node. Inputs: the nf-core output directory, the samplesheet used for the run,
and a Kraken2 database with Bracken files.

## Inputs and arguments

`submit_atacseq_qc.sh` runs all the scripts below; the paths and the few settings at its top
(`ATAC_DIR`, `SAMPLESHEET`, `KRAKEN_DB`, `KRAKEN_MAX_READS`, `MITO`, `GENOME_FAI`, `KEEP_CHR_REGEX`,
`BLACKLIST`, and the resources per job) are what you edit. `BLACKLIST` should be the blacklist
given to nf-core (its `params.yaml`); it is on by default, and empty turns it off. Run `<script> --help` for the full list. The IDR
results go to `<QC_DIR>/idr/<group>/`.

### `atacseq_encode_qc_metrics.py`

Reads the duplicate-marked BAMs and Picard insert-size files that nf-core publishes under
`bwa/merged_library` (needs `save_align_intermeds: true`), so it works on merged libraries.

| Argument | Meaning |
|---|---|
| `--bam-dir` (required) | Directory with `*.mLb.mkD.sorted.bam`, normally `<atac-dir>/bwa/merged_library`. |
| `--out` (required) | Output TSV, one row per library. |
| `--mito-name` | Mitochondrial contig name. Default `MT` (Ensembl); `chrM` for UCSC. |
| `--min-mapq` | MAPQ cutoff for ENCODE-aligned reads and NRF/PBC. Default 30. |
| `--unique-mapq` | Mapped reads with MAPQ at or above this are "unique", the rest "multi-mapped". Default 1. |
| `--picard-dir` | Picard `*.insert_size_metrics` files. Default `<bam-dir>/picard_metrics`. |
| `--nfr-upper`, `--mono-lower`, `--mono-upper` | Fragment-size cutoffs (bp) for the NFR and mono-nucleosome fractions. Defaults 150, 150, 300 (ENCODE). |
| `--nfr-peak-window`, `--mono-peak-window`, `--di-peak-window` `LO HI` | bp range searched for each fragment-size peak. Defaults 20 90, 120 250, 300 500. |
| `--cwt-width` | `find_peaks_cwt` width. Default 25. |
| `--threads` | samtools threads per library. Default 4. |
| `--jobs` | Libraries processed in parallel. Default 1. |

The ENCODE defaults should be kept for comparable numbers; the fragment-size options only exist
for experiments.

### `atacseq_qc_metrics.py` (combined table)

| Argument | Meaning |
|---|---|
| `--atac-dir` (required) | The nf-core/atacseq output directory. |
| `--samplesheet` (required) | The `samplesheet.csv` given to nf-core/atacseq for this run. Columns read: `sample`, `fastq_1`, `replicate` and optionally `fastq_2`. Rows with the same `sample` and `replicate` are technical replicates (see below). |
| `--encode-qc` (required) | The TSV written by `atacseq_encode_qc_metrics.py`. |
| `--idr-metrics` | One or more `concordance_metrics.tsv` files from `atacseq_idr_concordance.py`. Adds the `idr_*` and `overlap_N_optimal` columns to the library rows (joined by group = sample, so every replicate of a group gets the same values). Without it those columns are blank. |
| `--out` (required) | Output TSV. |
| `--qc-dir` | Where `atacseq_tools_qc.py` wrote the fq lint and kraken2/Bracken output. Default `<atac-dir>/qc_metrics`. Where it is missing, those columns are blank and a warning is printed. |
| `--peak-type` | `narrow_peak` (default) or `broad_peak`: which MultiQC, MACS2 and ataqv outputs to read. |
| `--species` | Expected species for the contamination check. Default `Homo sapiens`. |
| `--samples` | Only these libraries (`<sample>_REP<replicate>`). Default: all in the samplesheet. |

**Libraries and runs.** A library is `<sample>_REP<replicate>`. Samplesheet rows sharing both are
technical replicates (sequencing runs), which nf-core names `<library>_T1`, `_T2`, ... in sheet
order and merges after alignment. The table has a row per run (fq lint, read length, contamination,
FastQC) and a row per library (ENCODE, FRiP, TSS), joined on the `library` column.

### `atacseq_tools_qc.py` (fq lint, kraken2/Bracken)

Needs `fq`, `kraken2` and `bracken` on PATH (all in `namhub_atac_qc`). Skips output that already
exists, so a rerun resumes. Run it before the combined table, with the same `--qc-dir`.

| Argument | Meaning |
|---|---|
| `--samplesheet` (required) | The `samplesheet.csv` given to nf-core/atacseq (same file as for the combined table). |
| `--atac-dir` | The nf-core/atacseq output directory (for the FastQC tables, which give the read length that picks the Bracken k-mer file). Required except with `--list-libraries`. |
| `--qc-dir` | Output directory. Default `<atac-dir>/qc_metrics`. |
| `--peak-type` | `narrow_peak` (default) or `broad_peak`: which MultiQC output to read. |
| `--trimmed-fastq-dir` | Trim Galore FASTQs, if you kept them. Without it the trimmed fq lint columns are blank and contamination uses the raw reads. |
| `--kraken-db` | Kraken2 database directory containing Bracken `database*mers.kmer_distrib` files. Without it kraken2/Bracken is not run and the `contam_*` columns are blank. |
| `--kraken-max-reads` | Screen only the first N reads (pairs) of each file; 0 = all. Default 1,000,000. |
| `--threads` | kraken2 threads. Default 8. |
| `--samples` | Only these libraries (`<sample>_REP<replicate>`). Default: all. `submit_atacseq_qc.sh` runs one library per Slurm array task. |
| `--list-libraries` | Print the library names and exit. |

### `atacseq_idr_concordance.py`

Re-calls peaks the ENCODE way and measures IDR / peak concordance between the replicates of each
group. It reads the `*_REP<n>.mLb.clN.sorted.bam` files, treats the part of the name before
`_REP<n>` as the group, and works on groups with at least two replicates. It needs `samtools`,
`bedtools`, `bedClip` and `idr` on PATH (all in the `namhub_atac_qc` environment) and MACS2, which comes
from the cluster module through `macs2_cluster.sh` because the conda build does not load here.

| Argument | Meaning |
|---|---|
| `--bam-dir` (required) | Directory with the `*_REP<n>.mLb.clN.sorted.bam`, normally `<atac-dir>/bwa/merged_library`. |
| `--out-dir` (required) | Working directory; metrics go to `<out-dir>/concordance_metrics.tsv`. Use a different one per concurrent run. |
| `--chrom-sizes` (required) | 2-column chromosome sizes file with the same contig names as the BAMs (e.g. `cut -f1,2` of the genome `.fai`). |
| `--groups` | Only these groups. Default: all with at least 2 replicates. |
| `--min-mapq` | Use only reads with MAPQ at or above this; for paired-end data a pair is kept only if both mates pass (ENCODE's filter, with `fixmate -r`). Default 30, as ENCODE. nf-core's `clN` BAMs only drop MAPQ 0 (about 5% of reads had MAPQ 1-29 in the GSE195688 test). `0` = no extra filter. |
| `--macs2` | MACS2 executable. Default `macs2` on PATH; on this cluster use `macs2_cluster.sh`. |
| `--gensz` | MACS2 genome size. Default `hs`. |
| `--keep-chr-regex` | Keep IDR peaks on matching chromosomes. Default `^([0-9]+\|X\|Y)$` (Ensembl names); for UCSC names use `^chr([0-9]+\|X\|Y)$`; `''` keeps all. |
| `--blacklist` | Blacklist BED (same chromosome names as the BAMs). As in ENCODE, IDR and overlap peaks that overlap it by even 1 bp are dropped. nf-core already removes the reads in blacklisted regions from the `clN` BAMs, so only peaks reaching into a blacklisted region from outside are affected (about 0.02% of IDR peaks in the GSE195688 test). Optional for the script, which warns without it; `submit_atacseq_qc.sh` passes nf-core's blacklist by default. ENCODE's own hg38 blacklist may be a different version. |
| `--pval-thresh`, `--smooth-win`, `--cap-num-peak`, `--idr-thresh` | MACS2 p-value, MACS2 `--extsize`, peaks kept, and IDR threshold. Defaults 0.01, 150, 300000, 0.05 (ENCODE). |
| `--seed` | Pseudoreplication seed. 0 (default) = size of the tagAlign in bytes, as ENCODE. |
| `--chroms` | Restrict the BAMs to these contigs. For quick tests only. |
| `--threads`, `--sort-mem`, `--jobs` | samtools threads (default 4), samtools sort memory per thread (default 1G), and steps run in parallel (default 2). |

Steps whose output exists are skipped, so an interrupted run resumes in the same `--out-dir`. The
settings that change results (`--min-mapq`, the MACS2 and IDR settings, the seed, the filters) are
saved in `<out-dir>/<group>/settings.json`; a rerun with different settings, or into a folder made
before `--min-mapq` existed, stops with a message. Use a new `--out-dir` for that.

**Differences from ENCODE.** The input is nf-core's final `clN` BAM, not ENCODE's filtered BAM.
`--min-mapq 30` applies ENCODE's MAPQ filter, but nf-core also removes reads ENCODE keeps
(soft-clipped, more than 4 mismatches, inserts over 2 kb) and aligns with bwa, not bowtie2, so peak
counts are close to but not the same as ENCODE's. When a comparison has no IDR peaks, ENCODE's
pipeline stops with an error; this script reports the 0 and sets `reproducibility` to `fail` if
`Nt` and `Np` are both 0 (a ratio of 0/0 counts as above 2).

`submit_atacseq_qc.sh` runs one IDR job per group, each with its own `--out-dir`
(`<QC_DIR>/idr/<group>/`), which is the pattern to copy for running the script on its own.

## Test data

The pipeline has been tested on one public dataset so far:

- **[GSE195688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195688)** (PubMed
  [36920035](https://pubmed.ncbi.nlm.nih.gov/36920035/)): "Lineage-Specific Differences and
  Inference of Regulatory Networks Governing Human Chondrocyte Development". Human, in-vitro
  chondrocytes derived from H9 cells. The ATAC-seq part is 6 paired-end libraries: 2 treatment
  groups x 3 biological replicates, one sequencing run each.

  | Group (`sample`) | `replicate` | SRA experiment | GEO sample (title) |
  |---|---|---|---|
  | `BMP` | 1 | SRX13980378 | GSM5848465 (BMP1-1-BR2_S6) |
  | `BMP` | 2 | SRX13980379 | GSM5848466 (BMP3-1-BR1_S3) |
  | `BMP` | 3 | SRX13980380 | GSM5848467 (BMP3-A-BR3_S4) |
  | `TGF_beta` | 1 | SRX13980381 | GSM5848468 (TGF1-1-BR2_S4) |
  | `TGF_beta` | 2 | SRX13980382 | GSM5848469 (TGF3-1-BR1_S1) |
  | `TGF_beta` | 3 | SRX13980383 | GSM5848470 (TGF3-B-BR3_S3) |

  The replicate numbers are the samplesheet order, not the `BR<n>` in the GEO titles. The FASTQs
  were downloaded with nf-core/fetchngs, and nf-core/atacseq was run with bwa on Ensembl GRCh38
  (release 116, `mito_name: MT`, ENCODE hg38 blacklist v3 renamed to Ensembl contigs,
  `save_align_intermeds: true`). The two GEO mouse ATAC-seq samples (GSM5850925, GSM5850926) were
  not used.

Tested so far: the full pipeline end to end on this dataset (ENCODE metrics, IDR for both groups,
fq lint and kraken2/Bracken on all reads, combined table). Technical replicates (several runs per
library) have only been tested on hand-made inputs, and single-end data has not been tested.

**Example output.** [`test_outputs/`](test_outputs/) holds the result tables of this test, from a
run of `submit_atacseq_qc.sh` with the settings in this repo (kraken2/Bracken on all reads; IDR with
`--min-mapq 30`, nf-core's blacklist, and the other ENCODE defaults). They are reference copies for
seeing what the output looks like; the scripts never write there. Bulky intermediate files
(tagAlign, peaks, kraken2 reports) are not included.

| File | Written by |
|---|---|
| `atacseq_qc_metrics.tsv` | `atacseq_qc_metrics.py`, the combined table (12 rows: 6 runs, 6 libraries) |
| `atacseq_encode_qc_metrics.tsv` | `atacseq_encode_qc_metrics.py`, all ENCODE-style columns |
| `idr_BMP_concordance_metrics.tsv`, `idr_TGF_beta_concordance_metrics.tsv` | `atacseq_idr_concordance.py`, one per group (`<out-dir>/concordance_metrics.tsv`, renamed) |

## Things to know

- **Trimmed FASTQs are not kept** by nf-core/atacseq unless it is run with `save_trimmed`, so
  `fq_lint_trimmed_*` are blank and the contamination screen uses the raw reads. Pass
  `--trimmed-fastq-dir` to `atacseq_tools_qc.py` if you have them.
- **Contamination is on a subsample by default** in `atacseq_tools_qc.py` (`--kraken-max-reads`,
  first 1,000,000 pairs per file). nf-core/rnaseq screens all reads, so `submit_atacseq_qc.sh`
  sets it to 0 (all reads).
- **Technical replicates** (several samplesheet rows with the same `sample` and `replicate`) are
  supported. nf-core names the runs `<library>_T1`, `_T2`, ... in sheet order and merges them after
  alignment. The combined table has a row per run for the pre-merge columns (fq lint, read length,
  contamination, FastQC) and a row per library for the post-merge ones (ENCODE, FRiP, TSS); the
  other columns of each row are blank, and nothing is aggregated. Join the two on `library`.
- MultiQC here has no median read length, so it is computed from the FastQC zips.
