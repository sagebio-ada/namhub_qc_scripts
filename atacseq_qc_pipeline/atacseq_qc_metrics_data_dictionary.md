# Data dictionary: `atacseq_qc_metrics.py` output

Tab-separated, with two kinds of rows joined on `library`: one row per sequencing run (pre-merge columns: fq lint, read length, contamination, read counts and FastQC) and one row per library (post-merge columns: alignment, fragment size, FRiP, TSS enrichment). The columns of the other kind are blank; nothing is aggregated across runs. A library with one samplesheet row has one run row and one library row. Source: an nf-core/atacseq run (bwa, narrow peaks) plus the table from `atacseq_encode_qc_metrics.py`. Paths are relative to `--atac-dir`; `<peak>` = `narrow_peak` (`--peak-type`).

Conventions
- A library is `<sample>_REP<replicate>` from the samplesheet (e.g. `BMP_REP1`), as in nf-core's file names. Samplesheet rows of one library are technical replicates; run *n* (sheet order, from 1) is `<library>_T<n>`, as nf-core names them, and nf-core merges the runs after alignment.
- `r1` / `r2` = read 1 / read 2 (`r2` blank for single-end data). `raw` = before trimming; `trimmed` = after Trim Galore.
- A blank cell means the source was missing or unparsable (a warning is printed). A whole section is blank together, except fq lint, which is blank per stage.
- Fractions are 0-1 and percentages 0-100, as stated per column.
- Column groups that are the same as in the RNA-seq table are defined in [../rnaseq_qc_pipeline/rnaseq_qc_metrics_data_dictionary.md](../rnaseq_qc_pipeline/rnaseq_qc_metrics_data_dictionary.md): `fq_lint_*`, `contam_*`, `raw_reads`, `trimmed_reads`, `trimmed_pct_retained`, and all `fastqc_*`. Only what differs for ATAC is noted here.

## Identifier

| Column | Definition |
|---|---|
| `sample` | Row id: `<library>_T<n>` for a run row, the library name for a library row. |
| `library` | Join key: `<sample>_REP<replicate>` the row belongs to. |
| `row_type` | `run` or `library`. |

## fq lint — run by `atacseq_tools_qc.py` (not part of nf-core/atacseq); read from `<qc-dir>/fq_lint/{raw,trimmed}/`

| Column | Definition |
|---|---|
| `fq_lint_raw_status`, `fq_lint_raw_records`, `fq_lint_raw_issues` | fq lint on the samplesheet FASTQs: `pass` / `fail` / `incomplete`, read pairs counted, number of ERROR/WARN lines. |
| `fq_lint_trimmed_status`, `fq_lint_trimmed_records`, `fq_lint_trimmed_issues` | Same for the trimmed FASTQs. **Blank unless `--trimmed-fastq-dir` was given to `atacseq_tools_qc.py`**, because nf-core/atacseq does not keep trimmed FASTQs by default. |

## Read length — FastQC zips (`fastqc/zips`, `trimgalore/fastqc/zips`)

| Column | Definition |
|---|---|
| `read_length_r1`, `read_length_r2` | Median read length (bp) of the raw reads. |
| `read_length_trimmed_r1`, `read_length_trimmed_r2` | Median read length (bp) of the trimmed reads. |

The MultiQC FastQC table here has no median column, so the median is computed from FastQC's Sequence Length Distribution (a length range counts at its midpoint).

## Contamination — run by `atacseq_tools_qc.py` (kraken2 + Bracken); read from `<qc-dir>/contaminants/`

`contam_unclassified_pct`, `contam_top5_species`, `contam_expected_species`, `contam_expected_pct_of_classified`, `contam_expected_pct_of_all_reads`: as in the RNA-seq table. Differences: the reads screened are the trimmed FASTQs if `--trimmed-fastq-dir` was given to `atacseq_tools_qc.py`, otherwise the **raw** reads; and only the first `--kraken-max-reads` reads (pairs) of each file are screened (0 = all; default 1,000,000). ATAC reads are whole-genome and short (about 41 bp), and the unclassified percentage is much higher than in the RNA-seq table (about 80% vs 35%); do not compare the two directly.

## Read counts and FastQC — `multiqc/<peak>/multiqc_data/multiqc_fastqc*.txt`

`raw_reads`, `trimmed_reads`, `trimmed_pct_retained` and all `fastqc_{raw,trimmed}_*` columns: as in the RNA-seq table. FastQC rows are those of the run (`<library>_T<n>_{1,2}`).

## Alignment (from `atacseq_encode_qc_metrics.py`; BAM `bwa/merged_library/<library>.mLb.mkD.sorted.bam`)

| Column | Definition |
|---|---|
| `alignment_rate_raw` | Mapped primary reads / all primary reads in the raw BAM, all contigs including MT, before dedup or filtering. Fraction. Approximates ENCODE's mapped % (SAMstats on its raw bowtie2 BAM). |
| `alignment_rate_nonmt` | Mapped / all primary reads on non-MT contigs plus unplaced (fully unmapped) reads. Fraction. |
| `nonmt_nondup_aligned_reads` | Non-MT, non-duplicate reads that pass the ENCODE filter: properly paired, MAPQ >= 30 (`--min-mapq`), not unmapped / secondary / QC-fail / duplicate. |

## Library complexity (from the ENCODE script; non-MT, properly paired, MAPQ >= 30, counted on fragments before dedup)

Each fragment is a read pair identified by chromosome, leftmost position, fragment length and strand (of the leftmost mate, so identical fragments are not separated by which mate is read 1, as they are in ENCODE). Approximate relative to ENCODE's `bedtools bamtobed` method (see [atacseq_encode_qc_metrics_data_dictionary.md](atacseq_encode_qc_metrics_data_dictionary.md)).

| Column | Definition |
|---|---|
| `NRF` | Non-redundant fraction: distinct fragments / total fragments. Higher is better. |
| `PBC1` | PCR bottlenecking coefficient 1: fragments seen exactly once / distinct fragments. |
| `PBC2` | PCR bottlenecking coefficient 2: fragments seen once / fragments seen exactly twice. |

## Fragment size (from the ENCODE script; Picard insert-size histogram `bwa/merged_library/picard_metrics/<library>.mLb.clN.CollectMultipleMetrics.insert_size_metrics`)

Follows ENCODE's `encode_task_fraglen_stat_pe.py` (same cutoffs, peak windows and peak finder), but on nf-core's histogram of the whole final BAM, where ENCODE uses the first 5,000,000 reads of its own filtered BAM. Read-pair counts are divided by Picard's `READ_PAIRS`, since the histogram itself is trimmed. Reported only; no pass/fail.

| Column | Definition |
|---|---|
| `fragment_size_median` | Picard `MEDIAN_INSERT_SIZE` (bp). |
| `nfr_fraction` | Fraction of read pairs with fragment length < 150 bp (nucleosome-free region, NFR). |
| `mono_fraction` | Fraction of read pairs with fragment length > 150 and <= 300 bp (mono-nucleosome). |
| `nfr_over_mono` | NFR read pairs / mono-nucleosome read pairs. |
| `nfr_peak_found` | `True` if a fragment-length peak (scipy `find_peaks_cwt`, width 25) falls in 20-90 bp. |
| `nfr_peak_bp` | Fragment length of that peak (the tallest if several); blank if none. |
| `mono_peak_found` | `True` if a peak falls in 120-250 bp. |
| `mono_peak_bp` | Fragment length of the mono-nucleosome peak; blank if none. |

## Peak enrichment

| Column | Definition |
|---|---|
| `frip` | Fraction of reads in peaks: fraction of mapped reads overlapping MACS2 peaks (by at least 20% of the read), from `bwa/merged_library/macs2/<peak>/qc/<library>.mLb.clN_peaks.FRiP_mqc.tsv`. Fraction. Note: the peaks are nf-core's own MACS2 calls, which differ from the ENCODE-style calls used for IDR. |
| `ataqv_hqaa_overlapping_peaks_pct` | ataqv `hqaa_overlapping_peaks_percent`: % of high-quality autosomal alignments (HQAA) that overlap peaks. Percentage. From `bwa/merged_library/ataqv/<peak>/<library>.ataqv.json`. |
| `tss_enrichment` | ataqv `tss_enrichment`: signal at transcription start sites relative to flanking regions. From the same ataqv JSON. |

## Replicate concordance (from `atacseq_idr_concordance.py`; library rows only)

Copied from the IDR table (`--idr-metrics`); full definitions in [atacseq_idr_concordance_data_dictionary.md](atacseq_idr_concordance_data_dictionary.md). The IDR table has one row per **group** (the sample, i.e. the library name without `_REP<n>`), so all replicate libraries of a group get the same values. Blank if `--idr-metrics` is not given, or (with a warning) if the group has no row there (fewer than 2 replicates). IDR counts are of peaks after the blacklist / chromosome filter. The peaks are the ENCODE-style MACS2 re-calls made by that script, not nf-core's own.

| Column | Definition |
|---|---|
| `idr_Nt` | `Nt`: IDR peaks of the true-replicate pair with the most peaks. |
| `idr_Np` | `Np`: IDR peaks of pooled pseudoreplicate 1 vs 2. |
| `idr_N_optimal` | `N_optimal`: `idr_Nt` if it exceeds `idr_Np`, otherwise `idr_Np`. The number of peaks in the IDR (optimal) peak file, the ENCODE peak-count measure. |
| `idr_rescue_ratio` | `rescue_ratio`: max(`Np`, `Nt`) / min(`Np`, `Nt`). ENCODE flags values above 2. |
| `idr_self_consistency_ratio` | `self_consistency_ratio`: max / min of the self-pseudoreplicate IDR counts over the replicates. ENCODE flags values above 2. |
| `idr_reproducibility` | `reproducibility`: `pass` (both ratios <= 2), `borderline` (one above 2) or `fail` (both above 2, or `Nt` and `Np` both 0). |
| `overlap_N_optimal` | `overlap_N_optimal`: the same as `idr_N_optimal` for the naive-overlap peaks: the number of peaks in the replicated (overlap) peak file. |
