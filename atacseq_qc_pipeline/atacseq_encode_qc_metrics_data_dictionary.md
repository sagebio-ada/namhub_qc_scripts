# Data dictionary: `atacseq_encode_qc_metrics.py` output

One row per library, tab-separated, in the order of the BAM file names. Source: the duplicate-marked, coordinate-sorted library BAMs `<bam-dir>/<library>.mLb.mkD.sorted.bam` of an nf-core/atacseq run (published only with `save_align_intermeds: true`), and the Picard insert-size file `<picard-dir>/<library>.mLb.clN.CollectMultipleMetrics.insert_size_metrics` (default `<bam-dir>/picard_metrics`). The combined table written by `atacseq_qc_metrics.py` copies a subset of these columns (see [atacseq_qc_metrics_data_dictionary.md](atacseq_qc_metrics_data_dictionary.md)); this file defines all of them.

Conventions
- A library is `<sample>_REP<replicate>`, taken from the BAM file name (the name minus `.mLb.mkD.sorted.bam`).
- Counts are **reads** (each mate of a pair counts once), not pairs, unless the column says pairs or fragments. They come from `samtools view -c`.
- "Primary" reads exclude secondary and supplementary alignments (flags 0x100, 0x800). "Mapped" means the unmapped flag (0x4) is not set. "Non-MT" means all contigs except `--mito-name` (default `MT`), plus the unplaced `*` records (reads with no coordinates) wherever a column says so.
- Fractions are 0-1. A fraction whose denominator is 0 is written as `nan`.
- Defaults are the ENCODE values; `--min-mapq` (30), `--unique-mapq` (1), and the fragment-size options change the columns that say so.
- Fragment-size columns are blank (a warning is printed) if the Picard file for a library is missing.

## Identifier

| Column | Definition |
|---|---|
| `sample` | Library name, `<sample>_REP<replicate>`. Called `library` in the combined table's join key. |

## Alignment rate: raw BAM, all contigs including MT, duplicates included (approximates ENCODE)

| Column | Definition |
|---|---|
| `raw_primary_reads` | Primary reads in the BAM (mapped and unmapped). |
| `raw_mapped_reads` | `raw_primary_reads` minus primary reads flagged unmapped. |
| `alignment_rate_raw` | `raw_mapped_reads` / `raw_primary_reads`. Fraction. ENCODE reports the mapped % from SAMstats on its raw bowtie2 BAM; this is the same idea on nf-core's bwa BAM. |

## Mapping categories: non-MT contigs plus unplaced reads, primary reads only

| Column | Definition |
|---|---|
| `nonmt_primary_reads` | Primary reads on non-MT contigs plus the unplaced `*` records. |
| `nonmt_mapped_reads` | Of those, reads not flagged unmapped. |
| `nonmt_unique` | Mapped primary reads with MAPQ >= `--unique-mapq` (default 1; bwa gives MAPQ 0 to multi-mappers). |
| `nonmt_multimapped` | `nonmt_mapped_reads` minus `nonmt_unique`, i.e. mapped reads with MAPQ below `--unique-mapq`. |
| `nonmt_unmapped` | Primary reads flagged unmapped. |
| `nonmt_unique_prop`, `nonmt_multimapped_prop`, `nonmt_unmapped_prop` | The three counts above divided by `nonmt_primary_reads`. Fractions. |
| `alignment_rate_nonmt` | `nonmt_mapped_reads` / `nonmt_primary_reads`. Fraction. |

## Mitochondrial reads: raw BAM

| Column | Definition |
|---|---|
| `mito_reads` | Mapped primary reads on the mito contig. |
| `mito_fraction` | `mito_reads` / `raw_mapped_reads`. Fraction. Same idea as ENCODE's fraction of mitochondrial reads, but a different method: ENCODE counts MT reads (Rm) by aligning all reads again to a mito-only index and reports Rm / (Rn + Rm), Rn = non-MT mapped reads of the main alignment. Reads from nuclear copies of mito DNA can map to MT in that separate alignment, so ENCODE's value can be higher. |

## Duplication

| Column | Definition |
|---|---|
| `encode_filtered_reads` | Reads on all contigs, MT included, that are properly paired (0x2), not unmapped (0x4), mate-unmapped (0x8), secondary (0x100) or QC-fail (0x200), with MAPQ >= `--min-mapq` (default 30). Duplicates are kept. |
| `encode_duplicate_reads` | Of `encode_filtered_reads`, those with the duplicate flag (0x400). |
| `encode_duplication_rate` | `encode_duplicate_reads` / `encode_filtered_reads`. Fraction. Approximate: nf-core marked duplicates on the unfiltered BAM, ENCODE marks them on the filtered one. |
| `nonmt_duplicate_reads` | Mapped primary reads on non-MT contigs with the duplicate flag (the `*` records are not included). |
| `nonmt_duplication_rate` | `nonmt_duplicate_reads` / mapped primary reads on the same non-MT contigs. Fraction. |

## Usable reads: non-MT, non-duplicate

| Column | Definition |
|---|---|
| `nonmt_nondup_aligned_reads` | Non-MT reads that are properly paired, MAPQ >= `--min-mapq`, and not unmapped, mate-unmapped, secondary, QC-fail or duplicate (ENCODE "aligned" filter, flags 1804). |
| `nonmt_nondup_aligned_pairs` | The same filter restricted to first mates (0x40): approximately the number of read pairs. The filter is applied per read, so a first mate whose mate fails MAPQ is still counted. |

## Library complexity: non-MT, properly paired, MAPQ >= `--min-mapq`, duplicates kept

Each read pair is a fragment, identified by chromosome, leftmost position, fragment length (TLEN) and strand, and is counted once through its leftmost mate (TLEN > 0). The strand is that of the leftmost mate, which in a proper pair is almost always forward, so unlike ENCODE's key (R1 start, R2 end, both strands) it does not separate identical fragments by which mate is read 1. Approximate relative to ENCODE's `bedtools bamtobed -bedpe` method (see the script docstring: fragment totals about 0.15-0.21% higher, NRF/PBC1 within 0.0002, PBC2 within 0.004 on six test samples).

| Column | Definition |
|---|---|
| `frag_total` | Read pairs, duplicates included. |
| `frag_distinct` | Distinct fragments. |
| `frag_m1` | Fragments seen in exactly one read pair. |
| `frag_m2` | Fragments seen in exactly two read pairs. |
| `NRF` | Non-redundant fraction: `frag_distinct` / `frag_total`. Higher is better. |
| `PBC1` | PCR bottlenecking coefficient 1: `frag_m1` / `frag_distinct`. |
| `PBC2` | PCR bottlenecking coefficient 2: `frag_m1` / `frag_m2`. `nan` if `frag_m2` is 0. |

## Fragment size: Picard insert-size histogram of the final filtered BAM

Follows ENCODE's `encode_task_fraglen_stat_pe.py`: same size cutoffs, peak windows and peak finder. Read-pair counts are divided by `READ_PAIRS` from the Picard summary, because the histogram itself is trimmed (median + 10 x MAD). Reported only; no pass/fail (ENCODE's checks are NFR fraction >= 0.4 and NFR / mono >= 2.5). Sizes in bp.

Differences from ENCODE: ENCODE runs Picard on only the first 5,000,000 reads (`STOP_AFTER`, histogram width 1000) of its own filtered BAM (bowtie2, MAPQ >= 30) and divides by the histogram sum. Here the histogram is nf-core's, for the whole final `clN` BAM, whose filters are not identical to ENCODE's. The two denominators differed by 0.01% in the GSE195688 test.

| Column | Definition |
|---|---|
| `fragment_size_median` | Picard `MEDIAN_INSERT_SIZE`. |
| `nfr_fraction` | Read pairs with fragment length < `--nfr-upper` (default 150), divided by `READ_PAIRS`: nucleosome-free region. Fraction. |
| `mono_fraction` | Read pairs with `--mono-lower` < fragment length <= `--mono-upper` (default 150 < size <= 300), divided by `READ_PAIRS`: mono-nucleosome. Fraction. |
| `nfr_over_mono` | NFR read pairs / mono-nucleosome read pairs. |
| `nfr_peak_found` | `True` if a peak of the histogram, found with scipy `find_peaks_cwt` (width `--cwt-width`, default 25) starting from the first non-empty bin, lies in `--nfr-peak-window` (default 20-90 bp, inclusive); otherwise `False`. |
| `nfr_peak_bp` | Fragment length of that peak (the tallest if several fall in the window); blank if none. |
| `mono_peak_found`, `mono_peak_bp` | Same for `--mono-peak-window` (default 120-250 bp). |
| `di_peak_found`, `di_peak_bp` | Same for `--di-peak-window` (default 300-500 bp): di-nucleosome. |
