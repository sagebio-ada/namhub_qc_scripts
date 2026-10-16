# Data dictionary: `rnaseq_qc_metrics.py` output

One row per sample, tab-separated. Source: an nf-core/rnaseq run (aligner `star_salmon`). Paths below are relative to `--rnaseq-dir`; `<al>` = `star_salmon/` (the `--aligner` value) and `<mqc>` = `multiqc/<aligner>/multiqc_report_data/`.

Conventions
- `r1` / `r2` = read 1 / read 2. `r2` is blank for single-end data.
- `raw` = reads before trimming; `trimmed` = reads after adapter/quality trimming.
- A blank cell means the source file was missing or unparsable for that sample (a warning is printed to stderr). A whole section is blank together.
- Percentages are on a 0-100 scale unless noted (`picard_duplication_rate` is a 0-1 fraction).

## Identifier

| Column | Definition |
|---|---|
| `sample` | Sample name. Taken from the `sample` column of the nf-core samplesheet (`--samplesheet`; replicate rows collapse to one sample). Without a samplesheet, taken from the STAR files `<al>/log/<sample>.Log.final.out`. |

## fq-lint (file integrity) — `fq_lint/{raw,trimmed}/<sample>.fq_lint.txt`

| Column | Definition |
|---|---|
| `fq_lint_raw_status`, `fq_lint_trimmed_status` | `pass` = no ERROR/WARN lines and the log ends with "fq-lint end"; `fail` = at least one ERROR/WARN line; `incomplete` = the run never finished. |
| `fq_lint_raw_records`, `fq_lint_trimmed_records` | Read pairs (paired-end) or reads (single-end) counted by fq-lint. |
| `fq_lint_raw_issues`, `fq_lint_trimmed_issues` | Number of ERROR/WARN lines in the fq-lint log. |

## Read length — `<mqc>/multiqc_fastqc_fastqc_{raw,trimmed}.txt`

| Column | Definition |
|---|---|
| `read_length_r1`, `read_length_r2` | Median read length (bp) of the **raw** reads for read 1 / read 2, from the raw FastQC table. Blank if that table has no row for the sample; `r2` is blank for single-end data. |
| `read_length_trimmed_r1`, `read_length_trimmed_r2` | Median read length (bp) of the **trimmed** reads for read 1 / read 2, from the trimmed FastQC table. Blank if that table has no row for the sample; `r2` is blank for single-end data. |

## Strandedness — `<mqc>/multiqc_strand_check_summary_table.txt`

| Column | Definition |
|---|---|
| `strand_provided` | Strandedness given in the samplesheet (e.g. `forward`, `reverse`, `unstranded`, `auto`). |
| `strand_salmon` | Strandedness inferred by Salmon. |
| `strand_rseqc` | Strandedness inferred by RSeQC `infer_experiment`. |
| `strand_inferred_agree` | `True` if Salmon and RSeQC inferred the same strandedness and it is not `undetermined`; otherwise `False`. |
| `strand_status` | The pipeline's own strand-check result (e.g. pass/fail), reported as is. |

## Contamination — kraken2 report and Bracken TSV (`<al>/contaminants/`)

| Column | Definition |
|---|---|
| `contam_unclassified_pct` | % of all reads kraken2 could not classify (unclassified / (unclassified + classified)). |
| `contam_top5_species` | Five most abundant species by Bracken-estimated reads, as `name:reads(% of classified reads)`, separated by `; `. |
| `contam_expected_species` | The species expected for the sample (`--species`, default *Homo sapiens*). |
| `contam_expected_pct_of_classified` | Expected species' Bracken-estimated reads as % of classified (Bracken-estimated) reads. |
| `contam_expected_pct_of_all_reads` | Expected species' Bracken-estimated reads as % of all reads (classified + unclassified). |

## Read counts and retention — FastQC "Total Sequences" of read 1

| Column | Definition |
|---|---|
| `raw_reads` | Reads (read pairs for paired-end) before trimming. |
| `trimmed_reads` | Reads (read pairs) after trimming. |
| `trimmed_pct_retained` | `trimmed_reads / raw_reads × 100`. |

## FastQC per-read metrics — `{raw,trimmed}` × `{r1,r2}`

| Column | Definition |
|---|---|
| `fastqc_{raw,trimmed}_gc_pct_{r1,r2}` | FastQC %GC: average GC content of the reads. |
| `fastqc_{raw,trimmed}_dedup_pct_{r1,r2}` | FastQC "total deduplicated percentage": % of reads remaining if duplicates were removed (low = highly duplicated). |
| `fastqc_{raw,trimmed}_poor_quality_seqs_{r1,r2}` | Number of sequences FastQC flagged as poor quality. |

## FastQC module status — `fastqc_{raw,trimmed}_<module>`

Each cell is `pass`, `warn` or `fail`, written `r1/r2` for paired-end data (single value for single-end). Modules (one column each per stage):

| Module column suffix | FastQC module / what it checks |
|---|---|
| `basic_statistics` | File type, encoding, read count, length, %GC. |
| `per_base_sequence_quality` | Phred quality at each position along the read. |
| `per_sequence_quality_scores` | Distribution of mean quality per read. |
| `per_base_sequence_content` | A/C/G/T proportions at each position (bias suggests composition or priming artifacts). |
| `per_sequence_gc_content` | Distribution of per-read GC vs. a theoretical normal. |
| `per_base_n_content` | Fraction of N (uncalled) bases at each position. |
| `sequence_length_distribution` | Whether all reads are the same length. |
| `sequence_duplication_levels` | Degree of sequence duplication. |
| `overrepresented_sequences` | Sequences making up an unusually large share of reads. |
| `adapter_content` | Cumulative % of reads containing adapter sequence by position. |

## Salmon fragment length — `<al>/<sample>/logs/salmon_quant.log`

| Column | Definition |
|---|---|
| `salmon_fraglen_mean` | Mean fragment (insert) length (bp) Salmon estimated/used. |
| `salmon_fraglen_sd` | Standard deviation of that fragment length distribution (bp). |

## STAR alignment — `<al>/log/<sample>.Log.final.out`

| Column | Definition |
|---|---|
| `star_input_reads` | Number of input reads (pairs for paired-end) given to STAR. |
| `star_total_aligned` | Uniquely mapped + mapped to multiple loci (ENCODE definition; excludes "too many loci"). |
| `star_total_aligned_pct` | Percent of input reads in `star_total_aligned` (unique % + multi %). |
| `star_unique_n`, `star_unique_pct` | Reads mapping to exactly one locus; count and %. |
| `star_multi_n`, `star_multi_pct` | Reads mapping to multiple loci (up to `--outFilterMultimapNmax 20`); count and %. |
| `star_too_many_loci_n`, `star_too_many_loci_pct` | Reads mapping to more loci than the limit (reported separately, not counted as aligned); count and %. |
| `star_unmapped_n`, `star_unmapped_pct` | Unmapped reads: too many mismatches + too short + other; count and %. Percentages are summed from STAR's rounded values. |

## Qualimap — `<al>/qualimap/<sample>/rnaseq_qc_results.txt`

| Column | Definition |
|---|---|
| `qualimap_exonic_pct` | % of aligned reads falling in exons. |
| `qualimap_intronic_pct` | % of aligned reads falling in introns. |
| `qualimap_intergenic_pct` | % of aligned reads falling between annotated genes. |
| `qualimap_read_pairs_aligned` | Read pairs aligned, as reported by Qualimap. |

## featureCounts — `<al>/featurecounts/<sample>.featureCounts.tsv.summary`

nf-core runs featureCounts with `-p` but not `--countReadPairs`, so for paired-end data each mate is counted (about 2× the read pairs).

| Column | Definition |
|---|---|
| `featurecounts_assigned` | Alignments assigned to a gene (`Assigned`). |
| `featurecounts_no_features` | Alignments overlapping no annotated feature (`Unassigned_NoFeatures`). |
| `featurecounts_ambiguous` | Alignments overlapping more than one feature (`Unassigned_Ambiguity`). |
| `featurecounts_multimapping` | Multi-mapping alignments, not counted as assigned (`Unassigned_MultiMapping`). |
| `featurecounts_assigned_pct` | `assigned / all counted alignments × 100` (sum of all categories in the summary file). |

## rRNA — `<al>/featurecounts/<sample>.biotype_counts_*mqc.tsv`

| Column | Definition |
|---|---|
| `rrna_pct_nfcore` | nf-core's % rRNA (`percent_rRNA`); a percentage, based only on the exact `rRNA` biotype. |
| `rrna_reads` | Reads assigned to the `rRNA` biotype. |
| `mt_rrna_reads` | Reads assigned to the `Mt_rRNA` (mitochondrial rRNA) biotype. |
| `rrna_plus_mt_pct` | `(rRNA + Mt_rRNA) / all biotype-assigned reads × 100`. |

## Picard MarkDuplicates — `<al>/picard_metrics/<sample>.markdup.sorted.metrics.txt`

| Column | Definition |
|---|---|
| `picard_duplication_rate` | `PERCENT_DUPLICATION`: fraction (0-1, not %) of aligned reads that are duplicates. |
| `picard_estimated_library_size` | `ESTIMATED_LIBRARY_SIZE`: estimated number of distinct molecules in the library; blank for single-end data. |

## Replicate concordance — only with `--replicate-sheet`; `star_salmon/salmon.merged.gene_tpm.tsv`

Replicate sheet: a CSV with columns `sample, group, donor` and optionally `replicate` (a label the script does not use). `group` is the condition being compared; `donor` is the genetic source of the sample (donor, individual or cell line). Samples with the same `group` are compared pairwise: the same `donor` makes the pair **isogenic** (replicates from the same genetically identical source), different donors make it **anisogenic** (replicates from different donors, as in ENCODE). One group can have both kinds of pairs. `replicate_metadata_GSE195688.csv` is a temporary sheet for GSE195688 (donor = the H9 cell line for all four samples, from the sample titles).

| Column | Definition |
|---|---|
| `spearman_isogenic` | Median, over the other samples of this sample's group with the same donor, of the Spearman correlation of gene-level TPM between the two samples. Blank if there is no such partner or no replicate sheet. |
| `spearman_anisogenic` | Same, over the other samples of the group with a different donor. |

Genes used for each pair follow ENCODE's `MAD.R`: expressed in both samples (TPM > 0) and the mean of log2(TPM) of the pair above 0 (geometric mean TPM > 1). ENCODE applies this to RSEM FPKM; Salmon TPM is used here, so the cutoff is close to ENCODE's but not identical. Values are reported only, with no pass/fail; for reference, ENCODE's standard is > 0.9 for isogenic and > 0.8 for anisogenic replicates.

Every pair is also written to `<out>_replicate_pairs.tsv` with columns `sample_a, sample_b, group, relationship, n_genes, spearman` (`n_genes` = genes kept by the filter).
