# Data dictionary: CCHMC/MiraldiLab `IBD_multiome_qc_metrics.R` output

Transcribed from [`MiraldiLab/IBDmultiome`](https://github.com/MiraldiLab/IBDmultiome)'s
[`scripts/IBD_multiome_qc_metrics.R`](https://github.com/MiraldiLab/IBDmultiome/blob/main/scripts/IBD_multiome_qc_metrics.R)
(and the doublet calls it consumes from
[`IBD_multiome_qc_atac_doublets_ArchR.R`](https://github.com/MiraldiLab/IBDmultiome/blob/main/scripts/IBD_multiome_qc_atac_doublets_ArchR.R)/
[`IBD_multiome_qc_rna_doublets.R`](https://github.com/MiraldiLab/IBDmultiome/blob/main/scripts/IBD_multiome_qc_rna_doublets.R)),
as of this repo's review of that source. One row per cell barcode, per
sample; source file `meta_<file_save>_All.txt` (and its three QC-passing
subsets, see **Output barcode sets** below). Unlike this repo's own data
dictionaries (`rnaseq_qc_metrics_data_dictionary.md`,
`atacseq_qc_metrics_data_dictionary.md`), every threshold below is a literal
filtering cutoff in the source script, not just a reported number — see
`README.md`'s note on that design difference before adapting this.

## RNA metrics (Seurat, from CellBender-denoised counts)

| Column | Definition | Threshold applied |
|---|---|---|
| `nCount_RNA` | Total UMI count for the barcode. | `1000 <= nCount_RNA <= 20000` |
| `nFeature_RNA` | Number of genes detected (nonzero count) for the barcode. | `200 <= nFeature_RNA <= 8000` |
| `gex_pct_mito` | % of UMIs from genes matching `^MT-` (mitochondrial). | `<= 25` |
| `gex_pct_ribo` | % of UMIs from genes matching `^RP[SL]` (ribosomal protein). | none (reported only) |

## ATAC metrics (Signac `ChromatinAssay` + Cell Ranger ARC's own per-barcode metrics)

| Column | Definition | Threshold applied |
|---|---|---|
| `atac_fragments` | Total ATAC fragments for the barcode (Cell Ranger ARC `per_barcode_metrics.csv`). | `2000 <= atac_fragments <= 100000` |
| `TSS.enrichment` | Signac `TSSEnrichment()`: fold-enrichment of Tn5 insertions at annotated transcription start sites vs. flanking background. | `>= 3` |
| `nucleosome_signal` | Signac `NucleosomeSignal()`: ratio of mono-nucleosome-sized to nucleosome-free-sized fragments. Values near the expected nucleosome-spacing periodicity indicate intact chromatin. | `0.2 <= nucleosome_signal <= 2` |
| `atac_pct_mito` | % of ATAC reads from the mitochondrial genome: `atac_mitochondrial_reads / (atac_mitochondrial_reads + atac_fragments) × 100`. | `<= 10` |
| `atac_pct_blacklist` | % of ATAC counts overlapping the [extended hg38 ENCODE/maxATAC blacklist](https://github.com/MiraldiLab/maxATAC) (Cazares et al.), via Signac `FractionCountsInRegion()`. | `<= 3` |
| `atac_pct_reads_in_peaks` | % of ATAC fragments falling in called peak regions: `atac_peak_region_fragments / atac_fragments × 100` (both from Cell Ranger ARC's own metrics). | `>= 20` |

## Doublet calls (computed upstream, merged in by barcode)

| Column | Definition |
|---|---|
| `isDoubletATAC` | `1` if this barcode appears in the ArchR doublet-filtered-out list (`filterRatio = 0.5`, see `README.md`), else `0`. |
| `isDoubletRNA` | `1` if this barcode was classified `Doublet` by DoubletFinder (homotypic-adjusted expected rate = 0.8% per 1,000 recovered barcodes), else `0`. |

## Identifier and grouping columns (from sample metadata, not computed)

| Column | Definition |
|---|---|
| `barcode` | Cell Ranger ARC cell barcode (10x whitelist sequence + suffix). |
| `Sample` | Sample identifier, from the pipeline's own sample metadata table. |
| `Disease`, `Tissue`, `Status` | Study design metadata merged in per sample (e.g. `CD`/`UC`, `R`/`TI`, `NI`/`I` in the source study) — cohort-specific, not computed by this pipeline. |

## Output barcode sets

The metrics script exports three tab-separated barcode lists (all built from
the single metrics table above, filtered differently):

| File | Barcodes included |
|---|---|
| `meta_<file_save>_filtered_both.txt` | Pass **every** RNA and ATAC threshold above, **and** are not doublet-flagged in both modalities (`isDoubletRNA == 0 \| isDoubletATAC == 0`) — a barcode flagged a doublet in only one modality can still pass. |
| `meta_<file_save>_filtered_scatac.txt` | Pass all ATAC thresholds and `isDoubletATAC == 0`, unioned with every barcode already in `filtered_both`. |
| `meta_<file_save>_filtered_scrna.txt` | Pass all RNA thresholds and `isDoubletRNA == 0`, unioned with every barcode already in `filtered_both`. |

`meta_<file_save>_All.txt` is the unfiltered per-barcode metrics table (every
column above, every barcode in Cell Ranger ARC's filtered matrix, before any
threshold is applied) — the source for all three filtered sets above.
