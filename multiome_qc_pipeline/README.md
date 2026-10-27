# Multiome (snRNA + snATAC) QC Pipeline (planned)

**Status: reference documentation only — no code here yet.** This folder
documents a real, working QC pipeline shared with us by CCHMC/MiraldiLab
([`MiraldiLab/IBDmultiome`](https://github.com/MiraldiLab/IBDmultiome), the
codebase for Yang, Wayman, Denson, et al. 2026), so the approach and
thresholds are preserved as reference before any adaptation into this repo's
own conventions. Unlike `atacseq_qc_pipeline/`'s and `rnaseq_qc_pipeline/`'s
own planning-stage READMEs (hypothetical candidate checks for a tool not yet
chosen), this one transcribes an actual, already-run pipeline.

## Background, for readers new to this assay

10x Genomics' **Chromium Single Cell Multiome ATAC + Gene Expression** assay
profiles chromatin accessibility (snATAC) and gene expression (snRNA) from
the **same nucleus**, processed together by **Cell Ranger ARC**. That joint
measurement is what makes this QC genuinely different from everything else
in this repo:

- **QC operates per cell barcode, not per raw FASTQ file.** Every other
  pipeline here (`scrnaseq_qc_pipeline/`, `rnaseq_qc_pipeline/`,
  `atacseq_qc_pipeline/`) checks raw sequencing files or an aligned/quantified
  sample-level output. This pipeline starts from Cell Ranger ARC's own
  already-aligned, already-counted per-barcode output
  (`filtered_feature_bc_matrix.h5`, `atac_fragments.tsv.gz`,
  `per_barcode_metrics.csv`) and asks, for each of potentially tens of
  thousands of barcodes per sample, whether that barcode is a real,
  high-quality single nucleus. None of this repo's existing raw-FASTQ Tier 1
  checks (`fq_lint`, FastQC, Kraken2 — see
  [`CHECK_REUSABILITY.md`](../CHECK_REUSABILITY.md)) apply at this layer;
  they'd only apply further upstream, to the raw FASTQs Cell Ranger ARC
  itself consumes, which these scripts never touch.
- **A barcode can have quality signal in one modality but not the other.**
  A nucleus can yield clean RNA but poor ATAC (or vice versa), so this
  pipeline tracks RNA-passing, ATAC-passing, and jointly-passing cell sets
  separately rather than a single pass/fail per barcode.
- **Doublets (two nuclei captured under one barcode) are detected per
  modality, independently**, because the two underlying methods (see below)
  use different signal and have no shared ground truth to validate against
  each other.
- **The whole stack is R/Bioconductor** (Seurat v5, Signac, ArchR,
  DoubletFinder) — everything else in this repo is Python. Porting this
  pipeline's logic faithfully is a rewrite, not a direct translation; several
  of the R packages involved (ArchR, DoubletFinder) have no drop-in Python
  equivalent (the closest Python-ecosystem analogs would be
  [Scrublet](https://github.com/swolock/scrublet)/
  [scDblFinder](https://github.com/plger/scDblFinder)-style doublet detection
  and `scanpy`/`snapatac2` for the ATAC side — different algorithms, not
  reimplementations of ArchR's or DoubletFinder's specific methods).

## What the pipeline actually does

Four scripts, run in this order (file names as in the source repo):

1. **[`IBD_multiome_qc_atac_doublets_ArchR.R`](https://github.com/MiraldiLab/IBDmultiome/blob/main/scripts/IBD_multiome_qc_atac_doublets_ArchR.R)**
   — per sample, builds an ArchR Arrow file from the Cell Ranger ARC
   fragments file (restricted to barcodes in the filtered RNA matrix), scores
   ATAC doublets with ArchR's `addDoubletScores` (`k=10`, UMAP-embedding
   nearest neighbors, LSI method 2), filters with `filterDoublets` at
   `filterRatio = 0.5`, and exports the excluded-doublet cell list.
2. **[`IBD_multiome_qc_rna_doublets.R`](https://github.com/MiraldiLab/IBDmultiome/blob/main/scripts/IBD_multiome_qc_rna_doublets.R)**
   — per sample, loads CellBender-denoised RNA counts, applies the same
   RNA QC thresholds used in step 3 (so doublet detection only runs on
   already-plausible cells), runs a standard Seurat workflow
   (`NormalizeData`→`ScaleData`→`RunPCA`→`RunUMAP`→Louvain clustering at
   resolution 0.5), then runs
   [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder):
   `paramSweep`/`find.pK` to pick the optimal `pK`, a homotypic-adjusted
   expected-doublet count from an assumed rate of **0.8% per 1,000 recovered
   barcodes** (the 10x-published multiplet rate for this chemistry), and
   `doubletFinder` classification.
3. **[`IBD_multiome_qc_metrics.R`](https://github.com/MiraldiLab/IBDmultiome/blob/main/scripts/IBD_multiome_qc_metrics.R)**
   — the core step. Per sample: loads the ATAC peak matrix + fragments into a
   Signac `ChromatinAssay`, loads CellBender RNA counts for the same
   barcodes, computes every metric in the
   [data dictionary](multiome_qc_metrics_data_dictionary.md) below, merges in
   the doublet calls from steps 1–2, applies the fixed thresholds also listed
   there, and exports three QC-passing barcode lists: RNA-quality,
   ATAC-quality, and jointly-quality (excluding only barcodes doublet-flagged
   by **both** modalities — a barcode doublet-flagged in just one modality
   can still pass if the other modality's signal is clean).
4. **[`IBD_multiome_qc_filter.R`](https://github.com/MiraldiLab/IBDmultiome/blob/main/scripts/IBD_multiome_qc_filter.R)**
   — assembles filtered Seurat objects (joint multiome and ATAC-only) from
   the barcode lists step 3 already produced; no new QC logic of its own.

## Dependencies

R 4.2.3 in the source environment, with: Seurat (v5), Signac (v1.9.0), ArchR,
DoubletFinder, EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38, rtracklayer,
plus standard plotting packages (ggplot2, patchwork, ggpubr). Upstream of all
of this: Cell Ranger ARC v2.0.0 (10x's own joint snRNA+snATAC processing) and
[CellBender](https://github.com/broadinstitute/CellBender) v0.3.0 for ambient
RNA removal — this pipeline's RNA-side input is CellBender's *denoised*
output, not Cell Ranger's raw counts directly.

## A deliberate departure from this repo's own design rule, worth deciding on explicitly

Every pipeline in this repo so far has followed one rule: **the pipeline's
own code never asserts a PASS/FAIL verdict** — only a named tool's own
established criteria may (FastQC's own module thresholds, `fq lint`'s exit
code); everything else is reported as a plain number, because "what counts
as good" varies by sample type in ways this repo's own code has no basis to
adjudicate.

CCHMC's pipeline does the opposite, deliberately: every threshold in the
[data dictionary](multiome_qc_metrics_data_dictionary.md) is a **hard
filtering cutoff**, not just a reported number — a cell below `min_atac_tss_enrich`
is simply excluded from the output, full stop. That's a defensible choice
for this kind of per-cell QC (filtering individual low-quality nuclei out of
a large population is a different problem than judging one bulk sample's
overall quality), but it's a real design decision to make consciously if
this pipeline is ever adapted here, not something to carry over by default:
decide explicitly whether an adapted version keeps CCHMC's literal cutoffs
as real filters, or reports the same metrics INFO-only in this repo's usual
style and leaves filtering to whoever consumes the report.

## Notes for whenever this gets adapted

- This is a **new, fifth assay/modality** this repo doesn't track anywhere
  else (ATAC-seq, WGS, scRNA-seq, and bulk RNA-seq are the four in
  [`CHECK_REUSABILITY.md`](../CHECK_REUSABILITY.md); joint snRNA+snATAC
  multiome is a distinct fifth). None of that matrix's Tier 1 (raw-FASTQ)
  checks transfer to this pipeline's own scope, since it starts downstream
  of Cell Ranger ARC.
- The existing `scrnaseq_qc_pipeline/matrix_content_qc`/`matrix_download`
  checks are the closest existing precedent in spirit (both validate
  CellRanger-produced output), but far simpler: they sanity-check one
  matrix's file shape, not compute per-barcode statistical QC metrics across
  two paired modalities plus two independent doublet-detection models.
- Faithfully porting ArchR's and DoubletFinder's specific doublet-scoring
  algorithms into Python isn't realistic; adapting this would mean picking
  Python-ecosystem equivalents (see Background above) and accepting they
  won't reproduce CCHMC's exact calls, not translating this R code directly.
- All five file paths these scripts read (`dir_cellranger_arc`,
  `dir_cellbender`, the gene/blacklist reference files, sample metadata) are
  `/path/to/...` placeholders in the source repo — there's no runnable
  example here yet, only the documented logic and thresholds.
