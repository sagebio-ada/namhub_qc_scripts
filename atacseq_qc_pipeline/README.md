# ATAC-seq QC Pipeline (planned)

**Status: planning only — no code here yet.** This README is a reference list of
QC checks worth running for ATAC-seq data tracked in Synapse, in the same spirit
as [`scrnaseq_qc_pipeline/`](../scrnaseq_qc_pipeline/README.md): candidates to
evaluate and pick from when this pipeline actually gets built, not a committed
scope.

## Background, for readers new to this assay

ATAC-seq (**A**ssay for **T**ransposase-**A**ccessible **C**hromatin using
**seq**uencing) maps open/accessible chromatin genome-wide. A hyperactive Tn5
transposase simultaneously cuts open DNA and inserts sequencing adapters
("tagmentation") — so read density directly reflects chromatin accessibility,
with no separate antibody/pulldown step like ChIP-seq. A few things follow
directly from that mechanism and drive most of ATAC-seq's assay-specific QC:

- **Mitochondrial DNA has no chromatin** (no nucleosomes to block Tn5), so
  mitochondrial reads are always massively over-represented relative to their
  share of the genome — often 5-50%+ of raw reads depending on library prep.
  This is the single most commonly cited ATAC-seq QC number.
- **Fragment length encodes nucleosome positioning.** A fragment cut twice in
  a nucleosome-free region is short (<~100bp); a fragment spanning one
  nucleosome is ~200bp; two nucleosomes ~400bp; and so on. A healthy library's
  fragment-size distribution shows this periodicity as a visible sawtooth
  pattern — its *absence* is a strong signal something went wrong upstream
  (over-tagmentation, poor nuclei prep, degraded chromatin).
- **Signal should concentrate at regulatory elements**, promoters/transcription
  start sites (TSS) especially — most cell types keep a consistent core of
  promoters accessible, so TSS enrichment is a standard cross-sample quality
  bar even before looking at anything biological/comparative.
- Downstream analysis is peak-based (called with a tool like MACS2), so a
  meaningful chunk of QC only becomes possible after peak calling, not just
  from raw reads or alignment alone.

## Possible checks

| Tier | Check | Tool | What it would report | Why it matters for ATAC-seq specifically |
|---|---|---|---|---|
| Raw FASTQ | Structural validity | [`fq lint`](https://github.com/stjude-rust-labs/fq) | Malformed records, mismatched seq/quality lengths | Same generic file-integrity check as the scRNA pipeline — not assay-specific. |
| Raw FASTQ | Read quality / adapter content | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Per-base quality, adapter read-through, duplication, GC content | ATAC libraries use **Nextera** adapters (not TruSeq) — configure FastQC's adapter list accordingly or its Adapter Content module won't flag real contamination. |
| Raw FASTQ | Species/contamination screening | [Kraken2](https://github.com/DerrickWood/kraken2) | % reads classified per organism | Same tool `scrnaseq_qc_pipeline/kraken_qc.py` already wraps — directly reusable here, not assay-specific logic. |
| Alignment | Mapping rate, unique vs. multi-mapped | `samtools flagstat`, aligner's own log (e.g. Bowtie2) | % aligned, % uniquely aligned | A low unique-mapping rate flags degraded/contaminated libraries before any accessibility signal is even computable. |
| Alignment | **Mitochondrial read fraction** | `samtools idxstats` / a dedicated ATAC QC tool | % of reads mapping to chrM | *The* headline ATAC-seq QC number — extremely high values indicate poor nuclei isolation (whole-cell contamination) diluting the accessible-chromatin signal with mitochondrial noise. |
| Alignment | PCR duplicate rate | Picard `MarkDuplicates` / `samtools markdup` | % duplicate reads | Tagmentation + PCR amplification makes ATAC libraries prone to duplication at low input; distinguishing PCR duplicates from genuine independent tagmentation events at the same site matters more here than in most assays. |
| Alignment | Library complexity (NRF, PBC1, PBC2) | [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) metrics | Non-redundant fraction, PCR bottlenecking coefficients | Standard ENCODE complexity metrics — flags a library that's been PCR-over-amplified from too little starting material, a common ATAC-seq failure mode. |
| Alignment | **Fragment size distribution** | `samtools`/Picard `CollectInsertSizeMetrics` + custom periodicity check | Insert-size histogram; presence/absence of nucleosome-spacing periodicity | ATAC-seq's signature plot — the nucleosome-free/mono-/di-nucleosome sawtooth pattern is a direct, visual read on library quality that has no real analog in RNA-seq QC. |
| Peak-level | **TSS enrichment score** | [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) / [`ataqv`](https://github.com/ParkerLab/ataqv) | Fold-enrichment of signal at annotated TSS vs. flanking background | ENCODE's standard cross-sample accessibility-signal quality bar — low enrichment usually means high background/noise even when mapping rate looks fine. |
| Peak-level | FRiP (Fraction of Reads in Peaks) | Peak caller output (e.g. [MACS2](https://github.com/macs3-project/MACS)) + `bedtools`/`samtools` | % of reads falling inside called peaks | Signal-to-noise proxy — a library can map well and still have most of its reads outside any real accessible region. |
| Peak-level | Peak count / width distribution | MACS2 (or similar) | Number of peaks called, their width distribution | Sanity-checks that peak calling itself behaved reasonably (e.g. thousands of near-genome-length "peaks" signals a parameter or input problem, not real biology). |
| Peak-level | ENCODE blacklist overlap | [ENCODE blacklist](https://github.com/Boyle-Lab/Blacklist) + `bedtools` | % of peaks/reads overlapping known artifact-prone regions | These regions (satellite repeats, assembly gaps) produce spuriously high signal in nearly every ATAC/ChIP-seq library regardless of biology — a standard exclusion filter, not a judgment about this specific sample. |
| Metadata cross-check | Synapse annotations vs. GEO's record | `geo-synapse` (same package the scRNA pipeline already uses) | Platform, reference genome, library prep kit, read length declared on Synapse vs. GEO's SOFT record | Same pattern as `scrnaseq_qc_pipeline/geo_metadata_qc.py` — curator-entered values checked against an independent external source, not re-validating the curation schema itself. |

## Notes for whenever this gets built

- **`fq_lint_qc.py` and `kraken_qc.py` are almost certainly directly reusable
  as-is** — nothing about FASTQ structural validation or Kraken2 species
  screening is scRNA-specific. `common/fastq_layout.py` likely is too, for the
  same reason (paired-end ATAC libraries are common).
- Most of the assay-specific value here is **post-alignment and post-peak-
  calling** — unlike the current scRNA pipeline, which only ever needs raw
  FASTQ + a processed count matrix, an ATAC-seq pipeline needs to actually
  align reads (or accept a pre-aligned BAM already tracked in Synapse) to
  compute almost everything that matters (mito fraction, fragment periodicity,
  TSS enrichment, FRiP). Scope that decision — align in-pipeline vs. only QC
  an already-aligned BAM someone else produced — before writing code.
- Carry over the scRNA pipeline's core design rule: **this pipeline's own code
  should never assert a PASS/FAIL/WARN verdict.** Only a named tool's own
  established criteria may (FastQC's own module thresholds, `fq lint`'s exit
  code); everything else — including ENCODE's own published QC thresholds for
  TSS enrichment/FRiP/mito fraction — should be reported as a plain number,
  since "what counts as good" varies by cell type, input amount, and protocol
  version in ways this pipeline has no basis to adjudicate on its own.
