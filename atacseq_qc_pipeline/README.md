# ATAC-seq QC Pipeline (planned)

**Status: planning only — no code here yet.** This README is a reference list of
QC checks worth running for ATAC-seq data tracked in Synapse, in the same spirit
as [`scrnaseq_qc_pipeline/`](../scrnaseq_qc_pipeline/README.md): candidates to
evaluate and pick from when this pipeline actually gets built, not a committed
scope.

## NAMHub planning context

Per NAMHub's own assay-tracking table: ATAC-seq submissions so far are
**externally-deposited from GEO** (e.g. GSE195688) as **Bed, BedNarrowPeak, and
Tar** files from site BWH, **not reprocessed by NAMHub**, with an earliest due
date of 11/15/2026. No tool or pipeline is currently named for ATAC-seq in
either the NYU planning docs or the DCQC doc — this README is the first pass
at filling that gap.

**This matters for scope**: Bed/BedNarrowPeak means the submissions seen *so
far* arrive already **peak-called** — but raw FASTQ for this assay is still a
real possibility going forward (per internal discussion, not something to
assume away). So all three tiers below stay in scope: the peak-level tier is
what's directly runnable against today's actual BedNarrowPeak/Tar
submissions with no realignment needed, while the raw-FASTQ/alignment tiers
should be built out too so this pipeline is ready the moment raw data starts
arriving, not scrambling to add them after the fact. Worth confirming what's
actually bundled inside the `Tar` file (just BED outputs, or upstream
BAMs/FASTQs/logs too) — it may already contain more than the peak calls.

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

| Tier | Check | Candidate tools | What it would report | Why it matters for ATAC-seq specifically |
|---|---|---|---|---|
| Raw FASTQ¹ | Structural validity | [`fq lint`](https://github.com/stjude-rust-labs/fq) | Malformed records, mismatched seq/quality lengths | Same generic file-integrity check as the scRNA pipeline — not assay-specific. |
| Raw FASTQ¹ | Read quality / adapter content | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Per-base quality, adapter read-through, duplication, GC content | ATAC libraries use **Nextera** adapters (not TruSeq) — configure FastQC's adapter list accordingly or its Adapter Content module won't flag real contamination. |
| Raw FASTQ¹ | Species/contamination screening | [Kraken2](https://github.com/DerrickWood/kraken2) | % reads classified per organism | Same tool `scrnaseq_qc_pipeline/kraken_qc.py` already wraps — directly reusable here, not assay-specific logic. |
| Alignment¹ | Mapping rate, unique vs. multi-mapped | `samtools flagstat`; [Bowtie2](https://github.com/BenLangmead/bowtie2) or `bwa mem`'s own log | % aligned, % uniquely aligned | A low unique-mapping rate flags degraded/contaminated libraries before any accessibility signal is even computable. |
| Alignment¹ | **Mitochondrial read fraction** | `samtools idxstats`; [`ataqv`](https://github.com/ParkerLab/ataqv) | % of reads mapping to chrM | *The* headline ATAC-seq QC number — extremely high values indicate poor nuclei isolation (whole-cell contamination) diluting the accessible-chromatin signal with mitochondrial noise. |
| Alignment¹ | PCR duplicate rate | Picard `MarkDuplicates`; `samtools markdup`; `sambamba markdup` | % duplicate reads | Tagmentation + PCR amplification makes ATAC libraries prone to duplication at low input; distinguishing PCR duplicates from genuine independent tagmentation events at the same site matters more here than in most assays. |
| Alignment¹ | Library complexity (NRF, PBC1, PBC2) | [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) scripts; `ataqv` | Non-redundant fraction, PCR bottlenecking coefficients | Standard ENCODE complexity metrics — flags a library that's been PCR-over-amplified from too little starting material, a common ATAC-seq failure mode. |
| Alignment¹ | **Fragment size distribution** | Picard `CollectInsertSizeMetrics`; `ataqv`; custom (`numpy` over BAM insert sizes) | Insert-size histogram; presence/absence of nucleosome-spacing periodicity | ATAC-seq's signature plot — the nucleosome-free/mono-/di-nucleosome sawtooth pattern is a direct, visual read on library quality that has no real analog in RNA-seq QC. |
| Peak-level | **TSS enrichment score** | [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline); `ataqv`; [deepTools](https://deeptools.readthedocs.io/) `computeMatrix`/`plotProfile` + custom enrichment calc | Fold-enrichment of signal at annotated TSS vs. flanking background | ENCODE's standard cross-sample accessibility-signal quality bar — low enrichment usually means high background/noise even when mapping rate looks fine. |
| Peak-level | FRiP (Fraction of Reads in Peaks) | `bedtools intersect` + `samtools`; `ataqv`; ENCODE ATAC-seq pipeline | % of reads falling inside called peaks | Signal-to-noise proxy — a library can map well and still have most of its reads outside any real accessible region. |
| Peak-level | Peak count / width distribution | Custom code (`pandas`/`numpy` over the BED/narrowPeak file) | Number of peaks called, their width distribution | Sanity-checks that peak calling itself behaved reasonably (e.g. thousands of near-genome-length "peaks" signals a parameter or input problem, not real biology) — directly runnable on NAMHub's actual BedNarrowPeak submissions with no realignment needed. |
| Peak-level | ENCODE blacklist overlap | `bedtools intersect` + [ENCODE blacklist](https://github.com/Boyle-Lab/Blacklist) BED file | % of peaks overlapping known artifact-prone regions | These regions (satellite repeats, assembly gaps) produce spuriously high signal in nearly every ATAC/ChIP-seq library regardless of biology — a standard exclusion filter, not a judgment about this specific sample. Also directly runnable on a submitted BED/narrowPeak file alone. |
| Peak-calling (once raw data arrives) | Peak calling itself | [MACS2](https://github.com/macs3-project/MACS); [Genrich](https://github.com/jsh58/Genrich) (ATAC-aware mode, no shifting-model step required) | Called peaks, from a BAM | Becomes relevant as soon as NAMHub receives pre-peak-calling data (BAM/FASTQ) for this assay — worth having a candidate picked out now rather than deciding under time pressure once that data shows up. |
| Metadata cross-check | Synapse annotations vs. GEO's record | `geo-synapse` (same package the scRNA pipeline already uses) | Platform, reference genome, library prep kit, read length declared on Synapse vs. GEO's SOFT record | Same pattern as `scrnaseq_qc_pipeline/geo_metadata_qc.py` — curator-entered values checked against an independent external source, not re-validating the curation schema itself. Applicable regardless of which file type NAMHub actually receives. |
| File format/content validation | OME/format metadata validation | [`dcqc`](https://github.com/Sage-Bionetworks-Workflows/py-dcqc) | File-type-appropriate structural/metadata validation | NAMHub's own general-purpose file-content-validation tool (used across several other assays in the broader tracking table) — worth checking whether it already covers BED/narrowPeak/Tar before building anything custom here. |

¹ Needs raw FASTQ or an intermediate BAM for this assay — today's actual
submissions are already peak-called, but raw FASTQ is expected to be a real
possibility for this assay going forward (see "NAMHub planning context"
above), so these tiers are worth building out rather than deferring.

## Notes for whenever this gets built

- **`fq_lint_qc.py` and `kraken_qc.py` are almost certainly directly reusable
  as-is** once raw FASTQ starts arriving for this assay — nothing about FASTQ
  structural validation or Kraken2 species screening is scRNA-specific.
  `common/fastq_layout.py` likely is too, for the same reason (paired-end ATAC
  libraries are common).
- **Check `dcqc` first** before building custom peak-level checks — it's
  already NAMHub's named tool for file-content/format validation across other
  assays (OME-TIFF metadata for imaging, for instance); it may already have or
  be extensible to a BED/narrowPeak validator, which would cover some of the
  peak-level tier above without new code.
- Given raw FASTQ is expected but not here yet, a reasonable build order is
  **peak-level tier first** (directly usable against today's real
  BedNarrowPeak/Tar submissions) with the **raw-FASTQ/alignment tiers built
  out in parallel or immediately after**, so this pipeline doesn't lag behind
  once raw data actually starts showing up.
- Carry over the scRNA pipeline's core design rule: **this pipeline's own code
  should never assert a PASS/FAIL/WARN verdict.** Only a named tool's own
  established criteria may (FastQC's own module thresholds, `fq lint`'s exit
  code); everything else — including ENCODE's own published QC thresholds for
  TSS enrichment/FRiP/mito fraction — should be reported as a plain number,
  since "what counts as good" varies by cell type, input amount, and protocol
  version in ways this pipeline has no basis to adjudicate on its own. This
  also matches where NAMHub's own planning currently stands org-wide:
  thresholds for flagging low-quality files are explicitly called out as
  not-yet-set for every assay in the tracking table, not just this one.
