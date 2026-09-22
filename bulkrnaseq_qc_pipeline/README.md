# Bulk RNA-seq QC Pipeline (planned)

**Status: planning only — no code here yet.** This README is a reference list of
QC checks worth running for bulk RNA-seq data tracked in Synapse, in the same
spirit as [`scrnaseq_qc_pipeline/`](../scrnaseq_qc_pipeline/README.md):
candidates to evaluate and pick from when this pipeline actually gets built,
not a committed scope.

## Background, for readers new to this assay

Bulk RNA-seq sequences the pooled mRNA (or total RNA) of a whole sample —
unlike single-cell RNA-seq, there's no per-cell barcode/UMI structure and no
raw-vs-filtered matrix distinction; the end product is one gene/transcript
expression profile per sample. That has a few direct QC consequences:

- **No cell-barcode/UMI layer to inspect** — bulk libraries are usually
  standard paired- or single-end fragments, so raw-read QC looks much more
  like generic short-read sequencing QC than the scRNA pipeline's chemistry-
  aware FastQC interpretation notes.
- **Ribosomal RNA depletion (or poly-A selection) is a real failure point.**
  Most RNA-seq protocols either poly-A-select for mRNA or deplete ribosomal
  RNA directly — a failed depletion/selection step shows up as a large
  fraction of reads mapping to rRNA loci, which silently eats sequencing
  depth from everything else.
- **Where reads land across a gene's length and across genomic features
  (exon/intron/intergenic) is diagnostic**, not just how many reads there
  are — 3' bias signals degraded RNA, a high intronic/intergenic fraction
  signals genomic DNA contamination or immature-transcript capture, and the
  library's declared strandedness protocol can be verified directly against
  the data rather than trusted at face value.
- **QC has both a per-sample tier and a cohort tier** — sample correlation,
  PCA outlier detection, and batch-effect screening only make sense across a
  set of samples that are supposed to be comparable, not one file in
  isolation the way most of the scRNA pipeline's checks are.

## Possible checks

| Tier | Check | Tool | What it would report | Why it matters for bulk RNA-seq specifically |
|---|---|---|---|---|
| Raw FASTQ | Structural validity | [`fq lint`](https://github.com/stjude-rust-labs/fq) | Malformed records, mismatched seq/quality lengths | Same generic file-integrity check as the scRNA pipeline — not assay-specific. |
| Raw FASTQ | Read quality / adapter content | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Per-base quality, adapter read-through, duplication, GC content | Standard TruSeq-style adapters for most bulk protocols (vs. scRNA's 10x-specific artifacts) — the scRNA pipeline's "this WARN is expected for 10x chemistry" interpretation notes mostly don't transfer and would need re-deriving for whichever prep kit(s) this covers. |
| Raw FASTQ | Species/contamination screening | [Kraken2](https://github.com/DerrickWood/kraken2) | % reads classified per organism | Same tool `scrnaseq_qc_pipeline/kraken_qc.py` already wraps — directly reusable here. |
| Raw FASTQ | **rRNA contamination** | [SortMeRNA](https://github.com/sortmerna/sortmerna) or a lightweight k-mer screen against an rRNA reference | % of reads matching ribosomal RNA sequence | A failed poly-A selection or rRNA depletion step is one of the most common bulk RNA-seq library-prep failures, and it's invisible to FastQC/Kraken2 — a dedicated check catches it directly. |
| Alignment | Mapping rate, unique vs. multi-mapped | Aligner's own log (e.g. STAR's `Log.final.out`), `samtools flagstat` | % aligned, % uniquely aligned, % too-short/unmapped | Baseline signal that anything downstream is even trustworthy — RNA-seq aligners (STAR/HISAT2) report this natively, no extra tooling needed. |
| Alignment | Read distribution across genomic features | [RSeQC](http://rseqc.sourceforge.net/) `read_distribution.py` or Picard `CollectRnaSeqMetrics` | % reads in exons / introns / UTRs / intergenic | A high intronic or intergenic fraction points at genomic DNA contamination or substantial pre-mRNA capture, not a sequencing problem — a distinction only this kind of check can make. |
| Alignment | **5'-3' gene body coverage bias** | RSeQC `geneBody_coverage.py` or Picard `CollectRnaSeqMetrics` | Coverage profile across normalized transcript length | RNA degradation shows up as 3' (or 5', protocol-dependent) coverage bias well before it's visible in any raw-read metric — this is the closest bulk-RNA-seq analog to the scRNA pipeline's mitochondrial-fraction check: a per-cell/per-sample data-quality signal, not a metadata comparison. |
| Alignment | Strand specificity vs. declared protocol | RSeQC `infer_experiment.py` | Inferred library strandedness (unstranded / fwd / reverse) | Directly verifiable against the data, unlike most declared metadata — a mismatch between what's declared and what's inferred here is a genuine, checkable discrepancy (same spirit as the scRNA pipeline's GEO cross-checks, but derived from the data itself rather than a second external source). |
| Alignment | Insert size distribution | Picard `CollectInsertSizeMetrics` | Fragment size histogram (paired-end only) | Flags library-prep size-selection problems; also informs whether the declared fragment size in any metadata is plausible. |
| Alignment | Duplicate rate | Picard `MarkDuplicates` | % duplicate reads | Same caveat the scRNA pipeline already documents for its own duplication check: highly-expressed transcripts are legitimately sequenced many times over, so a high rate isn't automatically bad — report the number, not a verdict. |
| Quantification | Genes/transcripts detected | Quantifier output (e.g. [Salmon](https://github.com/COMBINE-lab/salmon), `featureCounts`, StringTie) | Count of genes/transcripts with nonzero expression | A basic complexity/depth sanity check on the final expression profile, independent of alignment-stage metrics. |
| Quantification | Sequencing depth / saturation | Quantifier output, or a downsampling curve | Whether additional depth would still be finding new expressed genes | Tells you whether a sample was sequenced deeply enough for its intended use (e.g. differential expression vs. just detecting major transcripts). |
| Cohort-level | Sample correlation / PCA outlier detection | Custom code (`numpy`/`scipy`) over a cohort's expression matrices | Pairwise sample correlation, PCA outlier flags | Only meaningful across a set of samples meant to be comparable — an outlier here (a sample that doesn't cluster with its expected group) is a strong, assay-generic quality signal, but requires a cohort-level entry point this pipeline's per-file design doesn't currently have. |
| Metadata cross-check | Synapse annotations vs. GEO's record | `geo-synapse` (same package the scRNA pipeline already uses) | Platform, reference genome/transcriptome build, library prep kit, declared strandedness vs. GEO's SOFT record | Same pattern as `scrnaseq_qc_pipeline/geo_metadata_qc.py` — curator-entered values checked against an independent external source. |

## Notes for whenever this gets built

- **`fq_lint_qc.py` and `kraken_qc.py` are almost certainly directly reusable
  as-is** — same reasoning as the ATAC-seq README: nothing about FASTQ
  structural validation or Kraken2 species screening is assay-specific.
  `common/fastq_layout.py` applies here too (bulk RNA-seq is commonly
  paired-end).
- Unlike the current scRNA pipeline (which only needs a raw FASTQ + a
  processed count matrix per sample), most of the assay-specific value here
  needs an **alignment step** (read-distribution, gene-body coverage, strand
  inference, insert size) — scope whether this pipeline aligns reads itself
  or only QCs an already-aligned BAM someone else produced, same open
  question as ATAC-seq.
- The **cohort-level tier is a genuinely different shape of check** than
  everything else in this repo so far — every existing check (scRNA included)
  evaluates one file/sample independently. Sample correlation/PCA outlier
  detection needs multiple samples' results at once, which doesn't fit the
  current per-entity `Finding` model without some extension (e.g. a
  cohort-level check that runs after every individual sample's checks and
  takes the whole group as input).
- Carry over the scRNA pipeline's core design rule: **this pipeline's own code
  should never assert a PASS/FAIL/WARN verdict.** Only a named tool's own
  established criteria may (FastQC's own module thresholds, `fq lint`'s exit
  code); everything else — read-distribution fractions, coverage bias,
  duplicate rate, rRNA %, cohort outlier flags — should be reported as a plain
  number, since "what counts as good" depends on tissue type, protocol, and
  study design in ways this pipeline has no basis to adjudicate on its own.
