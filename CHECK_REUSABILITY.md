# Check reusability across assay types

> **Note on repo state**: this document's first section was written while
> `atacseq_qc_pipeline/`/`bulkrnaseq_qc_pipeline/` were still README-only
> scaffolds (see the still-open PR #3). Since then, PR #5 merged a real,
> independently-built `atacseq_qc_pipeline/` and `rnaseq_qc_pipeline/` (note
> the name — not `bulkrnaseq_qc_pipeline/`) that take a different approach:
> summarizing already-processed nf-core/atacseq and nf-core/rnaseq output,
> not dockerized raw-FASTQ-level checks. The per-check table below still
> describes `scrnaseq_qc_pipeline/`'s own checks accurately; the "Universal
> vs. assay-specific" section further down reflects what's actually in the
> repo today, across both architectures.

The checks in [`scrnaseq_qc_pipeline/checks/`](scrnaseq_qc_pipeline/checks/README.md)
were built for that pipeline's own raw-FASTQ and processed-matrix tiers, but
most of them contain no scRNA-seq-specific logic at all — they operate on
generic FASTQ files or generic registry/metadata lookups. This table assesses
whether each one could be reused as-is for ATAC-seq, whole-genome sequencing
(WGS), or bulk RNA-seq, before writing any assay-specific pipeline code. The
**"Lives in py-dcqc?"** column is a separate question from assay reuse: it's
whether the check actually fits py-dcqc's own `ExternalBaseTest` contract
(3-way exit code, JSON metrics on stdout) well enough to belong there, versus
being this repo's own `repo_scripts/`-style check with its own
`--output findings.json` contract. The five raw-read-level checks fit it,
plus `matrix_content_qc` (a per-file structural test on the CellRanger matrix,
the same shape as the FASTQ-level tests just applied to a different file
type); everything else — fetch/resolve/download and GEO cross-checks — is
this repo's own orchestration logic, not a per-file structural test — see
**"Can py-dcqc express…"** below for why that boundary doesn't shift even
where a check is otherwise assay-agnostic.

| Check | Lives in py-dcqc? | ATAC-seq | WGS | scRNA-seq | Bulk RNA-seq | Why |
|---|---|---|---|---|---|---|
| `fq_lint` | ✓ | ✓ | ✓ | ✓ | ✓ | Pure FASTQ structural validation — doesn't know or care what the reads represent |
| `fastqc` | ✓ | ✓ | ✓ | ✓ | ✓ | Generic read-quality tool. Only the *interpretation* notes (10x hexamer bias, Nextera vs. TruSeq adapters) would need per-assay rewriting, not the check/container itself |
| `kraken2_qc` | ✓ | ✓ | ✓ | ✓ | ✓ | Species/contamination screening via k-mer matching — assay-blind by design |
| `fastq_aggregate` | ✓ | ✓ | ✓ | ✓ | ✓ | Compares FastQC output to an ENA registry record — works for any SRA/ENA-sourced FASTQ regardless of assay |
| `multiqc` | ✓ | ✓ | ✓ | ✓ | ✓ | Aggregates whatever FastQC/Kraken2 output exists — assay-blind |
| `synapse_fastq_download` | ✗ | ✓ | ✓ | ✓ | ✓ | Generic Synapse file fetch, no assay-specific logic |
| `fastq_download` | ✗ | ✓ | ✓ | ✓ | ✓ | Generic ENA URL fetch + size/md5 check |
| `fastq_resolution` | ✗ | ✓ | ✓ | ✓ | ✓ | Generic SRA→ENA link resolution via accession |
| `geo_metadata_fetch` | ✗ | ✓ | ~ | ✓ | ✓ | Works for any GEO-deposited series — but WGS is typically deposited to dbGaP/SRA, not GEO, so this wouldn't find a record for most WGS data |
| `geo_annotation_check` | ✗ | ✓ | ~ | ✓ | ✓ | Same caveat — depends on `geo_metadata_fetch` actually finding a GEO record |
| `matrix_content_qc` | ✓ | ✗ | ✗ | ✓ | ✗ | Hard-coded to CellRanger's specific 3-file sparse-matrix format (barcodes/features/matrix.mtx) — no equivalent concept in WGS or bulk RNA-seq's single-profile-per-sample output |
| `matrix_download` | ✗ | ✗ | ✗ | ✓ | ✗ | Same — downloads specifically the CellRanger 3-role file set |

**Totals: ATAC-seq 10/12, WGS 8/12 (+2 partial), scRNA-seq 12/12, Bulk RNA-seq 10/12. py-dcqc fit: 6/12.**

## Takeaway

**10 of the 12 checks are already assay-agnostic as written** — they operate
on raw FASTQ or generic file-fetch/metadata logic with zero scRNA-seq-specific
assumptions baked in. Only `matrix_content_qc`/`matrix_download` are genuinely
single-cell-specific (CellRanger's output format). `geo_metadata_fetch`/
`geo_annotation_check` are weaker for WGS specifically since that assay
usually lives in dbGaP/SRA rather than GEO — the same *pattern* (fetch
registry metadata, cross-check Synapse annotations against it) would just
need a dbGaP/SRA-biosample variant instead of a GEO one for WGS, not a
different architecture.

## Universal vs. assay-specific checks across the repo's actual pipelines

Pulling together everything now in this repo — `scrnaseq_qc_pipeline/`'s
dockerized checks, and PR #5's nf-core-summarizing `rnaseq_qc_pipeline/`/
`atacseq_qc_pipeline/` — the real picture splits into two tiers.

### Tier 1 — runs on every raw FASTQ, regardless of assay

| Check | Already built? |
|---|---|
| Structural validity (`fq lint`) | ✓ **genuinely shared** — `qc_common/fq_lint.py` runs and parses `fq lint`; `scrnaseq_qc_pipeline/checks/py_dcqc/fq_lint` and `qc_common/fastq_sections.fq_lint` are both thin wrappers over it |
| Read quality (FastQC) | ✓ both, but via different integration points, not shared code (`checks/py_dcqc/fastqc` runs FastQC itself per-file; `qc_common/fastq_sections.fastqc` reads MultiQC's already-aggregated tables from the nf-core pipeline's own upstream FastQC step) |
| Species/contamination (Kraken2 [+ Bracken]) | ✓ both, but not equivalent — `checks/py_dcqc/kraken2_qc` runs Kraken2 alone; `qc_common/fastq_sections.contaminants` runs Kraken2 **and** Bracken, re-estimating species abundance from Bracken's corrected counts (a real accuracy improvement the scrnaseq side lacks) |
| Aggregated report (MultiQC) | ✓ `scrnaseq_qc_pipeline/` only (`checks/py_dcqc/multiqc`) |
| SRA→ENA resolution / ENA registry size-md5-readcount check | ✓ `scrnaseq_qc_pipeline/` only (`checks/repo_scripts/fastq_resolution`, `fastq_download`, `checks/py_dcqc/fastq_aggregate`) |
| GEO metadata fetch + Synapse-annotation cross-check | ✓ `scrnaseq_qc_pipeline/` only — conditional: only applies when the series is actually GEO-deposited (true for scRNA-seq/bulk RNA-seq/ATAC-seq test datasets seen so far; **not** WGS, which is normally dbGaP/SRA) |

**fq lint has been consolidated** into one shared implementation
(`qc_common/fq_lint.py`); both architectures' own fq-lint checks are now
thin wrappers that each keep their own judgment logic (exit-code-based
PASS/FAIL for the py-dcqc contract vs. a TSV status column) on top of it.
Consolidating this also fixed a real bug on the nf-core side: its fq-lint
invocation never passed `--lint-mode log`, so it inherited `fq lint`'s
default "panic on first error" behavior — any file with more than one
structural issue only ever reported the first. The shared implementation
always uses log mode.

FastQC and Kraken2/Bracken remain genuinely divergent rather than merely
duplicated — FastQC because the two sides integrate with it at different
points (run-it-yourself vs. consume-nf-core's-own-output), and Kraken2
because the nf-core side's Bracken step is a real quality improvement, not
an equivalent reimplementation. Unifying either would mean changing one
side's behavior, not just refactoring — a decision to make deliberately,
not as a side effect of deduplication.

### Tier 2 — assay-specific, on top of Tier 1

| Assay | Checks | Already built? |
|---|---|---|
| **scRNA-seq** | CellRanger matrix QC (UMI/gene counts per cell, empty-barcode check, mito fraction); paired/interleaved/single-end layout detection (10x barcode+cDNA read structure) | ✓ (`checks/py_dcqc/matrix_content_qc`, `checks/repo_scripts/matrix_download`, `common/fastq_layout.py` on the still-open PR #2) |
| **Bulk RNA-seq** | STAR/Salmon mapping rate, Qualimap/Picard gene-body coverage & 3′ bias, strand-specificity inference, insert size, rRNA %, replicate concordance (Spearman) | ✓ (`rnaseq_qc_pipeline/rnaseq_qc_metrics.py`) |
| **ATAC-seq** | Mito fraction, NRF/PBC1/PBC2 complexity, fragment-size/nucleosome periodicity, TSS enrichment, FRiP, MACS2 peak calling, IDR replicate concordance, ENCODE blacklist | ✓ (`atacseq_qc_pipeline/atacseq_encode_qc_metrics.py`, `atacseq_idr_concordance.py`) |
| **WGS** | Mapping rate, coverage depth/uniformity (mosdepth/Picard `CollectWgsMetrics`), duplicate rate, insert size, Ti/Tv ratio, het/hom ratio, sex check from X/Y coverage, sample-contamination estimate (VerifyBamID) | ✗ not represented anywhere in this repo yet |

So: **Tier 1 exists twice, in two different styles, for three of four
assays; Tier 2 exists once per assay in whichever style that assay's
pipeline uses; WGS has nothing yet.**

### Tier 2, per check — `atacseq_qc_pipeline`/`rnaseq_qc_pipeline` (PR #5)

The same two-question analysis as the Tier 1 table above (reuse across
assays; fit with py-dcqc's contract), applied to the individual checks these
two pipelines actually implement — the rows the Tier 2 summary table above
collapses into one cell per assay. Checks already covered by the Tier 1
table (`fq_lint`, `fastqc`, contamination, read length) are not repeated
here even though both pipelines use them.

| Check | Lives in py-dcqc? | ATAC-seq | WGS | scRNA-seq | Bulk RNA-seq | Why |
|---|---|---|---|---|---|---|
| Strandedness (Salmon vs. RSeQC vs. declared) | ✗ | ✗ | ✗ | ~ | ✓ | Validates that a kit-dependent library-prep property (strand orientation) matches what was declared — ATAC/WGS have no such property; 10x scRNA kits use one fixed, already-known protocol, so there's nothing to infer and cross-check |
| Salmon fragment length | ✗ | ✗ | ✗ | ~ | ✓ | Salmon-quantification-specific; scRNA-seq is almost never quantified with Salmon/alevin in this repo's pipelines, it's CellRanger's own counting |
| STAR alignment (unique/multi/unmapped %) | ✗ | ~ | ~ | ~ | ✓ | The *concept* (mapped-read breakdown) is universal, but as built it's a parser for STAR's specific `Log.final.out` format; ATAC/WGS use different aligners (bwa) with different log formats, and CellRanger (which uses STAR internally) reports its own summary differently |
| Qualimap (exonic/intronic/intergenic %) | ✗ | ✗ | ✗ | ~ | ✓ | Needs a transcript/gene model to classify reads against; ATAC/WGS have no such framing. scRNA-seq could in principle run this on CellRanger's BAM, but exonic/intronic splitting is already handled differently there (intronic-inclusive counting modes) |
| featureCounts (assigned/no-feature/ambiguous/multimapping) | ✗ | ✗ | ✗ | ✗ | ✓ | Gene-level read assignment on pooled reads; scRNA-seq needs UMI-aware, per-cell-barcode counting (CellRanger/STARsolo's own job), not this bulk-style counter; ATAC/WGS have no expression-counting concept at all |
| rRNA % (biotype-based) | ✗ | ✗ | ✗ | ✓ | ✓ | rRNA-depletion efficiency is just as meaningful a QC signal for single-cell RNA-seq libraries as for bulk; ATAC/WGS profile DNA, not RNA, so there's no rRNA biotype to measure |
| Picard MarkDuplicates (dup rate, library size) | ✗ | ✓ | ✓ | ~ | ✓ | Generic BAM-level duplicate-rate QC, standard for any aligned DNA/RNA library; CellRanger's UMI-based dedup already resolves PCR duplicates before Picard would see the BAM, so flagging "duplicates" again there would double-count expected UMI collisions as a defect |
| Replicate concordance — Spearman on TPM | ✗ | ✗ | ✗ | ✗ | ✓ | TPM-correlation between replicates is specific to transcript/gene-level quantification; ATAC-seq already has its own ENCODE-correct replicate check for this (IDR on peaks, below) rather than a TPM analog; WGS has no expression values; scRNA-seq replicate concordance isn't assessed this way here |
| Alignment rate (raw / non-MT) [ENCODE] | ✗ | ✓ | ✓ | ~ | ~ | Mapped-fraction-from-a-BAM is a universal concept and standard WGS QC too; scRNA-seq's CellRanger already reports its own barcode-aware alignment rate, so re-deriving it from the pooled BAM the way this script does would lose that framing; bulk RNA-seq already gets the same information from STAR's own log, making this a redundant third measurement rather than something missing |
| Library complexity (NRF/PBC1/PBC2) [ENCODE] | ✗ | ✓ | ~ | ✗ | ✗ | ENCODE defines NRF/PBC specifically for DNA-based assays like ChIP/ATAC; WGS QC conventionally reports Picard's duplication rate instead (same underlying concern, different convention); RNA-seq (bulk or single-cell) is deliberately non-uniform — highly expressed transcripts are expected to produce many identical fragments, so a "low NRF" there is biology, not a library defect |
| Fragment size / nucleosome periodicity (NFR/mono) [ENCODE] | ✗ | ✓ | ✗ | ✗ | ✗ | This is specifically the signal that distinguishes ATAC-seq (periodicity from DNA wrapped around nucleosomes); none of the other three assays' fragment-length distributions carry that signal |
| Peak enrichment (FRiP, ataqv TSS enrichment, HQAA-overlap %) | ✗ | ✓ | ✗ | ✗ | ✗ | All three require a called-peak concept (MACS2); none of WGS, scRNA-seq or bulk RNA-seq call peaks |
| IDR replicate concordance | ✗ | ✓ | ✗ | ✗ | ✗ | Same peak-calling dependency as above, plus needs ≥2 replicate samples to compare against each other — the clearest case in this whole table of something py-dcqc's single-file, single-`FileType` `Suite` model structurally cannot express, regardless of the verdict question (see below) |

**py-dcqc fit: 0/14.** Every check here is either (a) a reported-only metric
with no tool-native pass/fail — consistent with this repo's own "no
self-asserted verdicts" rule, but a mismatch for py-dcqc's exit-code-based
`ExternalBaseTest` contract regardless of which way that rule gets applied —
or (b) inherently cross-file (joins a run row to a library row, compares
replicate samples to each other, needs a samplesheet's declared value to
check against) rather than a single-file structural test, which is the one
shape py-dcqc's `FileType`/`Suite` resolution can express at all (see
**"Can py-dcqc express…"** below). The one near-exception is
`idr_reproducibility`, which *is* a genuine tool-native verdict
(pass/borderline/fail, ENCODE's own thresholds) — but IDR is also the most
clearly cross-sample check in the table (it doesn't exist without ≥2
replicates to compare), so it fails on structural grounds even though it
would otherwise pass the verdict test that ruled out everything else.

This is a different conclusion from the Tier 1 table above: there, 6 of 12
checks fit py-dcqc because they're genuinely single-file structural tests on
a raw FASTQ. Here, zero of 14 fit, because every one of these checks exists
specifically *because* it needs more context than one file — an aligned BAM
plus a samplesheet's declared value, or one replicate's peaks plus another's.
That's not a gap to close; it's what Tier 2 checks are for.

## Can py-dcqc express "universal vs. assay-specific" today?

No — checked directly against py-dcqc's current source
(Sage-Bionetworks-Workflows/py-dcqc, verified by reading `src/dcqc/file.py`
and `src/dcqc/suites/suites.py`/`suite_abc.py`), not assumed. There is
exactly **one** generic `FASTQ` `FileType` and exactly **one** `FastqSuite`
(`add_tests = (tests.PairedFastqParityTest,)`, inherited from the base
`FileSuite`'s `FileExtensionTest`/`Md5ChecksumTest`) applied identically to
every FASTQ file regardless of assay. `SuiteABC.get_subclass_by_file_type()`
resolves a Suite strictly by `FileType.name` — a flat, collision-free
registry lookup with no assay concept anywhere: zero mentions of "assay" in
py-dcqc's entire codebase, docs, or issue/PR history (checked directly).

`File.metadata` (arbitrary extra CSV columns passed through from the input
manifest) already exists as an unused extensibility point — an `assay` key
could be stashed there today with zero effect, since nothing in Suite
resolution reads it. Making py-dcqc assay-aware would need new code along
one of two lines, neither of which exists today or is discussed anywhere in
py-dcqc's issue tracker:
1. New assay-specific `FileType` identifiers (e.g. a distinct type per
   assay) plus corresponding `Suite` subclasses, or
2. Extending `SuiteABC`'s resolution logic (and/or `File`/`FileType`) to
   also consult `File.metadata` when choosing which tests to bundle for a
   file, rather than `FileType` alone.

Until one of those exists, universal-vs-assay-specific scoping is something
this repo has to manage on its own — which is exactly what the two tables
above do.
