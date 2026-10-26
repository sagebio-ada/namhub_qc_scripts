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
(WGS), or bulk RNA-seq, before writing any assay-specific pipeline code.

| Check | ATAC-seq | WGS | scRNA-seq | Bulk RNA-seq | Why |
|---|---|---|---|---|---|
| `fq_lint` | ✓ | ✓ | ✓ | ✓ | Pure FASTQ structural validation — doesn't know or care what the reads represent |
| `fastqc` | ✓ | ✓ | ✓ | ✓ | Generic read-quality tool. Only the *interpretation* notes (10x hexamer bias, Nextera vs. TruSeq adapters) would need per-assay rewriting, not the check/container itself |
| `kraken2_qc` | ✓ | ✓ | ✓ | ✓ | Species/contamination screening via k-mer matching — assay-blind by design |
| `fastq_aggregate` | ✓ | ✓ | ✓ | ✓ | Compares FastQC output to an ENA registry record — works for any SRA/ENA-sourced FASTQ regardless of assay |
| `multiqc` | ✓ | ✓ | ✓ | ✓ | Aggregates whatever FastQC/Kraken2 output exists — assay-blind |
| `synapse_fastq_download` | ✓ | ✓ | ✓ | ✓ | Generic Synapse file fetch, no assay-specific logic |
| `fastq_download` | ✓ | ✓ | ✓ | ✓ | Generic ENA URL fetch + size/md5 check |
| `fastq_resolution` | ✓ | ✓ | ✓ | ✓ | Generic SRA→ENA link resolution via accession |
| `geo_metadata_fetch` | ✓ | ~ | ✓ | ✓ | Works for any GEO-deposited series — but WGS is typically deposited to dbGaP/SRA, not GEO, so this wouldn't find a record for most WGS data |
| `geo_annotation_check` | ✓ | ~ | ✓ | ✓ | Same caveat — depends on `geo_metadata_fetch` actually finding a GEO record |
| `matrix_content_qc` | ✗ | ✗ | ✓ | ✗ | Hard-coded to CellRanger's specific 3-file sparse-matrix format (barcodes/features/matrix.mtx) — no equivalent concept in WGS or bulk RNA-seq's single-profile-per-sample output |
| `matrix_download` | ✗ | ✗ | ✓ | ✗ | Same — downloads specifically the CellRanger 3-role file set |

**Totals: ATAC-seq 10/12, WGS 8/12 (+2 partial), scRNA-seq 12/12, Bulk RNA-seq 10/12.**

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
| Structural validity (`fq lint`) | ✓ both architectures (`scrnaseq_qc_pipeline/checks/py_dcqc/fq_lint`, `qc_common/fastq_sections.fq_lint`) |
| Read quality (FastQC) | ✓ both (`checks/py_dcqc/fastqc`, `qc_common/fastq_sections.fastqc`) |
| Species/contamination (Kraken2 [+ Bracken]) | ✓ both (`checks/py_dcqc/kraken2_qc`, `qc_common/fastq_sections.contaminants`) |
| Aggregated report (MultiQC) | ✓ `scrnaseq_qc_pipeline/` only (`checks/py_dcqc/multiqc`) |
| SRA→ENA resolution / ENA registry size-md5-readcount check | ✓ `scrnaseq_qc_pipeline/` only (`checks/repo_scripts/fastq_resolution`, `fastq_download`, `checks/py_dcqc/fastq_aggregate`) |
| GEO metadata fetch + Synapse-annotation cross-check | ✓ `scrnaseq_qc_pipeline/` only — conditional: only applies when the series is actually GEO-deposited (true for scRNA-seq/bulk RNA-seq/ATAC-seq test datasets seen so far; **not** WGS, which is normally dbGaP/SRA) |

**This tier is already duplicated between the two architectures** — one
dockerized, py-dcqc-contract-aligned implementation in
`scrnaseq_qc_pipeline/checks/py_dcqc/`, and one plain-Python,
nf-core-output-summarizing implementation in `qc_common/fastq_sections.py`.
Worth reconciling into one shared implementation before either grows
further, rather than letting both evolve independently.

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
