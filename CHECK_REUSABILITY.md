# Check reusability across assay types

The checks in [`scrnaseq_qc_pipeline/checks/`](scrnaseq_qc_pipeline/checks/README.md)
were built for that pipeline's own raw-FASTQ and processed-matrix tiers, but
most of them contain no scRNA-seq-specific logic at all — they operate on
generic FASTQ files or generic registry/metadata lookups. This table assesses
whether each one could be reused as-is for ATAC-seq, whole-genome sequencing
(WGS), or bulk RNA-seq, before writing any assay-specific pipeline code (see
[`atacseq_qc_pipeline/README.md`](atacseq_qc_pipeline/README.md) and
[`bulkrnaseq_qc_pipeline/README.md`](bulkrnaseq_qc_pipeline/README.md) for
those assays' own candidate-check lists, independently derived).

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

Practical implication for a future `atacseq_qc_pipeline/` or
`bulkrnaseq_qc_pipeline/` implementation: start by importing/reusing these
10 checks directly rather than re-deriving raw-FASTQ QC from scratch, and
only write new assay-specific checks for the processed-data tier (peak
calling/BED-level QC for ATAC-seq, alignment/quantification-level QC for
bulk RNA-seq — see each pipeline's own README for those candidate lists).
