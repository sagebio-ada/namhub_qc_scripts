# namhub_qc_scripts

Quality-control scripts for biomedical data pipelines that track datasets in
[Synapse](https://www.synapse.org), a data-sharing platform used by several
research consortia. Each pipeline here is a working example of what QC for
its data type can look like — see each pipeline's own README for background,
usage, and a full breakdown of what it checks.

## Pipelines

- [`scrnaseq_qc_pipeline/`](scrnaseq_qc_pipeline/README.md) — QC for
  single-cell RNA-seq data tracked in Synapse but originating from a public
  archive (GEO/SRA): cross-checks curated Synapse annotations against GEO's
  own record, verifies downloaded raw sequencing data against the archive's
  own record, runs read-quality/contamination screening, and sanity-checks
  processed (CellRanger) output matrices.

- [`rnaseq_qc_pipeline/`](rnaseq_qc_pipeline/README.md) — QC summary (one TSV row
  per sample) for nf-core/rnaseq output.
- [`atacseq_qc_pipeline/`](atacseq_qc_pipeline/README.md) — QC for nf-core/atacseq
  output (combined metrics table, ENCODE-style metrics, IDR concordance).

See [`CHECK_REUSABILITY.md`](CHECK_REUSABILITY.md) for which of
`scrnaseq_qc_pipeline/`'s individual checks are assay-agnostic enough to
reuse directly for ATAC-seq, whole-genome sequencing, and bulk RNA-seq —
independent from, and worth comparing against, `rnaseq_qc_pipeline/` and
`atacseq_qc_pipeline/` above, which take a different approach (summarizing
already-processed nf-core/rnaseq and nf-core/atacseq output rather than
dockerized raw-FASTQ-level checks).

## Shared code

[`qc_common/`](qc_common/__init__.py) holds code meant to be reused across
pipelines (standard library only):

- `util.py` — `to_num`, `pct`, `read_table`
- `samplesheet.py` — read sample names from an nf-core samplesheet
- `runner.py` — run a list of metric "sections" per sample (a missing or
  unparsable file blanks that section and warns) and write the TSV
- `fastq_sections.py` — sections for FASTQ-based data: fq-lint, read length,
  FastQC, kraken2/Bracken contamination
- `tools.py` — run fq lint and kraken2 + Bracken
- `fastqc_zip.py` — median read length from a FastQC zip

A pipeline script adds the repo root to `sys.path` and imports from
`qc_common` (see `rnaseq_qc_pipeline/rnaseq_qc_metrics.py`). Assay-specific
sections stay in the pipeline's own folder. Pipelines sit side by side, and
code is moved into `qc_common` when a second pipeline needs it.

Note: an unmerged branch for `scrnaseq_qc_pipeline/` independently started a
similarly-named shared package (`common/`, holding a `fastq_layout.py`
helper) — a naming collision worth resolving if/when these two efforts
converge.
