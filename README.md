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
