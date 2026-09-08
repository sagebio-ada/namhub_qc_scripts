# namhub_qc_scripts

QC scripts for NamHub/ADA data pipelines.

## Pipelines

- [`scrnaseq_qc_pipeline/`](scrnaseq_qc_pipeline/README.md) — QC pipeline for single-cell RNA-seq files under a Synapse project/folder: cross-checks curated Synapse annotations against GEO's own record, verifies downloaded raw FASTQ data against ENA's file report, runs FastQC/`fq lint`/Kraken2 content QC, and sanity-checks CellRanger processed matrices.
