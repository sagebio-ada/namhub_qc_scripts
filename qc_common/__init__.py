"""Code shared by the QC pipelines in namhub_qc_scripts (standard library only).

util            to_num, pct, read_table
samplesheet     read sample names from an nf-core samplesheet
runner          run a list of metric "sections" per sample and write the TSV
fastq_sections  sections for FASTQ-based data (fq-lint, read length, FastQC, kraken2/Bracken)
tools           run fq lint and kraken2+bracken (writes where fastq_sections reads)
fastqc_zip      median read length from a FastQC zip

Pipeline scripts import this with:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from qc_common import ...
"""
