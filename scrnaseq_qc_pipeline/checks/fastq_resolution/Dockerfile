# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/fastq_resolution/Dockerfile -t qc-fastq-resolution .
FROM python:3.11-slim
WORKDIR /app

RUN pip install --no-cache-dir requests "git+https://github.com/sagebio-ada/geo_dataset_creation.git"

# fastq_qc.py imports fq_lint_qc/kraken_qc at module load time even though
# this check only calls resolve_ena_links() -- both are copied so the import
# succeeds; see the dockerize-checks skill for why modules are copied whole.
COPY qc_common.py fastq_qc.py fq_lint_qc.py kraken_qc.py /app/
COPY checks/fastq_resolution/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
