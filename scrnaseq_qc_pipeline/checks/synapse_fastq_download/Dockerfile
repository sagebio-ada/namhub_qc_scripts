# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/synapse_fastq_download/Dockerfile -t qc-synapse-fastq-download .
# Run with SYNAPSE_AUTH_TOKEN set in the environment (docker run -e SYNAPSE_AUTH_TOKEN ...).
FROM python:3.11-slim
WORKDIR /app

RUN pip install --no-cache-dir synapseclient

COPY qc_common.py fastq_qc.py fq_lint_qc.py kraken_qc.py /app/
COPY checks/synapse_fastq_download/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
