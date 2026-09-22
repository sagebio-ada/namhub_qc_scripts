# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/fastq_aggregate/Dockerfile -t qc-fastq-aggregate .
FROM python:3.14-slim
WORKDIR /app

RUN pip install --no-cache-dir requests

# fastq_qc.py imports fq_lint_qc/kraken_qc at module load time even though
# this check only calls aggregate_fastq_checks() -- both are copied so the
# import succeeds.
COPY qc_common.py fastq_qc.py fq_lint_qc.py kraken_qc.py /app/
COPY checks/fastq_aggregate/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
