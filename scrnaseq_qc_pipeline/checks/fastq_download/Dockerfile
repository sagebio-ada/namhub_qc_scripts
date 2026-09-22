# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/fastq_download/Dockerfile -t qc-fastq-download .
FROM python:3.11-slim
WORKDIR /app

RUN pip install --no-cache-dir requests

# fastq_qc.py imports fq_lint_qc/kraken_qc at module load time even though
# this check only calls download_and_check_fastq() -- both are copied so the
# import succeeds; see the dockerize-checks skill for why modules are copied
# whole rather than further split.
COPY qc_common.py fastq_qc.py fq_lint_qc.py kraken_qc.py /app/
COPY checks/fastq_download/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
