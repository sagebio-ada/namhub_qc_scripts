# Build from the repo root (namhub_qc_scripts/), not scrnaseq_qc_pipeline/:
# this check's fq_lint_qc.py imports qc_common.fq_lint, and qc_common/ lives
# at the repo root, outside scrnaseq_qc_pipeline/ -- Docker can't COPY a path
# from outside the build context, so the context has to be the repo root.
#   docker build -f scrnaseq_qc_pipeline/checks/repo_scripts/synapse_fastq_download/Dockerfile -t qc-synapse-fastq-download .
# Run with SYNAPSE_AUTH_TOKEN set in the environment (docker run -e SYNAPSE_AUTH_TOKEN ...).
FROM python:3.14-slim
WORKDIR /app

RUN pip install --no-cache-dir synapseclient

# fq_lint_qc.py imports qc_common.fq_lint, so qc_common/ is copied alongside
# it (see fq_lint_qc.py's sys.path logic, which looks for a qc_common/
# directory next to itself).
COPY scrnaseq_qc_pipeline/scrnaseq_common.py scrnaseq_qc_pipeline/fastq_qc.py scrnaseq_qc_pipeline/fq_lint_qc.py scrnaseq_qc_pipeline/kraken_qc.py /app/
COPY qc_common/ /app/qc_common/
COPY scrnaseq_qc_pipeline/checks/repo_scripts/synapse_fastq_download/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
