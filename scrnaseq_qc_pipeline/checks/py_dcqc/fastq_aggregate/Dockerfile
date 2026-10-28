# Build from the repo root (namhub_qc_scripts/), not scrnaseq_qc_pipeline/:
# this check's fq_lint_qc.py imports qc_common.fq_lint, and qc_common/ lives
# at the repo root, outside scrnaseq_qc_pipeline/ -- Docker can't COPY a path
# from outside the build context, so the context has to be the repo root.
#   docker build -f scrnaseq_qc_pipeline/checks/py_dcqc/fastq_aggregate/Dockerfile -t qc-fastq-aggregate .
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     qc-fastq-aggregate /bin/bash -euo pipefail -c \
#     "( (python3 /app/run.py --ena-metadata-json 'ena_metadata.json' 'r1.summary.json' 'r2.summary.json') > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
FROM python:3.14-slim
WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir requests

# nf-dcqc runs the container as the user who runs the workflow, not as root,
# and that user has no home directory in the image.
ENV HOME=/tmp
ENV PYTHONDONTWRITEBYTECODE=1

# fastq_qc.py imports fq_lint_qc/kraken_qc at module load time even though
# this check only calls aggregate_fastq_checks() -- both are copied so the
# import succeeds. fq_lint_qc.py in turn imports qc_common.fq_lint, so
# qc_common/ is copied alongside it (see fq_lint_qc.py's sys.path logic,
# which looks for a qc_common/ directory next to itself).
COPY scrnaseq_qc_pipeline/scrnaseq_common.py scrnaseq_qc_pipeline/fastq_qc.py scrnaseq_qc_pipeline/fq_lint_qc.py scrnaseq_qc_pipeline/kraken_qc.py /app/
COPY qc_common/ /app/qc_common/
COPY scrnaseq_qc_pipeline/checks/py_dcqc/fastq_aggregate/run.py /app/run.py
RUN chmod 0755 /app/run.py
