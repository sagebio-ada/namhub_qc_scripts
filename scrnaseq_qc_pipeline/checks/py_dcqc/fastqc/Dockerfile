# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/py_dcqc/fastqc/Dockerfile -t qc-fastqc .
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     qc-fastqc /bin/bash -euo pipefail -c \
#     "( (python3 /app/run.py 'sample.fastq.gz') > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
FROM quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
WORKDIR /app

# The biocontainers single-package image only ships FastQC's own runtime
# deps (openjdk/perl, no python) -- it's Debian-based, so apt-get is the
# simplest way to add python3 and procps (for Nextflow's task-metrics `ps`).
# python3-requests is needed because fastq_qc.py (copied whole -- see the
# dockerize-checks skill) imports requests at module load time even though
# this check doesn't call the download functions.
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-requests procps && \
    rm -rf /var/lib/apt/lists/*

# nf-dcqc runs the container as the user who runs the workflow, not as root,
# and that user has no home directory in the image.
ENV HOME=/tmp
ENV PYTHONDONTWRITEBYTECODE=1

COPY qc_common.py fastq_qc.py fq_lint_qc.py kraken_qc.py /app/
COPY checks/py_dcqc/fastqc/run.py /app/run.py
RUN chmod 0755 /app/run.py
