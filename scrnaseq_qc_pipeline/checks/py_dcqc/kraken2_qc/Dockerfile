# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/py_dcqc/kraken2_qc/Dockerfile -t qc-kraken2 .
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     -v kraken2-db-cache:/db \
#     qc-kraken2 /bin/bash -euo pipefail -c \
#     "( (python3 /app/run.py --db-path /db 'sample.fastq.gz') > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
#
# KNOWN DEVIATION from py-dcqc's "What Makes a Good Test" rule ("Put all the
# reference data into the image. The tool must not download data when it
# runs."): run.py auto-downloads the ~5.5GB Kraken2 database on demand into
# --db-path rather than shipping it baked into the image, to avoid a
# multi-GB image build. Mount a *persistent* volume at --db-path (a named
# volume or host directory, never the container's own ephemeral
# filesystem) so that download only happens once, not on every run. Baking
# the database into the image is the fully-compliant alternative if this
# ever needs to be a literal py-dcqc submission.
FROM quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0
WORKDIR /app

# The biocontainers single-package image only ships kraken2's own runtime
# deps (perl/coreutils, no python) -- it's Debian-based, so apt-get is the
# simplest way to add python3 and procps (for Nextflow's task-metrics `ps`).
RUN apt-get update && apt-get install -y --no-install-recommends python3 procps && \
    rm -rf /var/lib/apt/lists/*

# nf-dcqc runs the container as the user who runs the workflow, not as root,
# and that user has no home directory in the image.
ENV HOME=/tmp
ENV PYTHONDONTWRITEBYTECODE=1

COPY qc_common.py kraken_qc.py /app/
COPY checks/py_dcqc/kraken2_qc/run.py /app/run.py
RUN chmod 0755 /app/run.py
