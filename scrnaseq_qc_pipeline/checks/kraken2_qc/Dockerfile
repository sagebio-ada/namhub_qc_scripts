# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/kraken2_qc/Dockerfile -t qc-kraken2 .
# Runs on demand: --db-path is auto-downloaded (Kraken2's own Standard-8
# build, ~5.5GB) if nothing is already staged there. Mount a *persistent*
# volume at --db-path -- a named volume or host directory, not the
# container's own ephemeral filesystem -- so that download only happens
# once, e.g.:
#   docker volume create kraken2-db-cache
#   docker run --rm -v kraken2-db-cache:/db -v $(pwd)/work:/data qc-kraken2 \
#       --db-path /db --fastq-path /data/sample.fastq.gz ...
FROM quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0
WORKDIR /app

# The biocontainers single-package image only ships kraken2's own runtime
# deps (perl/coreutils, no python) -- it's Debian-based, so apt-get is the
# simplest way to add a python3 to run run.py's own report-parsing logic.
RUN apt-get update && apt-get install -y --no-install-recommends python3 && \
    rm -rf /var/lib/apt/lists/*

COPY qc_common.py kraken_qc.py /app/
COPY checks/kraken2_qc/run.py /app/run.py

ENTRYPOINT ["python3", "run.py"]
