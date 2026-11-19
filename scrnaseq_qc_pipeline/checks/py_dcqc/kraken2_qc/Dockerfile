# Build from the repo root (namhub_qc_scripts/):
#   docker build -f scrnaseq_qc_pipeline/checks/py_dcqc/kraken2_qc/Dockerfile -t qc-kraken2-qc .
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     qc-kraken2-qc /bin/bash -euo pipefail -c \
#     "( (python3 /app/run.py 'r1.fastq.gz' 'r2.fastq.gz') > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
#
# Kraken2 and Python from one micromamba environment -- one package
# manager, no files copied between base images.
# Pinned for reproducible rebuilds: the base image by digest, and every
# conda package (with its dependencies) by exact build and md5 in
# conda-linux-64.lock. apt packages come from Debian stable, which only
# takes security fixes; they aren't version-pinned because Debian removes
# superseded versions, which would break the build.
#
# To change a package: build without the lock and regenerate it with
#   docker run --rm --entrypoint micromamba IMAGE env export -n base --explicit --md5
# Spec the lock was made from: python=3.10 kraken2=2.1.3
#   (channels: conda-forge, bioconda)
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.10-bookworm-slim@sha256:e3797091302382ea841498bc93a7b0a50f7c1448333d5e946d2d1608d0c5f43d

USER root
# procps: Nextflow's task metrics run `ps` inside the container. curl:
# downloads the Kraken2 database below.
RUN apt-get update && apt-get install -y --no-install-recommends procps curl ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# py-dcqc tests must not download data when they run, so Kraken2's
# "Standard-8" index (~6GB compressed) is baked in at build time. This layer
# depends only on the two ARGs below, not on any code in this repo, so
# editing the Python files or the conda lock never forces a re-download
# (it comes before both).
# The index site publishes no checksums, so KRAKEN2_DB_SHA256 is the
# SHA-256 of the tarball as downloaded on 2026-11-16 (two separate
# downloads agreed); the build fails if the file at that URL ever changes.
# Update both together.
ARG KRAKEN2_DB_URL=https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
ARG KRAKEN2_DB_SHA256=e592faea3e307f01b40d92e402a1707cefbfdd3f2b8acae886dd487eea5c73de
# The tarball is streamed: curl's output is hashed and extracted at the same
# time, so it's never stored (the step needs ~8GB of disk, not ~14GB). The
# hash is checked as soon as the stream ends, and a mismatch fails the
# build, so an unverified database never reaches the image.
RUN bash -euo pipefail -c ' \
      mkdir -p /opt/kraken2-db ; \
      mkfifo /tmp/kraken2-db.fifo ; \
      { sha256sum < /tmp/kraken2-db.fifo | cut -d" " -f1 > /tmp/kraken2-db.sha256 ; } & \
      hasher=$! ; \
      curl -fsSL "${KRAKEN2_DB_URL}" | tee /tmp/kraken2-db.fifo | tar -xz -C /opt/kraken2-db ; \
      wait "${hasher}" ; \
      actual="$(cat /tmp/kraken2-db.sha256)" ; \
      rm /tmp/kraken2-db.fifo /tmp/kraken2-db.sha256 ; \
      if [ "${actual}" != "${KRAKEN2_DB_SHA256}" ]; then \
        echo "Kraken2 database checksum mismatch: expected ${KRAKEN2_DB_SHA256}, got ${actual}" >&2 ; \
        exit 1 ; \
      fi ; \
      test -f /opt/kraken2-db/hash.k2d -a -f /opt/kraken2-db/opts.k2d -a -f /opt/kraken2-db/taxo.k2d ; \
      chmod -R a+rX /opt/kraken2-db'

COPY scrnaseq_qc_pipeline/checks/py_dcqc/kraken2_qc/conda-linux-64.lock /tmp/conda-linux-64.lock
RUN micromamba install -y -n base -f /tmp/conda-linux-64.lock && \
    micromamba clean -a -y && \
    rm /tmp/conda-linux-64.lock && \
    chmod -R a+rX /opt/conda

# The base image's ENTRYPOINT activates the environment; nf-dcqc supplies
# the full command, so clear it and put the environment on PATH instead.
ENTRYPOINT []
ENV PATH=/opt/conda/bin:${PATH}

# nf-dcqc runs the container as the user who runs the workflow, not as root,
# and that user has no home directory in the image.
ENV HOME=/tmp
ENV PYTHONDONTWRITEBYTECODE=1
WORKDIR /app

COPY scrnaseq_qc_pipeline/scrnaseq_common.py scrnaseq_qc_pipeline/kraken_qc.py /app/
COPY scrnaseq_qc_pipeline/checks/py_dcqc/kraken2_qc/run.py /app/run.py
RUN chmod 0755 /app/run.py
