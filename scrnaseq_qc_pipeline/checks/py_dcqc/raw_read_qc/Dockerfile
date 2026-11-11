# Build from the repo root (namhub_qc_scripts/), not scrnaseq_qc_pipeline/:
# this check's fq_lint_qc.py imports qc_common.fq_lint, and qc_common/ lives
# at the repo root, outside scrnaseq_qc_pipeline/ -- Docker can't COPY a path
# from outside the build context, so the context has to be the repo root.
#   docker build -f scrnaseq_qc_pipeline/checks/py_dcqc/raw_read_qc/Dockerfile -t qc-raw-read-qc .
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     qc-raw-read-qc /bin/bash -euo pipefail -c \
#     "( (python3 /app/run.py 'r1.fastq.gz' 'r2.fastq.gz') > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
#
# Merges fastqc + kraken2_qc + fastq_aggregate + multiqc into one test, per
# andrewelamb's PR #1 review (2026-11-03): py-dcqc runs each test in its own
# isolated container with only its own target file(s) staged, so
# fastq_aggregate's cross-mate comparison and multiqc's "scan the working
# directory for other tools' output" could never see fastqc/kraken2's
# output from a separate container. Putting all four in one image/test is
# the only way for that to work. See run.py's own docstring for the full
# contract and the per-check history this replaces.
#
# Three base images are combined here, none of which share a package
# manager or conda root, so each tool's own /usr/local is copied to its own
# non-colliding /opt path rather than merged into the final image's
# /usr/local (which belongs to the final FROM image's own conda install --
# see checks/py_dcqc/fastqc/Dockerfile for why merging conda roots directly
# is unsafe) :
#   - pydeps: python3 (for run.py itself and this file's own modules),
#     procps (Nextflow's task-metrics `ps`), and MultiQC (a pip package,
#     with its own matplotlib/numpy/jinja2 dependency tree) -- all on plain
#     debian:buster-slim, the same Debian release the final image's own
#     files are built on, for libc/ABI compatibility. buster is EOL, so its
#     sources.list points at archive.debian.org, not deb.debian.org.
#   - fastqc_src: just quay.io/biocontainers/fastqc, used only as a COPY
#     source for its own conda-installed /usr/local (fastqc + its
#     openjdk/perl runtime).
#   - final: quay.io/biocontainers/kraken2, whose own conda-installed
#     /usr/local already holds kraken2 + its runtime deps.
# --platform=linux/amd64 is pinned on every stage: the biocontainers images
# are amd64-only, and building pydeps/fastqc_src for a different
# architecture than the final stage would copy incompatible binaries.
FROM --platform=linux/amd64 debian:buster-slim AS pydeps
RUN echo "deb http://archive.debian.org/debian buster main" > /etc/apt/sources.list && \
    echo "deb http://archive.debian.org/debian-security buster/updates main" >> /etc/apt/sources.list && \
    apt-get -o Acquire::Check-Valid-Until=false update && \
    apt-get -o Acquire::Check-Valid-Until=false install -y --no-install-recommends \
      python3 python3-pip python3-setuptools python3-wheel python3-requests procps \
      fontconfig fonts-dejavu-core && \
    rm -rf /var/lib/apt/lists/* && \
    pip3 install --no-cache-dir --upgrade "pip<24" setuptools wheel && \
    pip3 install --no-cache-dir "multiqc==1.14" && \
    fc-cache -f && \
    chmod -R a+rX /etc/fonts /usr/share/fonts /var/cache/fontconfig

FROM --platform=linux/amd64 quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 AS fastqc_src

FROM --platform=linux/amd64 quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0
# The biocontainers base image bakes in its own ENTRYPOINT
# (/usr/local/env-execute), a conda-activation wrapper that RESETS PATH
# (and sets CONDA_*) in whatever command it execs -- silently wiping out
# the PATH extension below for /opt/pydeps-usrlocal and /opt/fastqc-conda.
# This has to be cleared, matching this file's own stated design ("no
# ENTRYPOINT/CMD -- nf-dcqc supplies the full command").
ENTRYPOINT []
WORKDIR /app

# python3 + procps -- see checks/py_dcqc/fastqc/Dockerfile for the full
# rationale on each copied path (apt-managed paths only, not /usr/local).
COPY --from=pydeps /usr/bin /usr/bin
COPY --from=pydeps /usr/lib /usr/lib
COPY --from=pydeps /lib /lib
COPY --from=pydeps /bin /bin
COPY --from=pydeps /etc/alternatives /etc/alternatives
# Without this, python3's urllib/ssl has no CA bundle to validate against --
# download_kraken2_db()'s HTTPS fetch fails with CERTIFICATE_VERIFY_FAILED.
# python3's own compiled-in default cert path is /usr/lib/ssl/certs (not
# /etc/ssl/certs, where the bundle actually landed), so SSL_CERT_FILE is set
# explicitly below rather than relying on that default resolving correctly.
COPY --from=pydeps /etc/ssl /etc/ssl
ENV SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt

# Without these, FastQC's --extract-free HTML/image report generation
# crashes: its per-base-quality plot needs AWT font metrics, which on
# modern OpenJDK means a real libfontconfig.so + /etc/fonts/fonts.conf +
# at least one actual font file -- none of which exist in either base
# image. That crash happens in FastQC's own background analysis thread
# (uncaught, so FastQC itself never notices or reports it), which means
# the main thread waiting on it never gets a completion signal and the
# whole process hangs indefinitely with its output zip stuck half-written
# -- reproduced identically on both Rosetta-emulated and genuine x86_64
# hardware, so this is the real root cause of that hang, not emulation.
COPY --from=pydeps /etc/fonts /etc/fonts
COPY --from=pydeps /usr/share/fonts /usr/share/fonts
COPY --from=pydeps /var/cache/fontconfig /var/cache/fontconfig

# MultiQC + its pip-installed dependency tree, and FastQC's own conda
# install, each kept in their own /opt path rather than merged into this
# image's own /usr/local (kraken2's conda root) -- see the file-level
# comment above for why.
COPY --from=pydeps /usr/local /opt/pydeps-usrlocal
COPY --from=fastqc_src /usr/local /opt/fastqc-conda
ENV PATH="/opt/pydeps-usrlocal/bin:/opt/fastqc-conda/bin:${PATH}"
ENV PYTHONPATH="/opt/pydeps-usrlocal/lib/python3.7/dist-packages"

# nf-dcqc runs the container as the user who runs the workflow, not as root,
# and that user has no home directory in the image.
ENV HOME=/tmp
ENV PYTHONDONTWRITEBYTECODE=1

# Only what download_kraken2_db() itself needs (via kraken_qc.py, which
# imports scrnaseq_common at module load) -- copied first and separately
# from the rest of this check's files so editing run.py later doesn't bust
# the cache on this layer and force a ~5.5GB re-download.
COPY scrnaseq_qc_pipeline/scrnaseq_common.py scrnaseq_qc_pipeline/kraken_qc.py /app/

# KNOWN DEVIATION resolved: earlier checks/py_dcqc/kraken2_qc/Dockerfile
# downloaded this database at container run time, which py-dcqc's own "What
# Makes a Good Test" rule disallows ("the tool must not download data when
# it runs"). Baked in at build time instead -- Kraken2's own full
# "Standard-8" index (~5.5GB compressed), so there's no coverage tradeoff
# against a reduced/custom reference.
RUN python3 -c "from pathlib import Path; from kraken_qc import download_kraken2_db, kraken2_db_staged; \
    download_kraken2_db(Path('/opt/kraken2-db')); \
    assert kraken2_db_staged(Path('/opt/kraken2-db')), 'Kraken2 database failed to stage'"

COPY scrnaseq_qc_pipeline/fastq_qc.py scrnaseq_qc_pipeline/fq_lint_qc.py scrnaseq_qc_pipeline/multiqc_qc.py /app/
COPY qc_common/ /app/qc_common/
COPY scrnaseq_qc_pipeline/checks/py_dcqc/raw_read_qc/run.py /app/run.py
RUN chmod 0755 /app/run.py
