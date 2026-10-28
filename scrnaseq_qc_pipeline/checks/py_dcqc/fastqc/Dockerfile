# Build from the repo root (namhub_qc_scripts/), not scrnaseq_qc_pipeline/:
# this check's fq_lint_qc.py imports qc_common.fq_lint, and qc_common/ lives
# at the repo root, outside scrnaseq_qc_pipeline/ -- Docker can't COPY a path
# from outside the build context, so the context has to be the repo root.
#   docker build -f scrnaseq_qc_pipeline/checks/py_dcqc/fastqc/Dockerfile -t qc-fastqc .
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     qc-fastqc /bin/bash -euo pipefail -c \
#     "( (python3 /app/run.py 'sample.fastq.gz') > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
# The biocontainers single-package image only ships FastQC's own runtime
# deps (openjdk/perl, no python). It looks Debian-based (buster) but is a
# conda-built "mulled" image with no package manager at all -- apt-get
# itself is missing, not just the packages -- so python3 can't be apt-get
# installed into it directly (verified by inspecting the image: dpkg and
# wget exist, apt-get doesn't). This stage installs python3, python3-requests
# (fastq_qc.py, copied whole -- see the dockerize-checks skill -- imports
# requests at module load time even though this check doesn't call the
# download functions) and procps (Nextflow's task-metrics `ps`) on plain
# debian:buster-slim -- the same Debian release the final image's files are
# built on -- so what gets copied into the final stage is ABI-compatible
# with its libc. buster is EOL, so its sources.list must point at
# archive.debian.org, not deb.debian.org, or apt-get update 404s.
# --platform=linux/amd64 is pinned on both stages: the final image is only
# published for amd64 (see the InvalidBaseImagePlatform warning this threw
# on an arm64 build host before this was pinned), and building pydeps for a
# different architecture than the final stage copies incompatible binaries.
FROM --platform=linux/amd64 debian:buster-slim AS pydeps
RUN echo "deb http://archive.debian.org/debian buster main" > /etc/apt/sources.list && \
    echo "deb http://archive.debian.org/debian-security buster/updates main" >> /etc/apt/sources.list && \
    apt-get -o Acquire::Check-Valid-Until=false update && \
    apt-get -o Acquire::Check-Valid-Until=false install -y --no-install-recommends \
      python3 python3-requests procps && \
    rm -rf /var/lib/apt/lists/*

FROM --platform=linux/amd64 quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
WORKDIR /app

# Only /usr/bin, /usr/lib, /lib and /etc/alternatives are copied -- not all
# of /usr -- because apt installs system packages there, not under
# /usr/local, and /usr/local is where the final image's own conda-installed
# fastqc/java/perl live; copying /usr wholesale conflicts with it (e.g.
# /usr/local/man is a real directory in one image and a symlink in the
# other). /lib is separate from /usr/lib on buster (no usr-merge) and holds
# python3's own shared-library dependencies (libexpat, libz, ...), so it has
# to be copied too, or python3 fails at runtime with "libexpat.so.1: cannot
# open shared object file". /bin is where procps installs the real `ps`
# binary (not /usr/bin) -- without it, `ps` still resolves to the image's
# existing BusyBox multi-call binary, whose `ps` doesn't support the output
# Nextflow's task-metrics step asks for.
COPY --from=pydeps /usr/bin /usr/bin
COPY --from=pydeps /usr/lib /usr/lib
COPY --from=pydeps /lib /lib
COPY --from=pydeps /bin /bin
COPY --from=pydeps /etc/alternatives /etc/alternatives

# nf-dcqc runs the container as the user who runs the workflow, not as root,
# and that user has no home directory in the image.
ENV HOME=/tmp
ENV PYTHONDONTWRITEBYTECODE=1

# fq_lint_qc.py imports qc_common.fq_lint, so qc_common/ is copied alongside
# it (see fq_lint_qc.py's sys.path logic, which looks for a qc_common/
# directory next to itself before falling back to the repo-checkout layout).
COPY scrnaseq_qc_pipeline/scrnaseq_common.py scrnaseq_qc_pipeline/fastq_qc.py scrnaseq_qc_pipeline/fq_lint_qc.py scrnaseq_qc_pipeline/kraken_qc.py /app/
COPY qc_common/ /app/qc_common/
COPY scrnaseq_qc_pipeline/checks/py_dcqc/fastqc/run.py /app/run.py
RUN chmod 0755 /app/run.py
