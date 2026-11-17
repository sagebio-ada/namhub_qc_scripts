# Build from the repo root (namhub_qc_scripts/), not scrnaseq_qc_pipeline/:
# fastq_qc.py imports fq_lint_qc, which imports qc_common/ from the repo root.
#   docker build -f scrnaseq_qc_pipeline/checks/py_dcqc/fastqc/Dockerfile -t qc-fastqc .
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     qc-fastqc /bin/bash -euo pipefail -c \
#     "( (python3 /app/run.py 'r1.fastq.gz' 'r2.fastq.gz') > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
#
# FastQC and MultiQC share this image because MultiQC can only summarize
# output written in its own container (py-dcqc runs each test in its own).
# Everything comes from one micromamba environment -- one package manager,
# no files copied between base images.
#
# MultiQC is pinned to 1.14: multiqc_qc.read_multiqc_data() reads the
# multiqc_data.json layout of that version. It still imports pkg_resources
# (setuptools < 70) and predates numpy 2.
#
# FastQC needs a real font to draw its plots. Without one its background
# thread crashes, it never finishes its output zip, and the run hangs until
# run_fastqc's timeout -- hence fontconfig + a DejaVu font.
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.10-bookworm-slim

USER root
# procps: Nextflow's task metrics run `ps` inside the container.
RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.10 fastqc=0.12.1 multiqc=1.14 "setuptools<70" "numpy<2" requests \
      fontconfig font-ttf-dejavu-sans-mono && \
    micromamba clean -a -y && \
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

COPY scrnaseq_qc_pipeline/scrnaseq_common.py scrnaseq_qc_pipeline/fastq_qc.py scrnaseq_qc_pipeline/fq_lint_qc.py \
     scrnaseq_qc_pipeline/kraken_qc.py scrnaseq_qc_pipeline/multiqc_qc.py /app/
COPY qc_common/ /app/qc_common/
COPY scrnaseq_qc_pipeline/checks/py_dcqc/fastqc/run.py /app/run.py
RUN chmod 0755 /app/run.py
