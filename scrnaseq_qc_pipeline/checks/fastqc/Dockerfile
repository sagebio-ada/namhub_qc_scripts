# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/fastqc/Dockerfile -t qc-fastqc .
FROM quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
WORKDIR /app

# The biocontainers single-package image only ships FastQC's own runtime
# deps (openjdk/perl, no python) -- it's Debian-based, so apt-get is the
# simplest way to add python3. python3-requests is needed because fastq_qc.py
# (copied whole -- see the dockerize-checks skill) imports requests at module
# load time even though this check doesn't call the download functions.
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-requests && \
    rm -rf /var/lib/apt/lists/*

COPY qc_common.py fastq_qc.py fq_lint_qc.py kraken_qc.py /app/
COPY checks/fastqc/run.py /app/run.py

ENTRYPOINT ["python3", "run.py"]
