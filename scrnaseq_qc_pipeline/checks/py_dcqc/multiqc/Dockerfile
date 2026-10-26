# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/py_dcqc/multiqc/Dockerfile -t qc-multiqc .
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it, from a directory already containing the
# FastQC/Kraken2 output to aggregate:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     qc-multiqc /bin/bash -euo pipefail -c \
#     "( (python /app/run.py) > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
#
# KNOWN DEVIATIONS -- see run.py's docstring: this needs an arbitrary number
# of already-staged files (beyond py-dcqc's SingleTarget/PairedTarget model)
# and writes its own report into the working directory (not just the three
# reserved output files), since the report is this check's entire purpose.
FROM quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0
WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

# nf-dcqc runs the container as the user who runs the workflow, not as root,
# and that user has no home directory in the image.
ENV HOME=/tmp
ENV PYTHONDONTWRITEBYTECODE=1

COPY qc_common.py multiqc_qc.py /app/
COPY checks/py_dcqc/multiqc/run.py /app/run.py
RUN chmod 0755 /app/run.py
