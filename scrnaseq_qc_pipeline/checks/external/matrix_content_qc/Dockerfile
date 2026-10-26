# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/external/matrix_content_qc/Dockerfile -t qc-matrix-content .
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     qc-matrix-content /bin/bash -euo pipefail -c \
#     "( (python3 /app/run.py 'barcodes.tsv.gz' 'features.tsv.gz' 'matrix.mtx.gz') > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
FROM python:3.14-slim
WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir numpy scipy

# nf-dcqc runs the container as the user who runs the workflow, not as root,
# and that user has no home directory in the image.
ENV HOME=/tmp
ENV PYTHONDONTWRITEBYTECODE=1

COPY qc_common.py matrix_qc.py /app/
COPY checks/external/matrix_content_qc/run.py /app/run.py
RUN chmod 0755 /app/run.py
