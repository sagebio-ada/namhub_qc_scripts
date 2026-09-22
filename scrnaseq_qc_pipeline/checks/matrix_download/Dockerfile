# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/matrix_download/Dockerfile -t qc-matrix-download .
# Run with SYNAPSE_AUTH_TOKEN set in the environment.
FROM python:3.11-slim
WORKDIR /app

# numpy is a hard module-level import in matrix_qc.py even though this check
# doesn't use it directly (scipy is only imported lazily inside
# qc_matrix_content, so it's not needed here).
RUN pip install --no-cache-dir synapseclient numpy

COPY qc_common.py matrix_qc.py /app/
COPY checks/matrix_download/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
