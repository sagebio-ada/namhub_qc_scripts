# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/matrix_content_qc/Dockerfile -t qc-matrix-content .
FROM python:3.11-slim
WORKDIR /app

RUN pip install --no-cache-dir numpy scipy

COPY qc_common.py matrix_qc.py /app/
COPY checks/matrix_content_qc/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
