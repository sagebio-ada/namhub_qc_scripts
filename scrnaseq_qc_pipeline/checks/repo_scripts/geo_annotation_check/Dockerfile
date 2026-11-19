# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/repo_scripts/geo_annotation_check/Dockerfile -t qc-geo-annotation-check .
FROM python:3.14-slim
WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends git && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir "git+https://github.com/sagebio-ada/geo_dataset_creation.git"

COPY scrnaseq_common.py geo_metadata_qc.py /app/
COPY checks/repo_scripts/geo_annotation_check/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
