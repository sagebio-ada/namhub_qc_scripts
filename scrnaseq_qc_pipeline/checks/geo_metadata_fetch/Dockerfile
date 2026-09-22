# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/geo_metadata_fetch/Dockerfile -t qc-geo-metadata-fetch .
FROM python:3.11-slim
WORKDIR /app

RUN pip install --no-cache-dir "git+https://github.com/sagebio-ada/geo_dataset_creation.git"

COPY qc_common.py geo_metadata_qc.py /app/
COPY checks/geo_metadata_fetch/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
