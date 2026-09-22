# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/multiqc/Dockerfile -t qc-multiqc .
FROM quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0
WORKDIR /app

COPY qc_common.py multiqc_qc.py /app/
COPY checks/multiqc/run.py /app/run.py

ENTRYPOINT ["python", "run.py"]
