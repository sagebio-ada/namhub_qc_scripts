# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/fq_lint/Dockerfile -t qc-fq-lint .
#
# No biocontainers/quay.io image exists for stjude-rust-labs/fq, so this
# fetches the prebuilt Linux release binary directly. Pin/verify the version
# against https://github.com/stjude-rust-labs/fq/releases on update.
FROM debian:bookworm-slim
WORKDIR /app

ARG FQ_VERSION=0.12.0
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates python3 && \
    rm -rf /var/lib/apt/lists/* && \
    curl -fsSL -o /tmp/fq.tar.gz \
      "https://github.com/stjude-rust-labs/fq/releases/download/v${FQ_VERSION}/fq-${FQ_VERSION}-x86_64-unknown-linux-gnu.tar.gz" && \
    tar -xzf /tmp/fq.tar.gz -C /tmp && \
    mv /tmp/fq-${FQ_VERSION}-x86_64-unknown-linux-gnu/fq /usr/local/bin/fq && \
    chmod +x /usr/local/bin/fq && \
    rm -rf /tmp/fq.tar.gz /tmp/fq-${FQ_VERSION}-x86_64-unknown-linux-gnu

COPY qc_common.py fq_lint_qc.py /app/
COPY checks/fq_lint/run.py /app/run.py

ENTRYPOINT ["python3", "run.py"]
