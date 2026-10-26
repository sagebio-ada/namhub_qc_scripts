# Build from the scrnaseq_qc_pipeline/ directory:
#   docker build -f checks/external/fq_lint/Dockerfile -t qc-fq-lint .
#
# No biocontainers/quay.io image exists for stjude-rust-labs/fq, so this
# fetches the prebuilt Linux release binary directly. Pin/verify the version
# against https://github.com/stjude-rust-labs/fq/releases on update.
#
# py-dcqc ExternalBaseTest image (Sage-Bionetworks-Workflows/py-dcqc#96):
# no ENTRYPOINT/CMD -- nf-dcqc supplies the full command. Test it the same
# way nf-dcqc would run it:
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" \
#     qc-fq-lint /bin/bash -euo pipefail -c \
#     "( (python3 /app/run.py 'sample.fastq.gz') > std_out.txt 2> std_err.txt; echo \$? > exit_code.txt ) || true"
FROM debian:bookworm-slim
WORKDIR /app

ARG FQ_VERSION=0.12.0
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates python3 procps && \
    rm -rf /var/lib/apt/lists/* && \
    curl -fsSL -o /tmp/fq.tar.gz \
      "https://github.com/stjude-rust-labs/fq/releases/download/v${FQ_VERSION}/fq-${FQ_VERSION}-x86_64-unknown-linux-gnu.tar.gz" && \
    tar -xzf /tmp/fq.tar.gz -C /tmp && \
    mv /tmp/fq-${FQ_VERSION}-x86_64-unknown-linux-gnu/fq /usr/local/bin/fq && \
    chmod +x /usr/local/bin/fq && \
    rm -rf /tmp/fq.tar.gz /tmp/fq-${FQ_VERSION}-x86_64-unknown-linux-gnu

# nf-dcqc runs the container as the user who runs the workflow, not as root,
# and that user has no home directory in the image.
ENV HOME=/tmp
ENV PYTHONDONTWRITEBYTECODE=1

COPY qc_common.py fq_lint_qc.py /app/
COPY checks/external/fq_lint/run.py /app/run.py
RUN chmod 0755 /app/run.py
