# `checks/external/` vs `checks/internal/`

Two categories, and the split is **not** simply "wraps a tool" vs. "doesn't."
It follows py-dcqc's own real distinction (see
[Sage-Bionetworks-Workflows/py-dcqc#96](https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96),
"Internal and External Tests"): **where does the business logic run** —
inside py-dcqc's own Python process, or as a separately-invoked
command/container.

## `external/` — genuine py-dcqc `ExternalBaseTest` candidates

A check belongs here only if it satisfies **all** of:
- It wraps an external tool, or is pure computation with no network/tool
  dependency at all.
- It needs nothing but the file(s) already staged on disk — no network
  call, no credentials, at runtime.
- Its business logic runs as its own process (a container nf-dcqc would
  invoke), judged only by exit code + stdout + stderr.

These follow the real py-dcqc contract: no `ENTRYPOINT`/`CMD` in the
Dockerfile, invoked with bare filenames relative to the working directory,
a 3-way exit code (`0` pass / `1` fail — a genuine tool-native verdict /
`2` error — the check itself couldn't run), metrics as one JSON object on
stdout, a short reason on stderr. See `qc_common.py`'s `DCQC_PASS`/
`DCQC_FAIL`/`DCQC_ERROR`/`emit_dcqc_metrics`/`emit_dcqc_reason`, and each
check's own `run.py` docstring for its specific exit-code mapping.

Current members: `fq_lint`, `fastqc`, `kraken2_qc`, `matrix_content_qc`,
`fastq_aggregate`, `multiqc`.

The last two don't wrap a tool-on-one-file the way the first four do, but
they still satisfy every bullet above, with a documented deviation each
(see their own `run.py`/`Dockerfile` docstrings for the specifics):
- **`fastq_aggregate`** is pure computation over already-staged JSON
  summaries — no network, no credentials — but takes a variable number of
  input files, beyond py-dcqc's tested SingleTarget/PairedTarget cardinality.
- **`multiqc`** wraps the `multiqc` binary on files already staged, but
  scans an arbitrary-N-file directory rather than one file/pair, and writes
  its own report into the working directory rather than confining itself to
  the three reserved output files — the report is this check's entire
  purpose, unlike a validator whose only job is a verdict.

## `internal/` — everything else that still needs to run as its own script

**Important: this is not the same thing as py-dcqc's own `InternalBaseTest`.**
A real py-dcqc internal test has *no* separate script or container at all —
its Python code is contributed directly into py-dcqc's own package, and
`dcqc compute-test` calls it as a plain method in-process. There is no
third py-dcqc category for "a standalone script that isn't merged into
py-dcqc's codebase" — if logic needs to run as its own process, py-dcqc
calls that an *external* test, full stop, regardless of whether it wraps a
tool.

Everything in this folder is a check that:
- Could, in principle, become a real `InternalBaseTest` (pure Python, no
  tool) **only by contributing that code into py-dcqc's own repository** —
  which this project has explicitly not done, so as a standalone script of
  ours it doesn't fit py-dcqc's internal category either; **or**
- Needs credentials or does its own file staging (`matrix_download`,
  `synapse_fastq_download`, `fastq_download`) — which no py-dcqc test of
  either kind is allowed to do. Under real py-dcqc, Synapse/URL file
  staging isn't a test at all; it's handled automatically by py-dcqc's own
  `File.stage()` before any test runs, keyed off the file's registered
  `syn://`/HTTPS URL. There's no py-dcqc test to write for this.

So `internal/` really means: **our own scripts, run by our own tooling,
living outside py-dcqc's world entirely** — not a py-dcqc test category.
They keep this repo's original `--output <findings.json>` contract (see
`qc_common.write_findings_json`), not the dcqc exit-code contract, since
they were never going to be judged by `dcqc compute-test` in the first
place. The plan is to run these via this repo's **own** Nextflow pipeline
(separate from nf-dcqc) once that's built — not to make them py-dcqc tests.

Current members: `matrix_download`, `synapse_fastq_download`,
`fastq_download`, `fastq_resolution`, `geo_metadata_fetch`,
`geo_annotation_check`.

## A note on `fastq_download`

This one is a hybrid: it does its own HTTPS fetch (why it's `internal/`,
not `external/`) *and* a size/md5 comparison against ENA's registry once
the file is local. The comparison half alone — given an already-staged file
plus the registered size/md5 as plain scalars, no network — would qualify
as a legitimate 6th external test. It hasn't been split out; this is a
natural follow-up if tighter py-dcqc alignment is wanted later.

## An open question: `geo_annotation_check`

By the same reasoning that moved `fastq_aggregate` here-to-external
(pure computation over already-staged JSON, no network/credentials at
runtime), `geo_annotation_check` arguably qualifies too — it only compares
two already-fetched JSON dicts (`annotations.json`, `geo_metadata.json`);
the network call that produced them already happened in the separate
`geo_metadata_fetch` step. It's kept in `internal/` for now since this
hasn't been explicitly decided either way — revisit it together with
`fastq_download`'s split if/when this repo tightens py-dcqc alignment
further.
