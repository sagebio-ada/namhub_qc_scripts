# `checks/py_dcqc/` vs `checks/repo_scripts/`

Two categories, and the split is **not** simply "wraps a tool" vs. "doesn't."
It follows py-dcqc's own real distinction between internal and external tests
(see [Sage-Bionetworks-Workflows/py-dcqc#96](https://github.com/Sage-Bionetworks-Workflows/py-dcqc/pull/96),
"Internal and External Tests"): **where does the business logic run** —
inside py-dcqc's own Python process, or as a separately-invoked
command/container. (Folder names here avoid py-dcqc's own "internal"/
"external" vocabulary on purpose — see why under `repo_scripts/` below.)

## `py_dcqc/` — genuine py-dcqc `ExternalBaseTest` candidates

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
stdout, a short reason on stderr. See `scrnaseq_common.py`'s `DCQC_PASS`/
`DCQC_FAIL`/`DCQC_ERROR`/`emit_dcqc_metrics`/`emit_dcqc_reason`, and each
check's own `run.py` docstring for its specific exit-code mapping. (This
module was renamed from `qc_common.py` — that name now belongs to the
repo-root `qc_common/` package shared across the other pipelines; the two
are unrelated and were colliding in name only.)

Current members: `fq_lint`, `raw_read_qc`, `matrix_content_qc`.

## Decided: `fastqc`/`kraken2_qc`/`fastq_aggregate`/`multiqc` merged into `raw_read_qc`

These four used to be separate checks. andrewelamb's PR #1 review
(2026-11-03) found that two of them could never actually work as
independent py-dcqc tests: py-dcqc runs every `ExternalBaseTest` in its own
isolated container with only its own target file(s) staged, so
- **`fastq_aggregate`** (needs every mate file's FastQC summary) and
- **`multiqc`** (scans the working directory for FastQC/Kraken2 output that
  other tools already wrote)

could never see another test's output — there was no file-passing mechanism
between separate py-dcqc containers to begin with, so the "documented
deviation" each carried before this wasn't just a cardinality quirk, it was
a hard blocker. The only way to keep this functionality inside py-dcqc at
all was to put the tools that *produce* the needed files (`fastqc`,
`kraken2_qc`) in the same container as the ones that *consume* them. The
same review also flagged `kraken2_qc`'s runtime database download as a
separate py-dcqc violation ("the tool must not download data when it
runs") — fixed in the same merge by baking Kraken2's full Standard-8 index
into the image at build time instead.

`raw_read_qc` now does what all four used to, in one container, on 1–2
mate fastq files (a py-dcqc SingleTarget/PairedTarget): FastQC + Kraken2
per mate, the ENA-registry aggregate comparison across mates, and a
MultiQC report over just this run's own FastQC/Kraken2 output (not a
cross-sample report the way `run_scrnaseq_qc.py`'s own local driver
produces one — a DCQC test is inherently per-target, so cross-sample
aggregation isn't reachable from inside one). See `raw_read_qc/run.py`'s
own docstring for the full exit-code contract.

## `repo_scripts/` — our own scripts, our own contract

**Important: this is not the same thing as py-dcqc's own `InternalBaseTest`.**
A real py-dcqc internal test has *no* separate script or container at all —
its Python code is contributed directly into py-dcqc's own package, and
`dcqc compute-test` calls it as a plain method in-process. There is no
third py-dcqc category for "a standalone script that isn't merged into
py-dcqc's codebase" — if logic needs to run as its own process, py-dcqc
calls that an *external* test, full stop, regardless of whether it wraps a
tool. That's exactly why this folder isn't called `internal/`: nothing in
it is a py-dcqc internal test, and naming it that way invited exactly that
confusion.

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

So `repo_scripts/` really means: **our own scripts, run by our own tooling,
living outside py-dcqc's world entirely** — not a py-dcqc test category.
They keep this repo's original `--output <findings.json>` contract (see
`scrnaseq_common.write_findings_json`), not the dcqc exit-code contract, since
they were never going to be judged by `dcqc compute-test` in the first
place. The plan is to run these via this repo's **own** Nextflow pipeline
(separate from nf-dcqc) once that's built — not to make them py-dcqc tests.

Current members: `matrix_download`, `synapse_fastq_download`,
`fastq_download`, `fastq_resolution`, `geo_metadata_fetch`,
`geo_annotation_check`.

## Decided: `fastq_download`'s comparison half and `geo_annotation_check` stay in `repo_scripts/`

Both satisfy the `py_dcqc/` membership bullets on paper — `fastq_download`'s
size/md5-against-ENA comparison, and `geo_annotation_check`'s diff of two
already-fetched JSON dicts, are each pure computation over already-staged
data with no network or credentials at runtime. But satisfying the
criteria isn't the same as being worth moving: neither check is actually
going to run through py-dcqc/nf-dcqc in practice, so reclassifying them
into `py_dcqc/` would just relocate the code without changing how either
one is invoked. Decision: leave both where they are. `fastq_download`
stays a single hybrid script (fetch + comparison, unsplit); `geo_annotation_check`
stays in `repo_scripts/`.
