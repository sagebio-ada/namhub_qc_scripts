---
name: dockerize-checks
description: >
  Splits a script that bundles multiple QC/validation checks into standalone,
  independently-dockerizable single-purpose scripts, following this repo's
  checks/<name>/{run.py, Dockerfile} convention. Use when asked to make a
  script's checks dockerizable, containerize a pipeline's checks, or split
  out checks in preparation for nf-dcqc/py-dcqc-style integration.
---

# Dockerize checks

A generalized procedure for turning a script that bundles multiple QC/validation
checks into a set of small, standalone, independently-dockerizable scripts — one
per check (or per tightly-coupled cluster of checks). It does **not** assume
anything about the target script's domain; the worked example is
`scrnaseq_qc_pipeline/`, linked at the end, but every step here applies to any
script in this repo that emits more than one distinct pass/fail/info-style
result from one big function.

This produces standalone scripts + Dockerfiles only. It does not integrate with
any external test framework (e.g. py-dcqc) or workflow engine (e.g. Nextflow) —
that is a distinct, later step this convention is meant to make easier, not
something this skill does itself.

## 1. Enumerate every check

Grep the target script for its finding/result-emission calls (in this repo,
`make_finding(...)` calls; in another script, whatever function reports one
verdict/fact). For each one, note:

- **What tool/library it depends on** — an external CLI binary (fastqc,
  samtools, a linter), a Python library, or nothing beyond the standard library.
- **What its inputs are** — a file already on disk, a network fetch (an API,
  a registry lookup), or both.
- **What it shares with neighboring checks** — the same downloaded file, the
  same subprocess's parsed output, an in-memory object (a parsed matrix, a
  dict of metadata) built by an earlier check in the same function.

This step is pure reading — don't refactor yet. The goal is a table like:
`check name | tool dep | python deps | input source | entangled with`.

## 2. Decide fusion vs. split boundaries, and which units are real external tests

Fusion/split is a rule, not a case-by-case guess:

- **Fuse** checks that share *expensive-to-reproduce* state: the same
  downloaded file, the same subprocess invocation, an in-memory object with no
  clean serialization (e.g. a sparse matrix). Splitting these would force an
  artifact hand-off that only exists to route around the split, for no benefit.
- **Split** along natural *external-tool* boundaries — each distinct
  binary/CLI tool gets its own function (and, later, its own container) — and
  along genuinely independent computations that don't share expensive state.
- When candidate units share state that's small and cheap to recompute (e.g.
  re-parsing a small file, redoing a cheap numeric pass), prefer splitting for
  clarity over fusing just to avoid the redundant recompute.

Then classify each resulting unit — this determines which contract it gets in
step 4, and it is **not** simply "wraps a tool" vs. "doesn't." Following
py-dcqc's own real distinction (Sage-Bionetworks-Workflows/py-dcqc#96,
"Internal and External Tests"): the question is *where the business logic
runs*, not what it depends on.

- **A genuine external-test candidate** needs nothing but the file(s) already
  on disk — no network call, no credentials, at runtime — regardless of
  whether it wraps a tool or is pure computation. These get the strict
  py-dcqc contract in step 4.
- **Everything else that still needs its own script/process** — a network
  API call, a Synapse/URL download, pure Python logic with no tool
  dependency — is *not* a py-dcqc internal test just because it has no tool
  dependency. Real py-dcqc internal tests have no separate script at all:
  the code is contributed directly into py-dcqc's own package and called
  in-process. A standalone script of ours that was never contributed
  upstream doesn't fit py-dcqc's internal category either — it's simply
  outside py-dcqc's two real categories, run by this repo's own tooling
  instead. See `scrnaseq_qc_pipeline/checks/README.md` for the fuller
  version of this distinction and why `matrix_download`-style
  Synapse-staging checks don't become py-dcqc tests of either kind (py-dcqc
  does file staging itself, before any test runs — there's no test to
  write for it).

Write down the resulting list of units, each tagged external or not, before
touching code — it's the scaffold for steps 3-4.

## 3. Refactor the source module in place

Break the tangled function into the small composable functions decided in
step 2, **in the module where the logic already lives** — never duplicate
logic into a new location. The refactored function should keep the exact same
external signature and behavior as before if anything still calls it as one
unit (e.g. an existing orchestrator script); it just becomes a thin wrapper
chaining the new pieces in the same order.

Preserve whatever verdict philosophy the target script already uses. In this
repo: no self-asserted PASS/FAIL/WARN — only a genuine tool-native verdict
(the tool's own exit code, or a threshold the tool itself publishes) may carry
one; everything else is `INFO` (a fact, no judgment) or `FAIL` only for a
real operational failure (an exception, not a data-quality opinion).

Verify the refactor didn't change behavior before moving on — run the
before/after versions against the same fixture input and diff their output.

## 4. Wrap each resulting unit in `checks/{py_dcqc,repo_scripts}/<name>/{run.py, Dockerfile}`

Directory layout, one pair per check, sorted into the category decided in
step 2. Name the two folders for what they actually are, not py-dcqc's own
"internal"/"external" vocabulary — a folder called `internal/` invites
exactly the confusion the note under the second folder below exists to
head off:

```
<pipeline_dir>/
  <existing modules, refactored per step 3 — unchanged otherwise>
  checks/
    README.md        # the py_dcqc/repo_scripts distinction, written out in full
    py_dcqc/
      <check_name>/{run.py, Dockerfile}
    repo_scripts/
      <check_name>/{run.py, Dockerfile}
```

### `py_dcqc/` — the real py-dcqc `ExternalBaseTest` contract

Follow Sage-Bionetworks-Workflows/py-dcqc#96 exactly — it is the canonical
spec, written by py-dcqc's own maintainers for exactly this purpose:

- **`run.py`** takes the file(s) as **bare filenames, positional, relative to
  the working directory** (nf-dcqc stages each file there as a symlink) —
  never a flag, never a directory prefix. Writes **one JSON object of
  metrics to stdout** (`scrnaseq_common.emit_dcqc_metrics`) and, on a non-pass
  result, a short reason to **stderr** (`scrnaseq_common.emit_dcqc_reason`).
  Nothing else goes to stdout — any other text makes the JSON invalid.
  (`scrnaseq_common.py` — formerly `qc_common.py`, renamed to stop colliding
  with the unrelated repo-root `qc_common/` package used by the other
  pipelines.)
- **Exit code is 3-way**, not 2-way (`scrnaseq_common.DCQC_PASS`/`DCQC_FAIL`/
  `DCQC_ERROR` = 0/1/2): `0` the file passed; `1` a genuine tool-native
  bad-file verdict; `2` the check itself couldn't run (crash, missing
  binary) — not a verdict about the data. If the check has no validated
  tool-native verdict at all (most of this pipeline's checks don't — see
  the INFO-only philosophy in step 3), it only ever resolves to `0` or `2`;
  say so explicitly in the docstring rather than leaving `1` a mystery.
- **`Dockerfile` has no `ENTRYPOINT`/`CMD`** — nf-dcqc supplies the full
  command every time. Needs `bash` and `procps` (`apt-get install procps`
  on Debian-based images — Nextflow uses `ps` for task metrics). Set
  `ENV HOME=/tmp` and `ENV PYTHONDONTWRITEBYTECODE=1` (the container runs as
  an arbitrary non-root UID with no home directory). `chmod 0755` the copied
  `run.py`. Comment the Dockerfile with the exact `docker run` invocation
  that tests it the way nf-dcqc would (see any `checks/py_dcqc/*/Dockerfile`
  in this repo for the template) — this is also how you manually verify the
  image before anything automated does.
- If a check needs an artifact from an earlier step (step 5) for something
  beyond the strict py-dcqc contract — this repo's own aggregate checks
  need each mate file's FastQC summary, which py-dcqc has no mechanism for
  passing between independent tests — add an explicit opt-in flag for it
  (e.g. `--summary-out`) and say in the docstring that it's an extension a
  literal py-dcqc submission would drop.
- **Base image**: prefer an existing community image
  (`quay.io/biocontainers/<tool>`, `ghcr.io/...`) for any check wrapping an
  external bioinformatics/CLI tool. Use a minimal custom image (e.g.
  `python:3.14-slim` + pip installs) only for pure-language logic, or a tool
  with no existing community image.

### `repo_scripts/` — our own scripts, our own contract

Not py-dcqc's own "internal test" category — see `checks/README.md`'s fuller
explanation of why that name was avoided (a real py-dcqc internal test has
no separate script at all; its code is contributed into py-dcqc's own
package). This folder is simply everything that still needs to run as its
own process but doesn't meet the `py_dcqc/` bar in step 2 (needs network or
credentials at runtime, or has no tool dependency and was never contributed
upstream).

Keep the simpler, original contract: `run.py` takes identifying flags (an
entity/sample id) and `--output PATH`, writes a **JSON array of Findings**
to it (`scrnaseq_common.write_findings_json`), exits `0` whenever it produced
valid output (including a FAIL finding — that's the check doing its job)
and non-zero only if the wrapper's own control flow crashed. `Dockerfile`
keeps its `ENTRYPOINT` — these were never going to be judged by `dcqc
compute-test`, so there's no reason to drop the convenience. Do not try to
force these into the `py_dcqc/` contract just for consistency; see
`scrnaseq_qc_pipeline/checks/README.md` for exactly which checks belong here
and why.

## 5. Define the artifact-passing contract

Containers can't share Python objects. For any state that must cross a
container boundary (identified as "too large/expensive to recompute" in step
2 — small, cheap-to-recompute state should just be recomputed, not passed):

- One small JSON file per artifact, written to a shared, host-mounted work
  directory.
- An explicit, minimal schema — document it once, next to the check that
  produces it. Don't invent a generic serialization framework; a plain `dict`
  written with `json.dump` is enough.
- Prefer passing through an existing data shape unchanged (e.g. if an
  upstream library call already returns exactly the structure a downstream
  check needs, write that structure to JSON as-is rather than reshaping it).

## 6. Add or update a local end-to-end driver

One script (plain language script calling into each container via its normal
process-invocation mechanism — e.g. `subprocess.run(["docker", "run", ...])`
in Python) that chains the new per-check containers together against one
shared mounted work directory, for the team's own local validation — not a
production orchestrator. Two invocation styles, matching step 4's two
contracts:

- **`repo_scripts/` checks**: unchanged — one shared mount, pass `--output`,
  read the Findings JSON back from the host.
- **`py_dcqc/` checks**: invoke them exactly the way nf-dcqc's
  `RUN_PROCESS` would (`modules/local/run_process.nf` in
  Sage-Bionetworks-Workflows/nf-dcqc) — a fresh task directory per
  invocation holding only the input file(s) symlinked in by their bare
  filename, the container run as the current host user (`-u
  "$(id -u):$(id -g)"`) with that directory as both the mount point and the
  working directory, and the result read back from the three files the
  wrapping bash produces (`std_out.txt`/`std_err.txt`/`exit_code.txt`) —
  never from a `--output` flag, since these containers don't have one. Write
  one small adapter function per external check that maps its 3-way exit
  code + stdout JSON back into this repo's own Finding shape for the
  driver's combined report; py-dcqc's one-test-one-verdict model doesn't
  preserve the original per-check Finding granularity 1:1 (e.g. a
  multi-module tool like FastQC collapses to one exit code), so don't expect
  a perfect byte-for-byte match against the pre-split script's output here —
  only that the real verdicts and numbers agree.

Other driver rules, for both styles:

- Mirror the original (pre-split) script's own control-flow branching, so its
  output can be diffed against that script's output as the validation step.
- Favor a real scripting language over a static tool like `docker-compose`
  whenever the pipeline has data-dependent branching or a fan-out that isn't
  known until an earlier step runs (e.g. "loop once per file discovered by
  step 1") — compose's static service graph fights that kind of pipeline.
- After each container invocation, read that step's output back on the
  host to decide the next step's arguments — never inside a container.
- If the original script produces a combined report (a CSV, a summary), the
  driver should produce the same shape from the concatenated per-step results,
  so the two can be diffed directly.

## Worked example

`scrnaseq_qc_pipeline/` in this repo applies every step above:
`fastq_qc.py` and `matrix_qc.py` hold the refactored composable functions
(step 3), `scrnaseq_qc_pipeline/checks/{py_dcqc,repo_scripts}/*/` holds the
per-check `run.py`/`Dockerfile` pairs (step 4, with
`scrnaseq_qc_pipeline/checks/README.md` writing out the py_dcqc/repo_scripts
distinction in full), and `scrnaseq_qc_pipeline/driver/run_dockerized_pipeline.py`
is the local end-to-end driver (step 6, with both invocation styles). Read
those as a concrete reference when applying this procedure to a new script.
