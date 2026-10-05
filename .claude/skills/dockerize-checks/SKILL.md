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

## 2. Decide fusion vs. split boundaries

This is the one real judgment call, and it's a rule, not a case-by-case guess:

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

Write down the resulting list of units before touching code — it's the
scaffold for steps 3-4.

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

## 4. Wrap each resulting unit in `checks/<name>/{run.py, Dockerfile}`

Directory layout, one pair per check, sitting alongside the script's home
directory:

```
<pipeline_dir>/
  <existing modules, refactored per step 3 — unchanged otherwise>
  checks/
    <check_name>/
      run.py
      Dockerfile
```

**`run.py`** — a thin argparse CLI:
- Imports the real function from its canonical module (never reimplements
  logic inline).
- Common flags: an identifying flag or two (whatever the target script uses to
  key its findings — an entity/sample id), `--output PATH` (required).
- Writes results as a **JSON array to `--output`**, never to stdout — stdout
  stays free for the wrapped tool's own subprocess logs. In this repo, use
  `qc_common.write_findings_json(findings, Path(args.output))`.
- Any additional artifact a check must hand off to a later step (see step 5)
  gets its own explicit flag (e.g. `--ena-metadata-out`).

**`Dockerfile`**:
- Build context is the pipeline's home directory (e.g.
  `docker build -f checks/<name>/Dockerfile -t qc-<name> scrnaseq_qc_pipeline/`).
- `COPY` only the module(s) this check actually needs, plus the shared
  schema/common module — vendor it via `COPY`, not a separately published
  package, unless that shared module is large enough to justify one. This is
  what keeps one source of truth: editing the canonical module once updates
  both the in-process pipeline and every container that `COPY`s it on next
  build.
- **Base image**: prefer an existing community image
  (`quay.io/biocontainers/<tool>`, `ghcr.io/...`) for any check wrapping an
  external bioinformatics/CLI tool — check quay.io/biocontainers first. Use a
  minimal custom image (e.g. `python:3.11-slim` + pip installs) only for
  pure-language logic, or a tool with no existing community image.

**Exit codes**: `0` whenever the script wrote valid output — including a
`FAIL` finding, since that's the check correctly doing its job. Non-zero only
when `run.py`'s own control flow crashes before producing output (a bug in
the wrapper itself, not a check that correctly reported a problem). This is
what lets an orchestrator (a local driver, or later a workflow engine)
distinguish "ran fine and found something wrong" from "the check is broken."

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
production orchestrator.

- Mirror the original (pre-split) script's own control-flow branching, so its
  output can be diffed against that script's output as the validation step.
- Favor a real scripting language over a static tool like `docker-compose`
  whenever the pipeline has data-dependent branching or a fan-out that isn't
  known until an earlier step runs (e.g. "loop once per file discovered by
  step 1") — compose's static service graph fights that kind of pipeline.
- After each container invocation, read that step's output JSON back on the
  host to decide the next step's arguments — never inside a container.
- If the original script produces a combined report (a CSV, a summary), the
  driver should produce the same shape from the concatenated per-step results,
  so the two can be diffed directly.

## Worked example

`scrnaseq_qc_pipeline/` in this repo applies every step above:
`fastq_qc.py` and `matrix_qc.py` hold the refactored composable functions
(step 3), `scrnaseq_qc_pipeline/checks/*/` holds the per-check
`run.py`/`Dockerfile` pairs (step 4), and
`scrnaseq_qc_pipeline/driver/run_dockerized_pipeline.py` is the local
end-to-end driver (step 6). Read those as a concrete reference when applying
this procedure to a new script.
