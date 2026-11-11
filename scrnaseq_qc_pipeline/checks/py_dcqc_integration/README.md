# py-dcqc test classes for the FASTQ checks

The py-dcqc side of `checks/py_dcqc/fq_lint` and `checks/py_dcqc/raw_read_qc`.
These files are written to go into
[py-dcqc](https://github.com/Sage-Bionetworks-Workflows/py-dcqc) unchanged.

## Apply

`py-dcqc.patch` holds every change below. It applies cleanly to py-dcqc's
`DPE-1842` branch (PR #99, the branch that makes `PairedTarget` work):

```
cd py-dcqc
git apply /path/to/namhub_qc_scripts/scrnaseq_qc_pipeline/checks/py_dcqc_integration/py-dcqc.patch
```

| File here | Goes to | What it does |
|---|---|---|
| `json_stdout_metrics.py` | `src/dcqc/tests/` | `compute_status_with_json_metrics()`: the `compute_status()` override from ARCHITECTURE.md ("Creating the Test Class"). It loads the container's stdout JSON into `metrics`, so the metrics reach `suites.json` and the `dcqc_metrics` CSV column. |
| `fq_lint_test.py` | `src/dcqc/tests/` | `FqLintTest` runs `python3 /app/run.py R1 R2`. Each file is linted on its own, and the paired read-name check is not done. |
| `raw_read_qc_test.py` | `src/dcqc/tests/` | `RawReadQcTest` runs FastQC, Kraken2 and MultiQC on the pair. |

The patch also:

- imports both classes in `src/dcqc/tests/__init__.py`,
- adds `FQ_LINT_CONTAINER` and `RAW_READ_QC_CONTAINER` to `src/dcqc/tests/constants.py`,
- sets `FastqSuite.add_tests` to `(PairedFastqParityTest, FqLintTest, RawReadQcTest)`.
  There is **no separate `FastqcTest`**. `RawReadQcTest` already runs FastQC,
  so a separate test would report the same FastQC FAIL two times.

### Unit tests in the patch

- `tests/test_external_tests.py`: `TestJsonMetricsFastqTests` runs for both
  classes. It checks the container and command, and the result for each
  case:
  - exit 0: pass, metrics read;
  - exit 1: fail, metrics and the stderr reason read;
  - exit 2: error, metrics not read;
  - stdout that is not JSON, or JSON that is not an object: error;
  - empty stdout: pass with no metrics.
- `tests/test_suites.py`: checks that `FastqSuite.add_tests` is exactly
  `(PairedFastqParityTest, FqLintTest, RawReadQcTest)`.
- `tests/test_parsers.py`: the existing "grouped FASTQ manifest gives a suite
  that runs" test now skips the external tests, as `run_pipeline` in
  `tests/test_main.py` already does. Without the skip, it would try to read
  container outputs that only nf-dcqc makes.

No Docker exit-code tests (`@docker_enabled_test`) were added. They need the
images to be published first.

## Before merging

- **Image names are placeholders.** The `ghcr.io/sage-bionetworks-workflows/qc-*:0.1.0`
  values in `constants.py` must match where the images are published.
- **Run py-dcqc's test suite on the patch.** The unit tests below were
  written but not run (see "How it was checked").
- `matrix_content_qc` has no class here. It needs three files
  (barcodes/features/matrix), which no py-dcqc target supports yet.

## How it was checked

- Both images were run with each class's command, the way nf-dcqc runs it:
  `( (cmd) > std_out.txt 2> std_err.txt; echo $? > exit_code.txt ) || true`.
- `compute_status_with_json_metrics()` was run against those outputs. Metrics
  were loaded for pass and fail. Stdout that is not JSON gave `error`, with
  stderr kept. Metrics were not read after an error exit.
- py-dcqc's own test suite, including the new unit tests, was **not** run
  against the patch: its dependencies could not be installed in the
  environment where this was written. All changed files compile, and the
  patch applies cleanly to `DPE-1842`.
