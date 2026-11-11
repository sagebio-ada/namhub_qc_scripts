# py-dcqc test classes for the scRNA-seq checks

The py-dcqc side of `checks/py_dcqc/fq_lint`, `checks/py_dcqc/raw_read_qc`
and `checks/py_dcqc/matrix_content_qc`.
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
| `matrix_content_qc_test.py` | `src/dcqc/tests/` | `MatrixContentQcTest` runs `python3 /app/run.py ARCHIVE` on one matrix archive. |

The patch also:

- imports the three classes in `src/dcqc/tests/__init__.py`,
- adds `FQ_LINT_CONTAINER`, `RAW_READ_QC_CONTAINER` and
  `MATRIX_CONTENT_QC_CONTAINER` to `src/dcqc/tests/constants.py`,
- adds a `10X-MATRIX` file type to `src/dcqc/file.py`: extensions `.tar.gz`,
  `.tgz` and `.tar`, EDAM `format_3981` ("TAR format"),
- adds `MatrixSuite` for `10X-MATRIX`, with `add_tests = (MatrixContentQcTest,)`
  and the default `SingleTarget`,
- sets `FastqSuite.add_tests` to `(PairedFastqParityTest, FqLintTest, RawReadQcTest)`.
  There is **no separate `FastqcTest`**. `RawReadQcTest` already runs FastQC,
  so a separate test would report the same FastQC FAIL two times.

### Unit tests in the patch

- `tests/test_external_tests.py`: `TestJsonMetricsTests` runs for all three
  classes. It checks the container and command, and the result for each
  case:
  - exit 0: pass, metrics read;
  - exit 1: fail, metrics and the stderr reason read;
  - exit 2: error, metrics not read;
  - stdout that is not JSON, or JSON that is not an object: error;
  - empty stdout: pass with no metrics.
- `tests/test_suites.py`: checks that `FastqSuite.add_tests` is exactly
  `(PairedFastqParityTest, FqLintTest, RawReadQcTest)`, and that a
  `10X-MATRIX` file gets `MatrixSuite`.
- `tests/test_parsers.py`: the existing "grouped FASTQ manifest gives a suite
  that runs" test now skips the external tests, as `run_pipeline` in
  `tests/test_main.py` already does. Without the skip, it would try to read
  container outputs that only nf-dcqc makes.

No Docker exit-code tests (`@docker_enabled_test`) were added. They need the
images to be published first.

## Before merging

- **The `10X-MATRIX` name is a proposal.** It is the value a manifest's
  file type column must use for a matrix archive.
- **Image names are placeholders.** The `ghcr.io/sage-bionetworks-workflows/qc-*:0.1.0`
  values in `constants.py` must match where the images are published.
- **Run py-dcqc's test suite on the patch.** The unit tests above were
  written but not run (see "How it was checked").

## Matrix archives

`MatrixContentQcTest` takes one tar archive (`.tar.gz`, `.tgz` or `.tar`)
that holds `barcodes.tsv`, `features.tsv` (or CellRanger v2's `genes.tsv`)
and `matrix.mtx`, each gzipped or not. The files are found by name at any
depth, so GEO's `<sample>_filtered_feature_bc_matrix.tar.gz` works as-is.
For three loose files, `scrnaseq_qc_pipeline/package_matrix.py` builds the
archive:

```
python3 package_matrix.py GSM123_barcodes.tsv.gz GSM123_features.tsv.gz \
    GSM123_matrix.mtx.gz -o GSM123_matrix.tar.gz
```

The archive must exist as its own file, with its own Synapse ID, for
py-dcqc to test it.

Exit codes: 0 pass, 1 when the archive can't be read, a member is missing
or duplicated (for example, a raw and a filtered matrix in one archive), a
file can't be parsed, or no barcode has any counts, and 2 when the input
file does not exist.

## How it was checked

- All three images were run with each class's command, the way nf-dcqc runs it:
  `( (cmd) > std_out.txt 2> std_err.txt; echo $? > exit_code.txt ) || true`.
- The matrix image was run on: a packaged archive, a GEO-style archive with
  a top-level folder, an uncompressed CellRanger v2 archive, and three loose
  files (all exit 0); an all-zero matrix, an archive missing a member, one
  with duplicate members, and truncated and garbage archives (all exit 1);
  and a file that does not exist (exit 2).
- `compute_status_with_json_metrics()` was run against the FASTQ outputs. Metrics
  were loaded for pass and fail. Stdout that is not JSON gave `error`, with
  stderr kept. Metrics were not read after an error exit.
- py-dcqc's own test suite, including the new unit tests, was **not** run
  against the patch: its dependencies could not be installed in the
  environment where this was written. All changed files compile, and the
  patch applies cleanly to `DPE-1842`.
