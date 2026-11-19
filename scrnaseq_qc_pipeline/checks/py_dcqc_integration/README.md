# py-dcqc test classes for the scRNA-seq checks

The py-dcqc side of `checks/py_dcqc/fq_lint`, `checks/py_dcqc/fastqc`,
`checks/py_dcqc/kraken2_qc` and `checks/py_dcqc/matrix_content_qc`.
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
| `fastqc_test.py` | `src/dcqc/tests/` | `FastqcTest` runs FastQC on the pair, then MultiQC on that output. Fails on a FastQC module FAIL. |
| `kraken2_qc_test.py` | `src/dcqc/tests/` | `Kraken2QcTest` runs Kraken2 on the pair. Metrics only: species never makes it fail. |
| `matrix_content_qc_test.py` | `src/dcqc/tests/` | `MatrixContentQcTest` runs `python3 /app/run.py ARCHIVE` on one matrix archive. |

The patch also:

- imports the four classes in `src/dcqc/tests/__init__.py`,
- adds `FQ_LINT_CONTAINER`, `FASTQC_CONTAINER`, `KRAKEN2_QC_CONTAINER` and
  `MATRIX_CONTENT_QC_CONTAINER` to `src/dcqc/tests/constants.py`,
- adds a `CELLRANGER-MATRIX-TAR` file type to `src/dcqc/file.py`: extensions `.tar.gz`,
  `.tgz` and `.tar`, EDAM `format_3981` ("TAR format"),
- adds `CellrangerMatrixTarSuite` for `CELLRANGER-MATRIX-TAR`, with `add_tests = (MatrixContentQcTest,)`
  and the default `SingleTarget`,
- sets `FastqSuite.add_tests` to
  `(PairedFastqParityTest, FqLintTest, FastqcTest, Kraken2QcTest)`. Each
  tool runs in one test only, so a FastQC FAIL is reported once.

### Unit tests in the patch

- `tests/test_external_tests.py`: `TestJsonMetricsTests` runs for all four
  classes. It checks the container and command, and the result for each
  case:
  - exit 0: pass, metrics read;
  - exit 1: fail, metrics and the stderr reason read;
  - exit 2: error, metrics not read;
  - stdout that is not JSON, or JSON that is not an object: error;
  - empty stdout: pass with no metrics.
- `tests/test_suites.py`: checks that `FastqSuite.add_tests` is exactly
  `(PairedFastqParityTest, FqLintTest, FastqcTest, Kraken2QcTest)`, and that a
  `CELLRANGER-MATRIX-TAR` file gets `CellrangerMatrixTarSuite`.
- `tests/test_parsers.py`: the existing "grouped FASTQ manifest gives a suite
  that runs" test now skips the external tests, as `run_pipeline` in
  `tests/test_main.py` already does. Without the skip, it would try to read
  container outputs that only nf-dcqc makes.

No Docker exit-code tests (`@docker_enabled_test`) were added. They need the
images to be published first.

## Before merging

- **File type name.** A manifest's file type column must use
  `CELLRANGER-MATRIX-TAR` for a matrix archive.
- **Image names are placeholders.** The `ghcr.io/sage-bionetworks-workflows/qc-*:0.1.0`
  values in `constants.py` must match where the images are published.
- **Run py-dcqc's test suite on the patch.** The unit tests above were
  written but not run (see "How it was checked").

## Matrix archives

`MatrixContentQcTest` takes one tar archive (`.tar.gz`, `.tgz` or `.tar`)
that holds `barcodes.tsv`, `features.tsv` (or CellRanger v2's `genes.tsv`)
and `matrix.mtx`, each gzipped or not. The files are found by name at any
depth, so GEO's `<sample>_filtered_feature_bc_matrix.tar.gz` works as-is.

The test never builds archives. That is a step before py-dcqc:
`scrnaseq_qc_pipeline/package_matrix.py` finds every barcodes/features/matrix
set among the file entities under a Synapse project or folder, packages each
set, and stores each archive in Synapse as its own file entity. Each archive
keeps the annotations its three source files share (for example `SampleId`),
and its provenance points at those three files:

```
SYNAPSE_AUTH_TOKEN=... python3 package_matrix.py --synapse-id syn123 \
    -o archives/ --upload-to syn456
```

The new entities are what the py-dcqc manifest lists, with file type
`CELLRANGER-MATRIX-TAR`. Without `--upload-to`, the archives are only
written to `archives/`. `--scan DIR` does the same for files already on
disk.

Exit codes: 0 pass, 1 when the archive can't be read, a member is missing
or duplicated (for example, a raw and a filtered matrix in one archive), a
file can't be parsed, or no barcode has any counts, and 2 when the input
file does not exist.

## How it was checked

- All four images were run with each class's command, the way nf-dcqc runs it:
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
