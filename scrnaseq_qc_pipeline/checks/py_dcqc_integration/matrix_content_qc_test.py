"""Goes in py-dcqc at src/dcqc/tests/matrix_content_qc_test.py.

Image: scrnaseq_qc_pipeline/checks/py_dcqc/matrix_content_qc in
namhub_qc_scripts.
"""

from dcqc.target import SingleTarget
from dcqc.tests.base_test import ExternalBaseTest, Process, TestStatus, TestTier
from dcqc.tests.constants import MATRIX_CONTENT_QC_CONTAINER
from dcqc.tests.json_stdout_metrics import compute_status_with_json_metrics


class MatrixContentQcTest(ExternalBaseTest[SingleTarget]):
    """Tests a single-cell count matrix packaged as one tar archive.

    The archive holds barcodes.tsv, features.tsv (or genes.tsv) and
    matrix.mtx, gzipped or not, at any depth, for example GEO's
    <sample>_filtered_feature_bc_matrix.tar.gz. Fails if the archive or a
    file in it can't be read, or if no barcode has any counts. Cell and gene
    counts, empty barcodes and the mitochondrial fraction are metrics.
    """

    tier = TestTier.INTERNAL_CONFORMANCE
    pass_code = 0
    fail_code = 1
    failure_reason_location = "std_err"

    def generate_process(self) -> Process:
        path = self.target.file.stage()
        command_args = ["python3", "/app/run.py", f"'{path.name}'"]
        return Process(container=MATRIX_CONTENT_QC_CONTAINER, command_args=command_args)

    def compute_status(self) -> TestStatus:
        return compute_status_with_json_metrics(self)
