"""Goes in py-dcqc at src/dcqc/tests/fastqc_test.py.

Image: scrnaseq_qc_pipeline/checks/py_dcqc/fastqc in namhub_qc_scripts
(FastQC + MultiQC in one container).
"""

from dcqc.target import PairedTarget
from dcqc.tests.base_test import ExternalBaseTest, Process, TestStatus, TestTier
from dcqc.tests.constants import FASTQC_CONTAINER
from dcqc.tests.json_stdout_metrics import compute_status_with_json_metrics


class FastqcTest(ExternalBaseTest[PairedTarget]):
    """Runs FastQC on each FASTQ file of a pair, then MultiQC on that output.

    Fails when a FastQC module reports FAIL or FastQC can't process a file
    (for example, a truncated gzip). FastQC's basic statistics, each
    module's result, the R1/R2 read counts and everything in the MultiQC
    report (tables, plot data and each section's notes) are reported as
    metrics, since nf-dcqc doesn't keep the HTML; a MultiQC error is a
    metric too and never changes the result.
    The ENA registry comparison is skipped: it needs --ena-metadata-json,
    which py-dcqc cannot supply.
    """

    tier = TestTier.SUBJECTIVE_CONFORMANCE
    pass_code = 0
    fail_code = 1
    failure_reason_location = "std_err"

    def generate_process(self) -> Process:
        paths = [file.stage() for file in self.target.files]
        command_args = ["python3", "/app/run.py", *(f"'{path.name}'" for path in paths)]
        return Process(container=FASTQC_CONTAINER, command_args=command_args)

    def compute_status(self) -> TestStatus:
        return compute_status_with_json_metrics(self)
