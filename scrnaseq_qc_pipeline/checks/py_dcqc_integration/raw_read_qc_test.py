"""Goes in py-dcqc at src/dcqc/tests/raw_read_qc_test.py.

Image: scrnaseq_qc_pipeline/checks/py_dcqc/raw_read_qc in namhub_qc_scripts
(FastQC + Kraken2 + MultiQC in one container). It replaces a separate
FastQC test, so FastqSuite must not also run one.
"""

from dcqc.target import PairedTarget
from dcqc.tests.base_test import ExternalBaseTest, Process, TestStatus, TestTier
from dcqc.tests.constants import RAW_READ_QC_CONTAINER
from dcqc.tests.json_stdout_metrics import compute_status_with_json_metrics


class RawReadQcTest(ExternalBaseTest[PairedTarget]):
    """Runs FastQC and Kraken2 on each FASTQ file of a pair, plus MultiQC.

    Fails only when a FastQC module reports FAIL. Kraken2 species, R1/R2
    read counts and MultiQC's General Statistics are reported as metrics.
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
        return Process(container=RAW_READ_QC_CONTAINER, command_args=command_args)

    def compute_status(self) -> TestStatus:
        return compute_status_with_json_metrics(self)
