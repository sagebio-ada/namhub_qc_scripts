"""Goes in py-dcqc at src/dcqc/tests/kraken2_qc_test.py.

Image: scrnaseq_qc_pipeline/checks/py_dcqc/kraken2_qc in namhub_qc_scripts
(Kraken2 with its Standard-8 database baked in).
"""

from dcqc.target import PairedTarget
from dcqc.tests.base_test import ExternalBaseTest, Process, TestStatus, TestTier
from dcqc.tests.constants import KRAKEN2_QC_CONTAINER
from dcqc.tests.json_stdout_metrics import compute_status_with_json_metrics


class Kraken2QcTest(ExternalBaseTest[PairedTarget]):
    """Classifies the reads of each FASTQ file of a pair with Kraken2.

    Reports only: the read count, the classified and unclassified
    percentages and the top species are metrics, and the species never
    makes the test fail. It fails only when Kraken2 can't process a file
    (for example, a truncated gzip).
    """

    tier = TestTier.SUBJECTIVE_CONFORMANCE
    pass_code = 0
    fail_code = 1
    failure_reason_location = "std_err"

    def generate_process(self) -> Process:
        paths = [file.stage() for file in self.target.files]
        command_args = ["python3", "/app/run.py", *(f"'{path.name}'" for path in paths)]
        return Process(container=KRAKEN2_QC_CONTAINER, command_args=command_args)

    def compute_status(self) -> TestStatus:
        return compute_status_with_json_metrics(self)
