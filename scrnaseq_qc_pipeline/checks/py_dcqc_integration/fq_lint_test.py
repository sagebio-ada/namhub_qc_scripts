"""Goes in py-dcqc at src/dcqc/tests/fq_lint_test.py.

Image: scrnaseq_qc_pipeline/checks/py_dcqc/fq_lint in namhub_qc_scripts.
"""

from dcqc.target import PairedTarget
from dcqc.tests.base_test import ExternalBaseTest, Process, TestStatus, TestTier
from dcqc.tests.constants import FQ_LINT_CONTAINER
from dcqc.tests.json_stdout_metrics import compute_status_with_json_metrics


class FqLintTest(ExternalBaseTest[PairedTarget]):
    """Tests that each FASTQ file of a pair is structurally valid (fq lint).

    Each file is linted on its own; the paired read-name check is not done.
    Exit codes: 0 pass, 1 a file is malformed, 2 the check could not run
    (for example, a missing file).
    """

    tier = TestTier.INTERNAL_CONFORMANCE
    pass_code = 0
    fail_code = 1
    failure_reason_location = "std_err"

    def generate_process(self) -> Process:
        paths = [file.stage() for file in self.target.files]
        command_args = ["python3", "/app/run.py", *(f"'{path.name}'" for path in paths)]
        return Process(container=FQ_LINT_CONTAINER, command_args=command_args)

    def compute_status(self) -> TestStatus:
        return compute_status_with_json_metrics(self)
