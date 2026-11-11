"""Read an external test's metrics from the JSON object its container
writes to stdout.

Goes in py-dcqc at src/dcqc/tests/json_stdout_metrics.py. Follows the
compute_status() example in py-dcqc's ARCHITECTURE.md ("Creating the Test
Class"): the base class never reads metrics from the container, so an
external test's metrics stay {} unless its class overrides compute_status().
"""

import json

from dcqc.tests.base_test import ExternalTestMixin, TestStatus


def compute_status_with_json_metrics(test: ExternalTestMixin) -> TestStatus:
    """Interpret the exit code as usual, then load stdout into test.metrics.

    Metrics are read only for passed and failed tests; after an error,
    stdout can hold partial output. Stdout that is not one JSON object makes
    the status error, with stderr kept in status_reason.
    """
    outputs = test._find_process_outputs()
    status = test._interpret_process_outputs(outputs)
    if status not in (TestStatus.PASS, TestStatus.FAIL):
        return status

    text = outputs["std_out"].read_text().strip()
    if not text:
        return status
    try:
        metrics = json.loads(text)
    except json.JSONDecodeError:
        metrics = None
    if not isinstance(metrics, dict):
        reason = "Tool wrote output that is not a JSON object"
        std_err = outputs["std_err"].read_text()
        test.status_reason = f"{reason}\n{std_err}" if std_err else reason
        return TestStatus.ERROR
    test.metrics = metrics
    return status
