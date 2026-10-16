"""Helpers shared by atacseq_qc_metrics.py (the combined table) and atacseq_tools_qc.py (fq lint, kraken2/Bracken)."""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qc_common.util import read_table  # noqa: E402


def read_samplesheet(path):
    """nf-core/atacseq samplesheet -> {library: [(fastq_1, fastq_2 or None), ...]} in sheet order.

    library = <sample>_REP<replicate>; several rows of one library are technical replicates,
    run n (1-based, sheet order) of which nf-core calls <library>_T<n>.
    """
    libs = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = {"sample", "fastq_1", "replicate"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"{path}: missing column(s) {', '.join(sorted(missing))}")
        for r in reader:
            lib = f"{r['sample'].strip()}_REP{r['replicate'].strip()}"
            libs.setdefault(lib, []).append((r["fastq_1"].strip(), (r.get("fastq_2") or "").strip() or None))
    if not libs:
        sys.exit(f"No samples found in {path}")
    return libs


def select_runs(libs, requested=None):
    """(libraries, runs): the requested libraries (default all) and their sequencing runs,
    named as nf-core does: {<library>_T<n>: (library, fastq_1, fastq_2)}."""
    samples = list(libs)
    if requested:
        unknown = set(requested) - set(samples)
        if unknown:
            sys.exit(f"Samples not found: {', '.join(sorted(unknown))}")
        samples = list(requested)
    runs = {f"{lib}_T{n}": (lib, r1, r2)
            for lib in samples for n, (r1, r2) in enumerate(libs[lib], 1)}
    return samples, runs


def fastqc_tables(mqc_dir):
    """{'raw': {id: row}, 'trimmed': {id: row}} from multiqc_fastqc*.txt (raw vs trimmed by file name)."""
    out = {"raw": {}, "trimmed": {}}
    for path in sorted(mqc_dir.glob("multiqc_fastqc*.txt")):
        rows = read_table(path)
        if not rows or "Filename" not in rows[0]:
            continue
        stage = "trimmed" if re.search(r"_val_\d|_trimmed", rows[0]["Filename"]) else "raw"
        out[stage].update({r["Sample"]: r for r in rows})
    for stage, tbl in out.items():
        if not tbl:
            print(f"WARNING: no {stage} FastQC table in {mqc_dir}", file=sys.stderr)
    return out
