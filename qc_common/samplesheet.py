"""nf-core samplesheet handling."""
import csv
import sys


def read_samples(path):
    """Unique values of the 'sample' column of an nf-core samplesheet CSV, in order of first
    appearance (replicate rows share a sample). Exits with a message if there are none."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if "sample" not in (reader.fieldnames or []):
            sys.exit(f"No 'sample' column in {path}")
        samples = list(dict.fromkeys(r["sample"].strip() for r in reader if r["sample"].strip()))
    if not samples:
        sys.exit(f"No samples found in {path}")
    return samples
