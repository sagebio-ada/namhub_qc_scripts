"""Run metric sections per sample and write the TSV.

A section is a (name, KEYS, fn) tuple. fn(ctx, sample) returns a dict with exactly KEYS. A
missing or unparsable file leaves that section's columns blank and prints a warning.
"""
import csv
import sys


def collect(sample, ctx, sections):
    row = {"sample": sample}
    for name, keys, fn in sections:
        try:
            vals = fn(ctx, sample)
        except FileNotFoundError as e:
            print(f"WARNING [{sample}] {name}: missing {e.filename}", file=sys.stderr)
            vals = {}
        except (KeyError, ValueError, IndexError, StopIteration) as e:
            print(f"WARNING [{sample}] {name}: could not parse ({type(e).__name__}: {e})",
                  file=sys.stderr)
            vals = {}
        row.update({k: vals.get(k, "") for k in keys})
    return row


def write_tsv(rows, sections, out):
    """Write one row per sample: 'sample' then each section's KEYS in order."""
    fieldnames = ["sample"] + [k for _, keys, _ in sections for k in keys]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return fieldnames
