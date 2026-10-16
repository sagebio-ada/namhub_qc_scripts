"""Median read length from a FastQC zip (for MultiQC tables that lack median_sequence_length)."""
import re
import zipfile


def median_read_length(zip_path):
    """Median of the 'Sequence Length Distribution' module of fastqc_data.txt.

    FastQC bins long reads into ranges ("100-104"); a bin counts at its midpoint.
    """
    with zipfile.ZipFile(zip_path) as z:
        name = next(n for n in z.namelist() if n.endswith("/fastqc_data.txt"))
        text = z.read(name).decode()
    block = re.search(r">>Sequence Length Distribution\t[^\n]*\n(.*?)>>END_MODULE", text, re.S)
    if not block:
        raise ValueError(f"no Sequence Length Distribution in {zip_path}")
    bins = []
    for ln in block[1].splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        length, count = ln.split("\t")
        lo, _, hi = length.partition("-")
        bins.append(((float(lo) + float(hi or lo)) / 2, float(count)))
    bins.sort()
    half = sum(c for _, c in bins) / 2
    seen = 0
    for length, count in bins:
        seen += count
        if seen >= half:
            return int(length) if length == int(length) else length
