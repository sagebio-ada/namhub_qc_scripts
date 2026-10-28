"""Classify a raw sequencing run's mate layout: paired (separate R1/R2
files), interleaved (one file, mates alternating within it), or single-end
(one file, no second mate present at all).

Pipeline-agnostic: works on any resolved fastq(.gz) file path(s), regardless
of where they came from (SRA/ENA, a direct upload, ...) or which pipeline in
this repo is calling it. Returns a plain FastqLayout result, not any one
pipeline's own Finding/report type -- each pipeline wraps this into its own
reporting schema (see scrnaseq_qc_pipeline/fastq_qc.py for an example).

Detection approach:
  - 2+ files resolved for one run -> paired. Mates are already split across
    separate files; nothing needs to be read to know that.
  - Exactly 1 file -> inspect its content, since ENA/a single upload can't
    tell you whether that lone file is truly interleaved (mates alternating
    within it) or single-end (only one mate ever archived). Samples the
    first `sample_records` reads and looks for two signals, in order:
      1. Read-header mate suffix (legacy "/1"/"/2", or Illumina CASAVA
         1.8+'s "<id> 1:N:0:..." / "<id> 2:N:0:..." second token) --
         alternating 1/2 every read means interleaved; a constant suffix
         throughout means single-end.
      2. If headers carry no usable mate suffix, fall back to read-length
         alternation (e.g. a 26bp barcode read alternating with a 91bp
         cDNA read is the classic 10x interleaved pattern) -- a weaker
         signal, since two genuinely single-end runs of different lengths
         concatenated would produce the same pattern, but useful when
         headers alone are inconclusive.
  - 0 files -> ambiguous (nothing to inspect).

This is a heuristic, not an authoritative source -- unusual header
conventions can fool it, and it can't match the certainty of directly
inspecting an SRA accession's own spot structure (e.g. via `vdb-dump`) for
a run that's genuinely SRA-sourced. Its advantage is that it works on any
fastq file already on disk, independent of where it came from.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

_MATE_SUFFIX_RE = re.compile(r"/([12])$")


@dataclass
class FastqLayout:
    layout: str  # "paired" | "interleaved" | "single_end" | "ambiguous"
    detail: str  # human-readable explanation of how this was determined


def _read_mate_number(read_id: str) -> Optional[str]:
    """Extract the mate number ('1'/'2') from a fastq read id (without the
    leading '@'), supporting the legacy '.../1'/'.../2' suffix and Illumina
    CASAVA 1.8+'s '<id> 1:N:0:...'/'<id> 2:N:0:...' second-token form."""
    first_token, _, rest = read_id.partition(" ")
    m = _MATE_SUFFIX_RE.search(first_token)
    if m:
        return m.group(1)
    if rest[:1] in ("1", "2") and (len(rest) == 1 or rest[1] == ":"):
        return rest[0]
    return None


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.name.endswith(".gz") else open(path, "rt")


def _sniff_single_file(path: Path, sample_records: int = 2000) -> FastqLayout:
    mate_numbers: List[Optional[str]] = []
    lengths: List[int] = []
    with _open_text(path) as fh:
        for i, line in enumerate(fh):
            if i >= sample_records * 4:
                break
            line = line.rstrip("\n")
            if i % 4 == 0:
                mate_numbers.append(_read_mate_number(line[1:]) if line.startswith("@") else None)
            elif i % 4 == 1:
                lengths.append(len(line))

    n = len(lengths)

    if mate_numbers and all(m is not None for m in mate_numbers):
        distinct_mates = set(mate_numbers)
        alternates = all(mate_numbers[i] != mate_numbers[i + 1] for i in range(len(mate_numbers) - 1))
        if distinct_mates == {"1", "2"} and alternates:
            return FastqLayout("interleaved", f"{n} reads sampled: header mate suffix alternates 1/2 every read")
        if len(distinct_mates) == 1:
            (only,) = distinct_mates
            return FastqLayout(
                "single_end", f"{n} reads sampled: every header carries mate suffix {only!r}, no alternation",
            )

    distinct_lengths = sorted(set(lengths))
    if len(distinct_lengths) == 2:
        a, b = distinct_lengths
        alternates_len = all((lengths[i] == a) != (lengths[i + 1] == a) for i in range(len(lengths) - 1))
        if alternates_len:
            return FastqLayout(
                "interleaved",
                f"{n} reads sampled: read length alternates {a}bp/{b}bp every read "
                f"(no usable mate suffix in headers)",
            )

    if len(distinct_lengths) == 1:
        return FastqLayout(
            "single_end", f"{n} reads sampled: uniform {distinct_lengths[0]}bp, no alternating mate pattern",
        )

    return FastqLayout(
        "ambiguous",
        f"{n} reads sampled: {len(distinct_lengths)} distinct length(s) {distinct_lengths}, "
        f"no clear alternating pattern",
    )


def detect_fastq_layout(file_paths: Sequence[Path]) -> FastqLayout:
    """Classify the mate layout of one raw run given its resolved fastq(.gz)
    file path(s). 2+ paths -> paired, no need to open them. 1 path -> sniff
    its content. 0 paths -> ambiguous."""
    paths = [Path(p) for p in file_paths]
    if len(paths) >= 2:
        names = ", ".join(p.name for p in paths)
        return FastqLayout("paired", f"paired across {len(paths)} separate files: {names}")
    if len(paths) == 1:
        return _sniff_single_file(paths[0])
    return FastqLayout("ambiguous", "no files given to inspect")
