"""Raw-read QC: resolve real FASTQ files behind Synapse's SRA/ENA links, then
QC the actual external data against its own external registry record.

Many ADA/MC2 projects store raw scRNA-seq reads as Synapse ExternalFileHandle
entities that just point at an SRA run on NCBI's public S3 bucket (the run
itself was never uploaded to Synapse). FastQC needs an actual FASTQ file, so
this module resolves the SRA run accession to a real, directly-downloadable
FASTQ via the ENA file report API (ENA mirrors SRA and re-encodes runs as
plain fastq.gz, unlike NCBI's .sra format) using the `geo_synapse.ena` helpers
already maintained in the geo_dataset_creation repo, downloads it, and checks
it two ways:

  1. Integrity: does the downloaded file's size/md5/read-count match what
     ENA's own file report registered for this run?
  2. Content quality: FastQC (per-base quality, GC content, adapter
     contamination, duplication, overrepresented sequences).

This module deliberately never compares against Synapse's curated
annotations — that's a separate, GEO-provenance-based check (see
geo_metadata_qc.py). Everything here is external data checked against its
own external registry record.
"""

from __future__ import annotations

import gzip
import hashlib
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

import fq_lint_qc
import kraken_qc
from qc_common import Finding, make_finding

try:
    from geo_synapse import ena as _ena
    HAS_GEO_SYNAPSE = True
except ImportError:
    HAS_GEO_SYNAPSE = False

GEO_SYNAPSE_INSTALL_HINT = (
    "geo-synapse is required to resolve SRA/ENA fastq download links. "
    "pip install git+https://github.com/sagebio-ada/geo_dataset_creation.git"
)

FASTQC_INSTALL_HINT = (
    "fastqc was not found on PATH. Install it with "
    "'brew install fastqc' or 'conda install -c bioconda fastqc'."
)

_DIGITS_RE = re.compile(r"^\d+$")


def fastqc_available() -> bool:
    return shutil.which("fastqc") is not None


def build_session() -> requests.Session:
    retry = Retry(
        total=3, connect=3, read=3, status=3, backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET"]), raise_on_status=False,
    )
    session = requests.Session()
    session.headers.update({"User-Agent": "ada-scrnaseq-qc/1.0"})
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def download_fastq(session: requests.Session, url: str, dest: Path, max_bytes: Optional[int] = None) -> int:
    """Stream a fastq(.gz) to disk. If max_bytes is set, truncates early for spot QC."""
    written = 0
    with session.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
                written += len(chunk)
                if max_bytes and written >= max_bytes:
                    break
    return written


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sanitize_capped_download(path: Path) -> None:
    """A byte-capped download ends mid-gzip-block, which FastQC refuses outright
    ("Ran out of data in the middle of a fastq entry"). Decode whatever complete
    lines came through before the truncation, drop any partial trailing FASTQ
    record (records are 4 lines), and rewrite as a small but valid fastq.gz.
    """
    lines: List[str] = []
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                lines.append(line)
    except (EOFError, OSError):
        pass  # expected: we cut the stream off mid-block
    complete = len(lines) - (len(lines) % 4)
    with gzip.open(path, "wt") as fh:
        fh.writelines(lines[:complete])


def run_fastqc(fastq_path: Path, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["fastqc", "--extract", "-o", str(outdir), str(fastq_path)],
        check=True, capture_output=True, text=True,
    )
    stem = fastq_path.name
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    data_txt = outdir / f"{stem}_fastqc" / "fastqc_data.txt"
    if not data_txt.exists():
        raise FileNotFoundError(f"Expected FastQC output not found: {data_txt}")
    return data_txt


def parse_fastqc_data(path: Path) -> Dict[str, object]:
    """Parse fastqc_data.txt into basic stats (dict) + per-module PASS/WARN/FAIL."""
    basic_stats: Dict[str, str] = {}
    modules: Dict[str, str] = {}
    current_module: Optional[str] = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">>END_MODULE"):
                current_module = None
                continue
            if line.startswith(">>"):
                name, status = line[2:].rsplit("\t", 1)
                current_module = name
                modules[name] = status
                continue
            if current_module == "Basic Statistics" and "\t" in line and not line.startswith("#"):
                key, val = line.split("\t", 1)
                basic_stats[key] = val
    return {"basic_statistics": basic_stats, "modules": modules}


def qc_raw_run(
    entity_id: str,
    entity_name: str,
    srr_accession: str,
    work_dir: Path,
    max_download_mb: Optional[float] = None,
    keep_files: bool = False,
    kraken2_db: Optional[str] = None,
    expected_species: Optional[str] = None,
    kraken2_threads: int = 4,
) -> List[Finding]:
    findings: List[Finding] = []

    def add(check: str, status: str, detail: str) -> None:
        findings.append(make_finding(entity_id, entity_name, check, status, detail))

    if not HAS_GEO_SYNAPSE:
        add("fastq_resolution", "FAIL", GEO_SYNAPSE_INSTALL_HINT)
        return findings
    if not fastqc_available():
        add("fastqc", "FAIL", FASTQC_INSTALL_HINT)
        return findings

    session = build_session()
    try:
        links = _ena.run_to_flat_links(session, srr_accession)
    except Exception as exc:
        add("fastq_resolution", "FAIL", f"Could not resolve ENA fastq links for {srr_accession}: {exc}")
        return findings
    if not links:
        add("fastq_resolution", "FAIL", f"No ENA fastq links found for {srr_accession} (run may not be public yet)")
        return findings
    add("fastq_resolution", "PASS", f"{len(links)} fastq file(s) resolved from ENA for {srr_accession}")

    run_dir = work_dir / entity_name
    run_dir.mkdir(parents=True, exist_ok=True)
    max_bytes = int(max_download_mb * 1024 * 1024) if max_download_mb else None

    observed_lengths: List[str] = []
    total_fastqc_seqs = 0
    any_capped = False
    per_file_seq_counts: Dict[str, int] = {}

    for item in links:
        url = item["link"]
        fname = url.rsplit("/", 1)[-1]
        dest = run_dir / fname
        registered_size_str = str(item.get("file_size", ""))
        already_cached = (
            not max_bytes and dest.exists() and _DIGITS_RE.match(registered_size_str)
            and dest.stat().st_size == int(registered_size_str)
        )
        if already_cached:
            size = dest.stat().st_size
            add("fastq_download", "PASS", f"{fname}: already downloaded ({size / 1e6:.1f} MB), reusing")
        else:
            try:
                size = download_fastq(session, url, dest, max_bytes=max_bytes)
            except Exception as exc:
                add("fastq_download", "FAIL", f"{fname}: {exc}")
                continue
            was_capped_now = bool(max_bytes and size >= max_bytes)
            truncated = " (capped for spot-QC, trimmed to whole records)" if was_capped_now else ""
            add("fastq_download", "PASS", f"{fname}: downloaded {size / 1e6:.1f} MB{truncated}")
        was_capped = bool(max_bytes and size >= max_bytes)
        any_capped = any_capped or was_capped

        if was_capped:
            try:
                sanitize_capped_download(dest)
            except Exception as exc:
                add("fastq_download", "FAIL", f"{fname}: could not sanitize capped download: {exc}")
                if not keep_files:
                    dest.unlink(missing_ok=True)
                continue
        else:
            # Integrity: does the downloaded byte stream match ENA's own
            # registered size/checksum for this file? Only meaningful for a
            # full (uncapped) download.
            registered_size = str(item.get("file_size", ""))
            if _DIGITS_RE.match(registered_size):
                status = "PASS" if int(registered_size) == size else "FAIL"
                add("fastq_size_vs_ena", status,
                    f"{fname}: downloaded {size} bytes vs ENA registered {registered_size} bytes")
            registered_md5 = item.get("md5", "")
            if registered_md5:
                actual_md5 = md5_file(dest)
                status = "PASS" if actual_md5 == registered_md5 else "FAIL"
                add("fastq_md5_vs_ena", status,
                    f"{fname}: md5={actual_md5} vs ENA registered md5={registered_md5}")

        findings.extend(fq_lint_qc.lint_fastq(entity_id, entity_name, dest))

        if kraken2_db:
            findings.extend(kraken_qc.qc_kraken2(
                entity_id, entity_name, dest, kraken2_db, run_dir / "kraken2",
                expected_species=expected_species, threads=kraken2_threads,
            ))

        try:
            data_txt = run_fastqc(dest, run_dir / "fastqc")
        except Exception as exc:
            add("fastqc", "FAIL", f"{fname}: FastQC failed: {exc}")
            if not keep_files:
                dest.unlink(missing_ok=True)
            continue

        parsed = parse_fastqc_data(data_txt)
        basic = parsed["basic_statistics"]
        modules = parsed["modules"]
        observed_lengths.append(basic.get("Sequence length", ""))
        total_seqs_str = basic.get("Total Sequences", "").replace(",", "")
        if _DIGITS_RE.match(total_seqs_str):
            total_fastqc_seqs += int(total_seqs_str)
            per_file_seq_counts[fname] = int(total_seqs_str)
        add("fastqc_basic_stats", "INFO",
            f"{fname}: {basic.get('Total Sequences', '?')} seqs, "
            f"len={basic.get('Sequence length', '?')}, GC={basic.get('%GC', '?')}%")
        for module_name, mod_status in modules.items():
            sev = {"pass": "PASS", "warn": "WARN", "fail": "FAIL"}.get(mod_status.lower(), "WARN")
            if sev != "PASS":
                add(f"fastqc_module:{module_name}", sev, f"{fname}: {mod_status}")

        if not keep_files:
            dest.unlink(missing_ok=True)

    if not any_capped and total_fastqc_seqs and links:
        registered_read_count = str(links[0].get("read_count", ""))
        if _DIGITS_RE.match(registered_read_count):
            status = "PASS" if total_fastqc_seqs == int(registered_read_count) else "WARN"
            add("read_count_vs_ena", status,
                f"FastQC total sequences (summed across {len(links)} file(s))={total_fastqc_seqs} "
                f"vs ENA registered read_count={registered_read_count}")

    if not any_capped and len(per_file_seq_counts) >= 2:
        # Split R1/R2 mate files (as opposed to a single interleaved fastq)
        # should have exactly the same number of reads. A mismatch means one
        # mate is truncated or the pair was assembled incorrectly.
        distinct_counts = set(per_file_seq_counts.values())
        if len(distinct_counts) == 1:
            add("paired_fastq_parity", "PASS",
                f"All {len(per_file_seq_counts)} mate file(s) have matching read counts "
                f"({distinct_counts.pop()})")
        else:
            detail = ", ".join(f"{fname}={count}" for fname, count in per_file_seq_counts.items())
            add("paired_fastq_parity", "FAIL", f"Mate files have mismatched read counts: {detail}")

    if not any_capped and observed_lengths and links:
        base_count = str(links[0].get("base_count", ""))
        read_count = str(links[0].get("read_count", ""))
        if _DIGITS_RE.match(base_count) and _DIGITS_RE.match(read_count) and int(read_count) > 0:
            mean_len = int(base_count) / int(read_count)
            # No pass/fail judgment -- just report both numbers side by side.
            add("mean_length_vs_registry", "INFO",
                f"ENA base_count/read_count implies mean length={mean_len:.1f}, "
                f"FastQC observed length(s)={observed_lengths}")

    return findings


def qc_direct_fastq(syn, entity_id: str, entity_name: str,
                     work_dir: Path, keep_files: bool = False,
                     kraken2_db: Optional[str] = None, expected_species: Optional[str] = None,
                     kraken2_threads: int = 4) -> List[Finding]:
    """QC for a FASTQ that was actually uploaded to Synapse (no SRA indirection).

    There's no external registry to check integrity against here (it's a
    direct Synapse upload, not a link to SRA/ENA) — this is content QC only.
    """
    findings: List[Finding] = []

    def add(check: str, status: str, detail: str) -> None:
        findings.append(make_finding(entity_id, entity_name, check, status, detail))

    if not fastqc_available():
        add("fastqc", "FAIL", FASTQC_INSTALL_HINT)
        return findings

    run_dir = work_dir / entity_name
    run_dir.mkdir(parents=True, exist_ok=True)
    try:
        ent = syn.get(entity_id, downloadLocation=str(run_dir), ifcollision="overwrite.local")
    except Exception as exc:
        add("fastq_download", "FAIL", f"Could not download from Synapse: {exc}")
        return findings
    add("fastq_download", "PASS", f"Downloaded {ent.name} from Synapse")

    findings.extend(fq_lint_qc.lint_fastq(entity_id, entity_name, Path(ent.path)))

    if kraken2_db:
        findings.extend(kraken_qc.qc_kraken2(
            entity_id, entity_name, Path(ent.path), kraken2_db, run_dir / "kraken2",
            expected_species=expected_species, threads=kraken2_threads,
        ))

    try:
        data_txt = run_fastqc(Path(ent.path), run_dir / "fastqc")
    except Exception as exc:
        add("fastqc", "FAIL", f"FastQC failed: {exc}")
        return findings

    parsed = parse_fastqc_data(data_txt)
    basic = parsed["basic_statistics"]
    modules = parsed["modules"]
    add("fastqc_basic_stats", "INFO",
        f"{basic.get('Total Sequences', '?')} seqs, "
        f"len={basic.get('Sequence length', '?')}, GC={basic.get('%GC', '?')}%")
    for module_name, mod_status in modules.items():
        sev = {"pass": "PASS", "warn": "WARN", "fail": "FAIL"}.get(mod_status.lower(), "WARN")
        if sev != "PASS":
            add(f"fastqc_module:{module_name}", sev, mod_status)

    if not keep_files:
        Path(ent.path).unlink(missing_ok=True)
    return findings
