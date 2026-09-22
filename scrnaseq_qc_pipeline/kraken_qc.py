"""Species/contamination screening via Kraken2 taxonomic classification.

Nothing else in this pipeline ever looks at whether the actual sequence
CONTENT matches the declared organism. GEO's own `sample_organism_ch1`
field is just as self-reported as Synapse's `species`/`referenceSet`
annotations, so cross-checking one against the other (geo_metadata_qc.py)
can't catch a mislabeled sample or cross-species contamination -- it would
just compare one unverified label to another. Kraken2 classifies reads
against a real reference database, giving an independent, content-based
read on what's actually in the file.

Kraken2 needs a pre-built database staged locally -- they run several GB.
`ensure_kraken2_db()` below will fetch one on demand (Kraken2's own
"Standard-8" build by default: https://benlangmead.github.io/aws-indexes/k2,
~5.5GB compressed; RefSeq archaea/bacteria/viral/plasmid/human/UniVec) if
nothing is already staged at the given path -- this is what lets the
dockerized `kraken2_qc` check run against a bare/empty mounted volume
without a human pre-downloading anything first. Point `--kraken2-db`/
`--db-path` at a *persistent* location (a host directory or named Docker
volume, not a container's own ephemeral filesystem) so the download only
happens once, not on every run.

Because classification runs at roughly a few million reads/minute even
multithreaded, running this against a full multi-hundred-million-read raw
file can take a long time -- pair it with --max-download-mb for routine QC
and only run it against the full file when you specifically need that.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

from qc_common import Finding, make_finding

KRAKEN2_INSTALL_HINT = (
    "kraken2 was not found on PATH. Install it with 'brew install kraken2', then "
    "download a reference database (e.g. Standard-8) from "
    "https://benlangmead.github.io/aws-indexes/k2 and pass its extracted "
    "directory via --kraken2-db."
)

_TOP_N_SPECIES = 5

# Kraken2's own prebuilt "Standard-8" index -- see the module docstring and
# https://benlangmead.github.io/aws-indexes/k2 for other prebuilt options.
_DEFAULT_DB_URL = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz"
_DB_MARKER_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")


def kraken2_available() -> bool:
    return shutil.which("kraken2") is not None


def kraken2_db_staged(db_path: Path) -> bool:
    """Whether db_path already has a usable Kraken2 database in it."""
    return all((db_path / f).exists() for f in _DB_MARKER_FILES)


def download_kraken2_db(db_path: Path, url: str = _DEFAULT_DB_URL) -> None:
    """Download and extract a Kraken2 database tarball into db_path.

    Stdlib-only (urllib/tarfile) so this also works inside the minimal
    kraken2_qc container without a new pip dependency just for this
    one-time fetch. Callers should check kraken2_db_staged() first --
    this doesn't skip an already-staged database itself.
    """
    import tarfile
    import urllib.request

    db_path.mkdir(parents=True, exist_ok=True)
    tarball = db_path / "_kraken2_db_download.tar.gz"
    with urllib.request.urlopen(url, timeout=60) as resp, open(tarball, "wb") as fh:
        shutil.copyfileobj(resp, fh, length=1 << 20)
    with tarfile.open(tarball) as tar:
        tar.extractall(db_path)
    tarball.unlink(missing_ok=True)


def ensure_kraken2_db(db_path: Path, url: Optional[str] = None) -> List[Finding]:
    """Make sure a usable Kraken2 database exists at db_path, downloading one
    (Kraken2's own Standard-8 build by default, or `url` if given) if it
    isn't already staged there -- what lets a kraken2_qc container run on
    demand against a bare/empty mounted volume instead of requiring a human
    to pre-stage the database. Emits one INFO/FAIL finding reporting what
    happened; callers should check kraken2_db_staged(db_path) before calling
    qc_kraken2() if this returns a FAIL.
    """
    if kraken2_db_staged(db_path):
        return [make_finding("kraken2_db", "kraken2_db", "kraken2_db", "INFO",
                              f"Database already staged at {db_path}")]
    try:
        download_kraken2_db(db_path, url=url or _DEFAULT_DB_URL)
    except Exception as exc:
        return [make_finding("kraken2_db", "kraken2_db", "kraken2_db", "FAIL",
                              f"Could not download Kraken2 database to {db_path}: {exc}")]
    if not kraken2_db_staged(db_path):
        return [make_finding("kraken2_db", "kraken2_db", "kraken2_db", "FAIL",
                              f"Downloaded to {db_path} but expected database files "
                              f"({', '.join(_DB_MARKER_FILES)}) are still missing")]
    return [make_finding("kraken2_db", "kraken2_db", "kraken2_db", "INFO",
                          f"Downloaded and extracted database to {db_path}")]


def run_kraken2(fastq_path: Path, db_path: str, outdir: Path, threads: int = 4) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    stem = fastq_path.name
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    report_path = outdir / f"{stem}.kraken2report.txt"
    cmd = [
        "kraken2", "--db", db_path, "--threads", str(threads),
        "--report", str(report_path), "--output", "/dev/null",
    ]
    if fastq_path.name.endswith(".gz"):
        cmd.append("--gzip-compressed")
    cmd.append(str(fastq_path))
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    if not report_path.exists():
        raise FileNotFoundError(f"Expected Kraken2 report not found: {report_path}")
    return report_path


def parse_kraken2_report(path: Path) -> List[dict]:
    """Parse a Kraken2 --report file (percent, clade_count, direct_count, rank, taxid, name)."""
    rows: List[dict] = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            pct, clade_count, direct_count, rank, taxid, name = parts[:6]
            try:
                rows.append({
                    "pct": float(pct), "clade_count": int(clade_count),
                    "direct_count": int(direct_count), "rank": rank,
                    "taxid": taxid, "name": name.strip(),
                })
            except ValueError:
                continue
    return rows


def qc_kraken2(
    entity_id: str,
    entity_name: str,
    fastq_path: Path,
    db_path: str,
    outdir: Path,
    expected_species: Optional[str] = None,
    threads: int = 4,
) -> List[Finding]:
    findings: List[Finding] = []

    def add(check: str, status: str, detail: str) -> None:
        findings.append(make_finding(entity_id, entity_name, check, status, detail))

    if not kraken2_available():
        add("kraken2", "FAIL", KRAKEN2_INSTALL_HINT)
        return findings

    try:
        report_path = run_kraken2(fastq_path, db_path, outdir, threads=threads)
    except Exception as exc:
        add("kraken2", "FAIL", f"Kraken2 classification failed: {exc}")
        return findings

    # No PASS/WARN/FAIL judgment here -- Kraken2 doesn't itself assert what
    # counts as "too much" unclassified or contamination, and we don't have a
    # validated threshold for this (RNA-seq reads classify against a
    # genome-based database at rates that don't resemble DNA-seq). Report the
    # numbers; let whoever reads the report judge them.
    rows = parse_kraken2_report(report_path)
    unclassified = next((r for r in rows if r["rank"] == "U"), None)
    unclassified_pct = unclassified["pct"] if unclassified else 0.0
    add("kraken2_unclassified", "INFO",
        f"{unclassified_pct:.1f}% of reads unclassified against this database")

    species_rows = sorted((r for r in rows if r["rank"] == "S"), key=lambda r: -r["pct"])
    top_desc = ", ".join(f"{r['name']} ({r['pct']:.1f}%)" for r in species_rows[:_TOP_N_SPECIES])
    add("kraken2_top_species", "INFO", f"Top classified species: {top_desc or 'none'}")

    if expected_species:
        match = next((r for r in species_rows if r["name"].lower() == expected_species.lower()), None)
        pct = match["pct"] if match else 0.0
        add("kraken2_species_match", "INFO",
            f"{pct:.1f}% of reads classified as expected species {expected_species!r}. Top hits: {top_desc}")

    return findings
