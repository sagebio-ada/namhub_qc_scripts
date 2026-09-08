"""Species/contamination screening via Kraken2 taxonomic classification.

Nothing else in this pipeline ever looks at whether the actual sequence
CONTENT matches the declared organism. GEO's own `sample_organism_ch1`
field is just as self-reported as Synapse's `species`/`referenceSet`
annotations, so cross-checking one against the other (geo_metadata_qc.py)
can't catch a mislabeled sample or cross-species contamination -- it would
just compare one unverified label to another. Kraken2 classifies reads
against a real reference database, giving an independent, content-based
read on what's actually in the file.

Kraken2 needs a pre-built database staged locally (not something this
pipeline downloads for you -- they run several GB). See
https://benlangmead.github.io/aws-indexes/k2 -- "Standard-8" (~5.5GB
compressed; RefSeq archaea/bacteria/viral/plasmid/human/UniVec) is a
reasonable general-purpose default. Pass its extracted directory via
--kraken2-db.

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


def kraken2_available() -> bool:
    return shutil.which("kraken2") is not None


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
