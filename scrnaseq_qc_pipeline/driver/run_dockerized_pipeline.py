#!/usr/bin/env python3
"""Local end-to-end driver for the dockerized scRNA-seq QC checks.

Chains the checks/*/{run.py, Dockerfile} containers together against one
shared, host-mounted work directory, mirroring run_scrnaseq_qc.py's own
control-flow branching (SRA vs. direct-upload, --geo-accession given or not,
loop over mate files, aggregate after) -- so its scrnaseq_qc_report.csv can
be diffed directly against a real run_scrnaseq_qc.py run on the same sample,
as the validation that dockerizing the checks didn't change their behavior.

This is a local testing tool, not a production orchestrator: Synapse-tree
walking/classification (who are the raw-run/processed-sample entities) stays
un-dockerized glue, reused directly from run_scrnaseq_qc.py.

Usage
-----
    # Build every checks/*/Dockerfile, then run:
    python driver/run_dockerized_pipeline.py --synapse-id syn73675040 \\
        --geo-accession GSE293390 --output-dir output/docker_qc --build

    # Images already built, skip fastqc's kraken2 tier:
    python driver/run_dockerized_pipeline.py --synapse-id syn73675040 \\
        --output-dir output/docker_qc

Requires Docker running locally. Images are tagged `qc-<check-name>`
(underscores -> hyphens), e.g. checks/fastq_download -> qc-fastq-download.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

_HERE = Path(__file__).resolve().parent
_PIPELINE_DIR = _HERE.parent
sys.path.insert(0, str(_PIPELINE_DIR))

from qc_common import Finding, summarize, write_report_csv  # noqa: E402
import run_scrnaseq_qc as native  # noqa: E402  -- reuses un-dockerized Synapse-tree-walking glue

CHECKS_DIR = _PIPELINE_DIR / "checks"
IMAGE_PREFIX = "qc-"


def image_tag(check_name: str) -> str:
    return f"{IMAGE_PREFIX}{check_name.replace('_', '-')}"


def build_all_images() -> None:
    for check_dir in sorted(CHECKS_DIR.iterdir()):
        dockerfile = check_dir / "Dockerfile"
        if not dockerfile.exists():
            continue
        tag = image_tag(check_dir.name)
        print(f"Building {tag} …")
        subprocess.run(
            ["docker", "build", "-f", str(dockerfile), "-t", tag, str(_PIPELINE_DIR)],
            check=True,
        )


def docker_run(check_name: str, args: List[str], work_dir: Path,
                extra_mounts: Optional[Dict[str, str]] = None, env: Optional[Dict[str, str]] = None) -> None:
    cmd = ["docker", "run", "--rm", "-v", f"{work_dir}:/data"]
    for host_path, container_path in (extra_mounts or {}).items():
        cmd += ["-v", f"{host_path}:{container_path}"]
    for key, val in (env or {}).items():
        cmd += ["-e", f"{key}={val}"]
    cmd.append(image_tag(check_name))
    cmd += args
    subprocess.run(cmd, check=True)


def read_findings(work_dir: Path, rel_path: str) -> List[Finding]:
    path = work_dir / rel_path
    if not path.exists():
        return []
    with open(path) as fh:
        return json.load(fh)


def read_json(path: Path) -> dict:
    with open(path) as fh:
        return json.load(fh)


# --- Per-tier orchestration (mirrors run_scrnaseq_qc.py's main(), one docker
# run per step instead of one in-process function call) -----------------


def run_geo_tier(work_dir: Path, geo_accession: str, raw_runs: Dict[str, dict]) -> "tuple[List[Finding], dict]":
    findings: List[Finding] = []
    docker_run("geo_metadata_fetch", [
        "--geo-accession", geo_accession,
        "--output", "/data/findings/geo_metadata_fetch.json",
        "--geo-metadata-out", "/data/geo_metadata.json",
    ], work_dir)
    findings.extend(read_findings(work_dir, "findings/geo_metadata_fetch.json"))
    geo_metadata_path = work_dir / "geo_metadata.json"
    geo_metadata = read_json(geo_metadata_path) if geo_metadata_path.exists() else {"samples": {}}

    for entity_id, info in raw_runs.items():
        ann_rel = f"annotations/{entity_id}.json"
        ann_path = work_dir / ann_rel
        ann_path.parent.mkdir(parents=True, exist_ok=True)
        with open(ann_path, "w") as fh:
            json.dump(info["annotations"], fh, default=str)
        out_rel = f"findings/geo_annotation_check_{entity_id}.json"
        docker_run("geo_annotation_check", [
            "--entity-id", entity_id, "--sample-name", info["name"],
            "--annotations-json", f"/data/{ann_rel}",
            "--geo-metadata-json", "/data/geo_metadata.json",
            "--output", f"/data/{out_rel}",
        ], work_dir)
        findings.extend(read_findings(work_dir, out_rel))

    return findings, geo_metadata


def run_fastq_sra(work_dir: Path, entity_id: str, info: dict, expected_species: Optional[str],
                   kraken2_db: Optional[str], max_download_mb: Optional[float], keep_files: bool) -> List[Finding]:
    findings: List[Finding] = []
    sample_dir = f"raw/{info['name']}"
    ena_meta_rel = f"{sample_dir}/ena_metadata.json"

    docker_run("fastq_resolution", [
        "--entity-id", entity_id, "--sample-name", info["name"], "--srr-accession", info["srr"],
        "--output", f"/data/findings/fastq_resolution_{entity_id}.json",
        "--ena-metadata-out", f"/data/{ena_meta_rel}",
    ], work_dir)
    findings.extend(read_findings(work_dir, f"findings/fastq_resolution_{entity_id}.json"))
    ena_meta_path = work_dir / ena_meta_rel
    if not ena_meta_path.exists():
        return findings
    links = read_json(ena_meta_path).get("links", [])

    max_bytes = int(max_download_mb * 1024 * 1024) if max_download_mb else None
    any_capped = False
    summary_rels: List[str] = []

    for item in links:
        url = item["link"]
        fname = url.rsplit("/", 1)[-1]
        dest_rel = f"{sample_dir}/{fname}"

        dl_args = [
            "--entity-id", entity_id, "--sample-name", info["name"],
            "--url", url, "--dest", f"/data/{dest_rel}",
            "--registered-size", str(item.get("file_size", "")),
            "--registered-md5", item.get("md5", ""),
            "--output", f"/data/findings/fastq_download_{entity_id}_{fname}.json",
        ]
        if max_download_mb:
            dl_args += ["--max-download-mb", str(max_download_mb)]
        if keep_files:
            dl_args.append("--keep-files")
        docker_run("fastq_download", dl_args, work_dir)
        findings.extend(read_findings(work_dir, f"findings/fastq_download_{entity_id}_{fname}.json"))

        fastq_path = work_dir / dest_rel
        if not fastq_path.exists():
            continue
        was_capped = bool(max_bytes) and fastq_path.stat().st_size >= max_bytes
        any_capped = any_capped or was_capped

        docker_run("fq_lint", [
            "--entity-id", entity_id, "--sample-name", info["name"],
            "--fastq-path", f"/data/{dest_rel}",
            "--output", f"/data/findings/fq_lint_{entity_id}_{fname}.json",
        ], work_dir)
        findings.extend(read_findings(work_dir, f"findings/fq_lint_{entity_id}_{fname}.json"))

        if kraken2_db:
            kraken2_args = [
                "--entity-id", entity_id, "--sample-name", info["name"],
                "--fastq-path", f"/data/{dest_rel}", "--db-path", "/db",
                "--outdir", f"/data/{sample_dir}/kraken2",
                "--output", f"/data/findings/kraken2_{entity_id}_{fname}.json",
            ]
            if expected_species:
                kraken2_args += ["--expected-species", expected_species]
            docker_run("kraken2_qc", kraken2_args, work_dir, extra_mounts={kraken2_db: "/db"})
            findings.extend(read_findings(work_dir, f"findings/kraken2_{entity_id}_{fname}.json"))

        summary_rel = f"{sample_dir}/fastqc/{fname}.summary.json"
        docker_run("fastqc", [
            "--entity-id", entity_id, "--sample-name", info["name"],
            "--fastq-path", f"/data/{dest_rel}", "--outdir", f"/data/{sample_dir}/fastqc",
            "--label", fname,
            "--output", f"/data/findings/fastqc_{entity_id}_{fname}.json",
            "--summary-out", f"/data/{summary_rel}",
        ], work_dir)
        findings.extend(read_findings(work_dir, f"findings/fastqc_{entity_id}_{fname}.json"))
        if (work_dir / summary_rel).exists():
            summary_rels.append(summary_rel)

        if not keep_files:
            fastq_path.unlink(missing_ok=True)

    if summary_rels:
        agg_args = [
            "--entity-id", entity_id, "--sample-name", info["name"],
            "--ena-metadata-json", f"/data/{ena_meta_rel}",
            "--output", f"/data/findings/fastq_aggregate_{entity_id}.json",
        ]
        for rel in summary_rels:
            agg_args += ["--fastqc-summary", f"/data/{rel}"]
        if any_capped:
            agg_args.append("--any-capped")
        docker_run("fastq_aggregate", agg_args, work_dir)
        findings.extend(read_findings(work_dir, f"findings/fastq_aggregate_{entity_id}.json"))

    return findings


def run_fastq_direct(work_dir: Path, entity_id: str, info: dict, expected_species: Optional[str],
                      kraken2_db: Optional[str], keep_files: bool) -> List[Finding]:
    findings: List[Finding] = []
    sample_dir = f"raw/{info['name']}"

    docker_run("synapse_fastq_download", [
        "--entity-id", entity_id, "--sample-name", info["name"],
        "--work-dir", "/data/raw",
        "--output", f"/data/findings/fastq_download_{entity_id}.json",
    ], work_dir, env={"SYNAPSE_AUTH_TOKEN": os.environ.get("SYNAPSE_AUTH_TOKEN", "")})
    findings.extend(read_findings(work_dir, f"findings/fastq_download_{entity_id}.json"))

    sample_path = work_dir / sample_dir
    if not sample_path.exists():
        return findings
    candidates = [p for p in sample_path.iterdir()
                  if p.is_file() and p.name.lower().endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))]
    if not candidates:
        return findings
    fastq_path = candidates[0]
    dest_rel = f"{sample_dir}/{fastq_path.name}"

    docker_run("fq_lint", [
        "--entity-id", entity_id, "--sample-name", info["name"],
        "--fastq-path", f"/data/{dest_rel}",
        "--output", f"/data/findings/fq_lint_{entity_id}.json",
    ], work_dir)
    findings.extend(read_findings(work_dir, f"findings/fq_lint_{entity_id}.json"))

    if kraken2_db:
        kraken2_args = [
            "--entity-id", entity_id, "--sample-name", info["name"],
            "--fastq-path", f"/data/{dest_rel}", "--db-path", "/db",
            "--outdir", f"/data/{sample_dir}/kraken2",
            "--output", f"/data/findings/kraken2_{entity_id}.json",
        ]
        if expected_species:
            kraken2_args += ["--expected-species", expected_species]
        docker_run("kraken2_qc", kraken2_args, work_dir, extra_mounts={kraken2_db: "/db"})
        findings.extend(read_findings(work_dir, f"findings/kraken2_{entity_id}.json"))

    docker_run("fastqc", [
        "--entity-id", entity_id, "--sample-name", info["name"],
        "--fastq-path", f"/data/{dest_rel}", "--outdir", f"/data/{sample_dir}/fastqc",
        "--output", f"/data/findings/fastqc_{entity_id}.json",
    ], work_dir)
    findings.extend(read_findings(work_dir, f"findings/fastqc_{entity_id}.json"))

    if not keep_files:
        fastq_path.unlink(missing_ok=True)

    return findings


def run_matrix(work_dir: Path, sample_id: str, files: Dict[str, str], cellranger_output_class: str) -> List[Finding]:
    findings: List[Finding] = []
    sample_dir = f"processed/{sample_id}"

    docker_run("matrix_download", [
        "--sample-id", sample_id,
        "--barcodes-entity-id", files.get("barcodes", ""),
        "--features-entity-id", files.get("features", ""),
        "--matrix-entity-id", files.get("matrix", ""),
        "--work-dir", "/data/processed",
        "--output", f"/data/findings/matrix_download_{sample_id}.json",
    ], work_dir, env={"SYNAPSE_AUTH_TOKEN": os.environ.get("SYNAPSE_AUTH_TOKEN", "")})
    findings.extend(read_findings(work_dir, f"findings/matrix_download_{sample_id}.json"))

    sample_path = work_dir / sample_dir
    if not sample_path.exists():
        return findings
    barcodes = next((p for p in sample_path.iterdir() if "barcodes" in p.name.lower()), None)
    features = next((p for p in sample_path.iterdir() if "features" in p.name.lower() or "genes" in p.name.lower()), None)
    matrix = next((p for p in sample_path.iterdir() if "matrix" in p.name.lower()), None)
    if not (barcodes and features and matrix):
        return findings

    docker_run("matrix_content_qc", [
        "--sample-id", sample_id,
        "--barcodes-path", f"/data/{sample_dir}/{barcodes.name}",
        "--features-path", f"/data/{sample_dir}/{features.name}",
        "--matrix-path", f"/data/{sample_dir}/{matrix.name}",
        "--cellranger-output-class", cellranger_output_class or "",
        "--output", f"/data/findings/matrix_content_{sample_id}.json",
    ], work_dir)
    findings.extend(read_findings(work_dir, f"findings/matrix_content_{sample_id}.json"))
    return findings


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--synapse-id", required=True)
    p.add_argument("--geo-accession", default=None)
    p.add_argument("--output-dir", "-d", default="output/docker_qc")
    p.add_argument("--work-dir", default=None)
    p.add_argument("--kraken2-db", default=None, help="Host path to an extracted Kraken2 database directory")
    p.add_argument("--max-download-mb", type=float, default=None)
    p.add_argument("--keep-files", action="store_true")
    p.add_argument("--skip-fastqc", action="store_true")
    p.add_argument("--skip-matrix", action="store_true")
    p.add_argument("--build", action="store_true", help="Build every checks/*/Dockerfile before running")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    if args.build:
        build_all_images()

    out_dir = Path(args.output_dir)
    work_dir = Path(args.work_dir) if args.work_dir else out_dir / "docker_work"
    (work_dir / "findings").mkdir(parents=True, exist_ok=True)

    print("Logging into Synapse …")
    syn = native._synapse_login()

    print(f"Scanning Synapse tree under {args.synapse_id} …")
    children = native.walk_file_entities(syn, args.synapse_id)
    print(f"  {len(children)} file(s) found.")

    raw_runs: Dict[str, dict] = {}
    processed_groups: Dict[str, Dict[str, str]] = {}
    processed_meta: Dict[str, dict] = {}

    for child in children:
        entity_id = child["id"]
        ent = syn.get(entity_id, downloadFile=False)
        annotations = dict(ent.annotations)
        role = native.classify_entity(ent)
        if role == "raw_sra":
            srr = native._SRA_RUN_URL_RE.search(ent._file_handle.get("externalURL", "")).group(1).upper()
            raw_runs[entity_id] = dict(name=ent.name, annotations=annotations, srr=srr)
        elif role == "raw_fastq_direct":
            raw_runs[entity_id] = dict(name=ent.name, annotations=annotations, srr=None)
        elif role.startswith("processed_"):
            sample_id = native.first_value(annotations.get("SampleId")) or ent.name
            sub_role = role.split("_", 1)[1]
            processed_groups.setdefault(sample_id, {})[sub_role] = entity_id
            processed_meta[sample_id] = annotations

    print(f"  {len(raw_runs)} raw run(s), {len(processed_groups)} processed sample group(s).")

    all_findings: List[Finding] = []
    expected_species_by_name: Dict[str, str] = {}

    if args.geo_accession:
        geo_findings, geo_metadata = run_geo_tier(work_dir, args.geo_accession, raw_runs)
        all_findings.extend(geo_findings)
        for gsm, row in geo_metadata.get("samples", {}).items():
            species = row.get("sample_organism_ch1")
            if species:
                expected_species_by_name[gsm] = str(species)
    else:
        print("No --geo-accession given; skipping Synapse-annotations-vs-GEO check.")

    if args.skip_fastqc:
        print("Skipping raw FastQC tier (--skip-fastqc).")
    else:
        for entity_id, info in raw_runs.items():
            print(f"Running dockerized fastq QC for {info['name']} …")
            expected_species = expected_species_by_name.get(info["name"])
            if info["srr"]:
                all_findings.extend(run_fastq_sra(
                    work_dir, entity_id, info, expected_species, args.kraken2_db,
                    args.max_download_mb, args.keep_files,
                ))
            else:
                all_findings.extend(run_fastq_direct(
                    work_dir, entity_id, info, expected_species, args.kraken2_db, args.keep_files,
                ))

    if args.skip_matrix:
        print("Skipping processed matrix QC tier (--skip-matrix).")
    else:
        for sample_id, files in processed_groups.items():
            print(f"Running dockerized matrix QC for sample {sample_id} …")
            all_findings.extend(run_matrix(
                work_dir, sample_id, files,
                native.first_value(processed_meta[sample_id].get("CellrangerOutputClass")),
            ))

    report_path = out_dir / "scrnaseq_qc_report.csv"
    write_report_csv(all_findings, report_path)
    print(f"\nQC complete: {summarize(all_findings)}")
    print(f"Report written to {report_path}")
    print("Diff this against a real run_scrnaseq_qc.py run on the same --synapse-id "
          "to validate that dockerizing the checks didn't change their behavior.")

    return 1 if any(f["status"] == "FAIL" for f in all_findings) else 0


if __name__ == "__main__":
    sys.exit(main())
