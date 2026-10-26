#!/usr/bin/env python3
"""Local end-to-end driver for the dockerized scRNA-seq QC checks.

Chains the checks/{external,internal}/*/{run.py, Dockerfile} containers
together against one shared, host-mounted work directory, mirroring
run_scrnaseq_qc.py's own control-flow branching (SRA vs. direct-upload,
--geo-accession given or not, loop over mate files, aggregate after) -- so
its scrnaseq_qc_report.csv can be diffed against a real run_scrnaseq_qc.py
run on the same sample, as a validation that dockerizing the checks didn't
change their behavior.

Two invocation styles, matching the two check categories:
  - checks/external/*: invoked exactly the way nf-dcqc's RUN_PROCESS would
    (see Sage-Bionetworks-Workflows/py-dcqc#96) -- a fresh task directory
    holding only the symlinked input file(s), run as the current host user,
    judged only by exit_code.txt/std_out.txt/std_err.txt. This driver's own
    report is rebuilt from that raw 3-way result, since py-dcqc's one-test-
    one-verdict model doesn't preserve this pipeline's original per-module
    Finding granularity 1:1 -- see the adapter functions below.
  - checks/internal/*: invoked with this repo's own --output-a-Findings-JSON
    contract, unchanged from before this split (these were never going to
    be real py-dcqc external tests -- see the dockerize-checks skill).

This is a local testing tool, not a production orchestrator: Synapse-tree
walking/classification (who are the raw-run/processed-sample entities) stays
un-dockerized glue, reused directly from run_scrnaseq_qc.py.

Usage
-----
    # Build every checks/{external,internal}/*/Dockerfile, then run:
    python driver/run_dockerized_pipeline.py --synapse-id syn73675040 \\
        --geo-accession GSE293390 --output-dir output/docker_qc --build

    # Images already built:
    python driver/run_dockerized_pipeline.py --synapse-id syn73675040 \\
        --output-dir output/docker_qc

Requires Docker running locally. Images are tagged `qc-<check-name>`
(underscores -> hyphens), e.g. checks/external/fastq_download -> qc-fastq-download.
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

from qc_common import Finding, make_finding, summarize, write_report_csv  # noqa: E402
import run_scrnaseq_qc as native  # noqa: E402  -- reuses un-dockerized Synapse-tree-walking glue

CHECKS_DIR = _PIPELINE_DIR / "checks"
IMAGE_PREFIX = "qc-"


def image_tag(check_name: str) -> str:
    return f"{IMAGE_PREFIX}{check_name.replace('_', '-')}"


def build_all_images() -> None:
    for category in ("external", "internal"):
        for check_dir in sorted((CHECKS_DIR / category).iterdir()):
            dockerfile = check_dir / "Dockerfile"
            if not dockerfile.exists():
                continue
            tag = image_tag(check_dir.name)
            print(f"Building {tag} ({category}) …")
            subprocess.run(
                ["docker", "build", "-f", str(dockerfile), "-t", tag, str(_PIPELINE_DIR)],
                check=True,
            )


# --- checks/internal/* invocation (unchanged --output-JSON contract) ------


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


# --- checks/external/* invocation (real nf-dcqc RUN_PROCESS semantics) ----


def docker_run_external(
    check_name: str, task_dir: Path, input_files: Dict[str, Path],
    extra_args: Optional[List[str]] = None, extra_mounts: Optional[Dict[str, str]] = None,
) -> "tuple[int, dict, str]":
    """Invoke an external/ check exactly the way nf-dcqc's RUN_PROCESS does
    (modules/local/run_process.nf in Sage-Bionetworks-Workflows/nf-dcqc): a
    fresh task directory holding only the staged input file(s) -- symlinked
    in by their bare filename, matching how nf-dcqc stages a file -- run as
    the current host user with that directory as the container's working
    directory, judged only by the exit code and the two captured streams.

    Returns (exit_code, metrics_dict_or_empty, reason_text). metrics is {}
    if stdout wasn't valid JSON (py-dcqc would call that an error regardless
    of the exit code; this driver just reports an empty metrics dict and
    lets the caller's own exit-code handling decide the status).
    """
    task_dir.mkdir(parents=True, exist_ok=True)
    arg_names = []
    for fname, src in input_files.items():
        link = task_dir / fname
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(src.resolve())
        arg_names.append(fname)

    quoted_names = " ".join(f"'{n}'" for n in arg_names)
    inner_command = " ".join(["python3", "/app/run.py", *(extra_args or []), quoted_names]).strip()
    task_abs = str(task_dir.resolve())

    cmd = [
        "docker", "run", "--rm",
        "-u", f"{os.getuid()}:{os.getgid()}",
        "-v", f"{task_abs}:{task_abs}", "-w", task_abs,
    ]
    for host_path, container_path in (extra_mounts or {}).items():
        cmd += ["-v", f"{host_path}:{container_path}"]
    cmd.append(image_tag(check_name))
    cmd += [
        "/bin/bash", "-euo", "pipefail", "-c",
        f'( ({inner_command}) > std_out.txt 2> std_err.txt; echo $? > exit_code.txt ) || true',
    ]
    subprocess.run(cmd, check=True)

    exit_code = int((task_dir / "exit_code.txt").read_text().strip())
    stdout_text = (task_dir / "std_out.txt").read_text().strip()
    stderr_text = (task_dir / "std_err.txt").read_text().strip()
    try:
        metrics = json.loads(stdout_text) if stdout_text else {}
    except json.JSONDecodeError:
        metrics = {}
    return exit_code, metrics, stderr_text


def run_fq_lint_external(work_dir: Path, entity_id: str, sample_name: str, fastq_path: Path) -> List[Finding]:
    task_dir = work_dir / "_tasks" / f"fq_lint_{entity_id}_{fastq_path.name}"
    exit_code, _metrics, reason = docker_run_external("fq_lint", task_dir, {fastq_path.name: fastq_path})

    if exit_code == 0:
        return [make_finding(entity_id, sample_name, "fq_lint", "PASS",
                              f"{fastq_path.name}: no structural issues found")]
    if exit_code == 1:
        return [make_finding(entity_id, sample_name, "fq_lint", "FAIL",
                              reason or f"{fastq_path.name}: fq lint reported issues")]
    # exit_code == 2 (or anything else): an operational failure, not fq's verdict.
    return [make_finding(entity_id, sample_name, "fq_lint", "FAIL",
                          reason or f"{fastq_path.name}: fq lint could not run")]


def run_fastqc_external(
    work_dir: Path, entity_id: str, sample_name: str, fastq_path: Path,
) -> "tuple[List[Finding], Optional[dict]]":
    task_dir = work_dir / "_tasks" / f"fastqc_{entity_id}_{fastq_path.name}"
    exit_code, metrics, reason = docker_run_external(
        "fastqc", task_dir, {fastq_path.name: fastq_path}, extra_args=["--summary-out", "fastqc_summary.json"],
    )

    if exit_code == 2:
        return [make_finding(entity_id, sample_name, "fastqc", "FAIL", reason or "FastQC failed")], None

    findings: List[Finding] = [make_finding(
        entity_id, sample_name, "fastqc_basic_stats", "INFO",
        f"{fastq_path.name}: {metrics.get('total_sequences', '?')} seqs, "
        f"len={metrics.get('sequence_length', '?')}, GC={metrics.get('pct_gc', '?')}%",
    )]
    module_details = metrics.get("module_details") or {}
    for module_name, module_status in (metrics.get("modules") or {}).items():
        if module_status == "PASS":
            continue
        detail = module_details.get(module_name, f"{fastq_path.name}: {module_status.lower()}")
        findings.append(make_finding(entity_id, sample_name, f"fastqc_module:{module_name}", module_status, detail))

    summary_path = task_dir / "fastqc_summary.json"
    summary = read_json(summary_path) if summary_path.exists() else None
    return findings, summary


def run_kraken2_external(
    work_dir: Path, entity_id: str, sample_name: str, fastq_path: Path,
    kraken2_db: str, expected_species: Optional[str], threads: int,
) -> List[Finding]:
    task_dir = work_dir / "_tasks" / f"kraken2_{entity_id}_{fastq_path.name}"
    extra_args = ["--db-path", "/db", "--threads", str(threads)]
    if expected_species:
        extra_args += ["--expected-species", expected_species]
    exit_code, metrics, reason = docker_run_external(
        "kraken2_qc", task_dir, {fastq_path.name: fastq_path},
        extra_args=extra_args, extra_mounts={kraken2_db: "/db"},
    )

    if exit_code != 0:
        # Always an operational failure for this check -- it has no real
        # FAIL verdict of its own (see checks/external/kraken2_qc/run.py).
        return [make_finding(entity_id, sample_name, "kraken2", "FAIL", reason or "Kraken2 classification failed")]

    findings = [
        make_finding(entity_id, sample_name, "kraken2_unclassified", "INFO", metrics.get("unclassified", "")),
        make_finding(entity_id, sample_name, "kraken2_top_species", "INFO", metrics.get("top_species", "")),
    ]
    if "species_match" in metrics:
        findings.append(make_finding(entity_id, sample_name, "kraken2_species_match", "INFO", metrics["species_match"]))
    return findings


def run_matrix_content_external(
    work_dir: Path, sample_id: str, barcodes: Path, features: Path, matrix: Path, cellranger_output_class: str,
) -> List[Finding]:
    task_dir = work_dir / "_tasks" / f"matrix_content_{sample_id}"
    input_files = {barcodes.name: barcodes, features.name: features, matrix.name: matrix}
    exit_code, metrics, reason = docker_run_external(
        "matrix_content_qc", task_dir, input_files,
        extra_args=["--cellranger-output-class", cellranger_output_class or ""],
    )

    if exit_code == 2:
        return [make_finding(sample_id, sample_id, "matrix_integrity", "FAIL", reason or "Could not parse matrix")]

    findings: List[Finding] = []
    for check, key in (
        ("matrix_integrity", "matrix_integrity"), ("matrix_dimensions", "matrix_dimensions"),
        ("matrix_metrics", "matrix_metrics"), ("empty_barcode_check", "empty_barcode_check"),
    ):
        if metrics.get(key):
            status = "FAIL" if (key == "matrix_metrics" and exit_code == 1) else "INFO"
            findings.append(make_finding(sample_id, sample_id, check, status, metrics[key]))
    if metrics.get("mito_fraction"):
        findings.append(make_finding(sample_id, sample_id, "mito_fraction", "INFO", metrics["mito_fraction"]))
    return findings


# --- Per-tier orchestration (mirrors run_scrnaseq_qc.py's main()) ---------


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

        findings.extend(run_fq_lint_external(work_dir, entity_id, info["name"], fastq_path))

        if kraken2_db:
            findings.extend(run_kraken2_external(
                work_dir, entity_id, info["name"], fastq_path, kraken2_db, expected_species, threads=4,
            ))

        fastqc_findings, summary = run_fastqc_external(work_dir, entity_id, info["name"], fastq_path)
        findings.extend(fastqc_findings)
        if summary is not None:
            summary_rel = f"{sample_dir}/fastqc/{fname}.summary.json"
            summary_path = work_dir / summary_rel
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            with open(summary_path, "w") as fh:
                json.dump(summary, fh, indent=2)
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

    findings.extend(run_fq_lint_external(work_dir, entity_id, info["name"], fastq_path))

    if kraken2_db:
        findings.extend(run_kraken2_external(
            work_dir, entity_id, info["name"], fastq_path, kraken2_db, expected_species, threads=4,
        ))

    fastqc_findings, _summary = run_fastqc_external(work_dir, entity_id, info["name"], fastq_path)
    findings.extend(fastqc_findings)

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

    findings.extend(run_matrix_content_external(work_dir, sample_id, barcodes, features, matrix, cellranger_output_class))
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
    p.add_argument("--build", action="store_true", help="Build every checks/{external,internal}/*/Dockerfile before running")
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
