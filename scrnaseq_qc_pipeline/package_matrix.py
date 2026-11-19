#!/usr/bin/env python3
"""Package barcodes/features/matrix files into one .tar.gz per sample, so
the matrix content check (py-dcqc file type CELLRANGER-MATRIX-TAR) can run
on each sample as one file.

This is a step done BEFORE py-dcqc, not part of any py-dcqc test: the
archives it makes are the files py-dcqc is given.

Three modes:

  From Synapse file entities: find every set under a Synapse project or
  folder, download it, write one archive per set into OUTDIR, and (with
  --upload-to) store each archive in Synapse as a new file entity.
    python3 package_matrix.py --synapse-id syn123 -o archives/ --upload-to syn456

  From a local folder: the same, for files already on disk.
    python3 package_matrix.py --scan downloads/ -o archives/

  One sample, files given explicitly:
    python3 package_matrix.py GSM123_barcodes.tsv.gz GSM123_features.tsv.gz \\
        GSM123_matrix.mtx.gz -o GSM123_matrix.tar.gz

How sets are found: a file belongs to a set by its folder (Synapse parent)
and the part of its name before the standard CellRanger name, so all of
these work:
  GSM123_barcodes.tsv.gz, GSM123_features.tsv.gz, GSM123_matrix.mtx.gz
      -> GSM123_matrix.tar.gz
  sample1/filtered_feature_bc_matrix/{barcodes.tsv.gz,features.tsv.gz,matrix.mtx.gz}
      -> sample1_filtered_feature_bc_matrix.tar.gz (Synapse: named after the
         three files' SampleId annotation if they share one, else the folder)
  GSM5_raw_barcodes.tsv.gz ... and GSM5_filtered_barcodes.tsv.gz ...
      -> GSM5_raw_matrix.tar.gz and GSM5_filtered_matrix.tar.gz
CellRanger v2's genes.tsv counts as the features file. A set missing a file,
or with two candidates for one, is reported on stderr and skipped. Archives
already present (e.g. GEO's <sample>_filtered_feature_bc_matrix.tar.gz) are
left alone: the check reads them as they are.

Inside each archive the files have their standard names (barcodes.tsv[.gz],
features.tsv[.gz], matrix.mtx[.gz]), whatever they're called on disk. Each
file is stored as-is; a .gz file stays gzipped inside the archive.

Synapse upload (--upload-to FOLDER): each archive becomes a file entity in
FOLDER, with the annotations the three source files share (e.g. SampleId)
and provenance "used" pointing at the three source entities. Storing an
archive whose name is already in FOLDER makes a new version of that entity.
Needs SYNAPSE_AUTH_TOKEN (or a cached Synapse login).

Output: one line per archive written (with its Synapse ID when uploaded),
on stdout. Exit codes: 0 if every set found was complete; 1 if any was
skipped (the complete ones are still written); 2 for bad arguments or a
missing input file.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Callable, Dict, Hashable, Iterable, List, Optional, Tuple

# Standard names for each role, without .gz. genes.tsv is the CellRanger v2
# name for features.tsv. (Same names matrix_qc.extract_matrix_archive
# looks for inside an archive.)
ROLE_NAMES = {
    "barcodes": ("barcodes.tsv",),
    "features": ("features.tsv", "genes.tsv"),
    "matrix": ("matrix.mtx",),
}
ROLES = tuple(ROLE_NAMES)


def package_matrix(barcodes: Path, features: Path, matrix: Path, out_path: Path) -> Path:
    """Write barcodes/features/matrix into out_path (.tar.gz) under their
    standard names. Returns out_path."""
    members = (("barcodes.tsv", barcodes), ("features.tsv", features), ("matrix.mtx", matrix))
    for _, path in members:
        if not path.is_file():
            raise FileNotFoundError(f"{path}: file not found")
    with tarfile.open(out_path, "w:gz") as tar:
        for standard_name, path in members:
            suffix = ".gz" if path.name.endswith(".gz") else ""
            tar.add(path, arcname=standard_name + suffix, recursive=False)
    return out_path


def file_role(name: str) -> Optional[Tuple[str, str]]:
    """("GSM123_barcodes.tsv.gz") -> ("barcodes", "GSM123"): the role and
    the name prefix before the standard name (separators stripped), or None
    if the name isn't one of the three files."""
    base = name[: -len(".gz")] if name.endswith(".gz") else name
    for role, names in ROLE_NAMES.items():
        for standard in names:
            if base == standard or (base.endswith(standard) and base[-len(standard) - 1] in "_.-"):
                return role, base[: -len(standard)].rstrip("_.-")
    return None


def archive_name(stem: str) -> str:
    """The archive file name for a set called `stem`."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("_") or "matrix"
    return f"{stem}.tar.gz" if stem.lower().endswith("matrix") else f"{stem}_matrix.tar.gz"


def group_matrix_files(
    items: Iterable[Tuple[Hashable, str, object]],
    name_for: Callable[[Hashable, str, Dict[str, object]], str],
    label_for: Callable[[Hashable, str], str],
) -> "Tuple[Dict[str, Dict[str, object]], List[str]]":
    """Group (folder key, file name, item) triples into barcodes/features/
    matrix sets, by folder and name prefix.

    name_for(folder, prefix, {role: item}) names a complete set's archive;
    label_for(folder, prefix) names a set in problem messages. Returns
    ({archive name: {role: item}} for complete sets, [problem messages]).
    """
    groups: Dict[Tuple[Hashable, str], Dict[str, List[object]]] = defaultdict(lambda: defaultdict(list))
    for folder, name, item in items:
        match = file_role(name)
        if match:
            role, prefix = match
            groups[(folder, prefix)][role].append((name, item))

    sets: Dict[str, Dict[str, object]] = {}
    problems: List[str] = []
    for (folder, prefix), by_role in sorted(groups.items(), key=lambda kv: (str(kv[0][0]), kv[0][1])):
        label = label_for(folder, prefix)
        missing = [role for role in ROLES if not by_role.get(role)]
        duplicated = [role for role in ROLES if len(by_role.get(role, [])) > 1]
        if missing or duplicated:
            details = []
            if missing:
                details.append(f"missing {', '.join(missing)}")
            for role in duplicated:
                details.append(f"more than one {role} file ({', '.join(n for n, _ in by_role[role])})")
            problems.append(f"{label}: {'; '.join(details)}")
            continue
        chosen = {role: by_role[role][0][1] for role in ROLES}
        name = name_for(folder, prefix, chosen)
        if name in sets:
            problems.append(f"{label}: archive name {name} is already used by another set")
            continue
        sets[name] = chosen
    return sets, problems


# --- Local folder ---------------------------------------------------------


def find_local_sets(root: Path) -> "Tuple[Dict[str, Dict[str, Path]], List[str]]":
    """The barcodes/features/matrix sets among the files under root."""
    items = [(path.parent, path.name, path) for path in sorted(root.rglob("*")) if path.is_file()]

    def name_for(folder: Path, prefix: str, _files) -> str:
        if prefix:
            return archive_name(prefix)
        # No prefix (CellRanger's own layout): name it after the folder path.
        rel = folder.relative_to(root)
        return archive_name("_".join(rel.parts) if rel.parts else root.resolve().name)

    def label_for(folder: Path, prefix: str) -> str:
        return str(folder.relative_to(root) / (prefix or "(no prefix)"))

    return group_matrix_files(items, name_for, label_for)


# --- Synapse --------------------------------------------------------------


def _first(value):
    return value[0] if isinstance(value, (list, tuple)) and value else value


def _walk_synapse_files(syn, root_id: str) -> List[dict]:
    """Every file entity under a Synapse project/folder: {"id", "name",
    "parentId"}."""
    out: List[dict] = []
    stack = [root_id]
    while stack:
        parent = stack.pop()
        for child in syn.getChildren(parent, includeTypes=["folder", "file"]):
            if child["type"].endswith("Folder"):
                stack.append(child["id"])
            elif child["type"].endswith("FileEntity"):
                out.append({"id": child["id"], "name": child["name"], "parentId": parent})
    return out


def find_synapse_sets(syn, root_id: str) -> "Tuple[Dict[str, Dict[str, dict]], List[str]]":
    """The barcodes/features/matrix sets among the file entities under
    root_id. Each item is {"id", "name", "parentId"}."""
    entities = _walk_synapse_files(syn, root_id)
    folder_names: Dict[str, str] = {}

    def folder_name(folder_id: str) -> str:
        if folder_id not in folder_names:
            folder_names[folder_id] = syn.get(folder_id, downloadFile=False).name
        return folder_names[folder_id]

    def name_for(folder_id: str, prefix: str, files: Dict[str, dict]) -> str:
        if prefix:
            return archive_name(prefix)
        # No prefix: use the SampleId the three files share, if any.
        sample_ids = {
            _first(syn.get(f["id"], downloadFile=False).annotations.get("SampleId")) for f in files.values()
        }
        if len(sample_ids) == 1 and None not in sample_ids:
            return archive_name(str(sample_ids.pop()))
        return archive_name(folder_name(folder_id))

    def label_for(folder_id: str, prefix: str) -> str:
        return f"{folder_id}/{prefix or '(no prefix)'}"

    items = [(e["parentId"], e["name"], e) for e in entities]
    return group_matrix_files(items, name_for, label_for)


def shared_annotations(entities: List) -> Dict[str, object]:
    """The annotations every entity has with the same value."""
    common = dict(entities[0].annotations)
    for ent in entities[1:]:
        common = {k: v for k, v in common.items() if k in ent.annotations and ent.annotations[k] == v}
    return common


def package_from_synapse(syn, root_id: str, out_dir: Path, upload_to: Optional[str]) -> int:
    import synapseclient

    sets, problems = find_synapse_sets(syn, root_id)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, files in sets.items():
        with tempfile.TemporaryDirectory() as tmp:
            entities = {}
            for role, f in files.items():
                # Each file in its own folder, so two files with the same name
                # can't overwrite each other.
                entities[role] = syn.get(f["id"], downloadLocation=os.path.join(tmp, role), ifcollision="overwrite.local")
            out_path = package_matrix(
                Path(entities["barcodes"].path), Path(entities["features"].path), Path(entities["matrix"].path),
                out_dir / name,
            )
        line = str(out_path)
        if upload_to:
            stored = syn.store(
                synapseclient.File(
                    str(out_path), parent=upload_to, annotations=shared_annotations(list(entities.values())),
                ),
                used=[ent.id for ent in entities.values()],
                activityName="Package CellRanger matrix",
                activityDescription="barcodes/features/matrix packaged into one tar for py-dcqc "
                                    "(scrnaseq_qc_pipeline/package_matrix.py)",
            )
            line += f"\t{stored.id}"
        print(line)
    for problem in problems:
        print(f"skipped {problem}", file=sys.stderr)
    if not sets and not problems:
        print(f"No barcodes/features/matrix files found under {root_id}", file=sys.stderr)
    return 1 if problems else 0


# --- CLI ------------------------------------------------------------------


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="*", type=Path,
                   help="barcodes, features (or genes) and matrix file, in that order (one-sample mode)")
    source = p.add_mutually_exclusive_group()
    source.add_argument("--synapse-id", metavar="SYN_ID",
                        help="Find every set among the file entities under this Synapse project/folder")
    source.add_argument("--scan", type=Path, metavar="DIR", help="Find every set among the files under DIR")
    p.add_argument("-o", "--output", type=Path, required=True,
                   help="With --synapse-id or --scan: the folder to write archives into. "
                        "Otherwise: the archive to write, e.g. <sample>_matrix.tar.gz")
    p.add_argument("--upload-to", metavar="SYN_ID",
                   help="With --synapse-id: store each archive as a file entity in this Synapse folder")
    args = p.parse_args()

    if args.upload_to and not args.synapse_id:
        print("--upload-to needs --synapse-id", file=sys.stderr)
        return 2
    if (args.synapse_id or args.scan) and args.files:
        print("Give --synapse-id, --scan, or three files, not more than one", file=sys.stderr)
        return 2

    if args.synapse_id:
        import synapseclient

        syn = synapseclient.Synapse()
        syn.login(authToken=os.environ.get("SYNAPSE_AUTH_TOKEN"), silent=True)
        return package_from_synapse(syn, args.synapse_id, args.output, args.upload_to)

    if args.scan:
        if not args.scan.is_dir():
            print(f"{args.scan}: not a folder", file=sys.stderr)
            return 2
        sets, problems = find_local_sets(args.scan)
        args.output.mkdir(parents=True, exist_ok=True)
        for name, files in sets.items():
            print(package_matrix(files["barcodes"], files["features"], files["matrix"], args.output / name))
        for problem in problems:
            print(f"skipped {problem}", file=sys.stderr)
        if not sets and not problems:
            print(f"No barcodes/features/matrix files found under {args.scan}", file=sys.stderr)
        return 1 if problems else 0

    if len(args.files) != 3:
        print(f"Expected 3 files (barcodes, features, matrix), --synapse-id or --scan, got {len(args.files)} file(s)",
              file=sys.stderr)
        return 2
    if not args.output.name.endswith((".tar.gz", ".tgz")):
        print(f"Output name must end in .tar.gz or .tgz: {args.output}", file=sys.stderr)
        return 2
    try:
        package_matrix(*args.files, args.output)
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 2
    print(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
