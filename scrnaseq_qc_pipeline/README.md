# scRNA-seq QC Pipeline

QC pipeline for single-cell RNA-seq files living under a Synapse project/folder
(built against `syn73675040`, GSE293390). Walks the Synapse tree, classifies
every file as a raw sequencing run or a CellRanger-style processed sample, and
runs a battery of checks split into two kinds:

- **External data vs. its own external registry/provenance** (ENA, GEO) — the
  main focus of this pipeline. Curator-entered Synapse annotations are cross-
  checked against GEO's own SOFT record, and downloaded FASTQ content is
  cross-checked against ENA's file report — never against each other, and
  never against the Synapse curation schema (that's already enforced by the
  curator tooling at entry time, so this pipeline doesn't re-check it).
- **Actual data content QC** — FastQC, `fq lint`, Kraken2 species/contamination
  screening, and CellRanger matrix sanity checks.

## Usage

```bash
python run_scrnaseq_qc.py --synapse-id syn73675040 --output-dir output/qc

# Also cross-check curated annotations against the GEO record:
python run_scrnaseq_qc.py --synapse-id syn73675040 --geo-accession GSE293390

# Fast pass, no FastQC / no matrix download:
python run_scrnaseq_qc.py --synapse-id syn73675040 --skip-fastqc --skip-matrix

# Spot-check raw reads without downloading full (multi-GB) fastq files:
python run_scrnaseq_qc.py --synapse-id syn73675040 --max-download-mb 25

# Also screen for species/contamination (needs a local Kraken2 database):
python run_scrnaseq_qc.py --synapse-id syn73675040 --geo-accession GSE293390 \
    --max-download-mb 25 --kraken2-db /path/to/k2_standard_08gb
```

| Flag | Purpose |
|---|---|
| `--synapse-id` | Synapse project or folder ID to scan (required) |
| `--geo-accession` | GEO series (e.g. `GSE293390`) to cross-check curated annotations against. Omit to skip that tier. |
| `--output-dir`, `-d` | Where `scrnaseq_qc_report.csv` is written (default `output`) |
| `--work-dir` | Where downloaded fastq/matrix files land (default `<output-dir>/scrnaseq_qc_work`) |
| `--max-download-mb` | Cap raw fastq downloads for a fast spot-check instead of downloading the full file. Strongly recommended whenever `--kraken2-db` is set — classification is slow on full multi-hundred-million-read files. |
| `--kraken2-db` | Path to an extracted Kraken2 database. If omitted, species/contamination screening is skipped entirely. |
| `--kraken2-threads` | Threads for Kraken2 classification (default 4) |
| `--skip-fastqc` | Skip the whole raw-read tier (fq lint + FastQC + Kraken2 + ENA integrity checks) |
| `--skip-matrix` | Skip the processed CellRanger matrix tier |
| `--keep-downloads` | Keep downloaded fastq/matrix files instead of deleting them after QC |

## Dependencies

Python: `synapseclient`, `numpy`, `scipy`, `requests` (see repo-root `requirements.txt`).

External tools, none of which are pip-installable:

| Tool | Install | Used for |
|---|---|---|
| `fastqc` | `brew install fastqc` / `conda install -c bioconda fastqc` | Content quality QC |
| `fq` (stjude-rust-labs) | Download a release binary from [github.com/stjude-rust-labs/fq/releases](https://github.com/stjude-rust-labs/fq/releases) and put it on PATH. No native ARM64 macOS build exists — the `x86_64-apple-darwin` build runs fine under Rosetta 2. | `fq lint` structural validation |
| `kraken2` | `brew install kraken2` | Species/contamination screening |
| Kraken2 database | Pre-built databases at [benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2) — "Standard-8" (~5.5GB compressed; RefSeq archaea/bacteria/viral/plasmid/human/UniVec) is a reasonable general-purpose default. Extract it and pass the directory via `--kraken2-db`. | Species/contamination screening |
| `geo-synapse` | `pip install git+https://github.com/sagebio-ada/geo_dataset_creation.git` | Resolving SRA run accessions to real ENA fastq download links, and fetching GEO sample metadata |

**Note:** installing `kraken2` via Homebrew pulls in `python@3.14` as a
dependency, which Homebrew may link as the generic `python3` — silently
shadowing whatever Python previously resolved to `python3` (e.g. the system
Python with your other dependencies already installed). If your dependencies
suddenly go missing after installing `kraken2`, check `which python3` and
`brew unlink python@3.14` if needed (`kraken2` itself is a Perl script with no
Python runtime dependency, so this is safe).

## Architecture

| Module | Role |
|---|---|
| `run_scrnaseq_qc.py` | Orchestrator/CLI: walks the Synapse tree, classifies files, wires the other modules together, writes the report |
| `geo_metadata_qc.py` | Synapse annotations vs. GEO's own SOFT record |
| `fastq_qc.py` | Resolves SRA→ENA fastq links, downloads, checks integrity against ENA's registry, calls `fq_lint_qc` and `kraken_qc`, runs FastQC |
| `fq_lint_qc.py` | Wraps `fq lint` for FASTQ structural validation |
| `kraken_qc.py` | Wraps Kraken2 for species/contamination screening |
| `matrix_qc.py` | CellRanger matrix.mtx.gz + barcodes/features sanity checks |
| `qc_common.py` | Shared finding format, CSV report writer |

## Metrics generated

Every row in `scrnaseq_qc_report.csv` has `entity_id`, `sample`, `check`,
`status` (PASS/WARN/FAIL), and `detail`. 27 distinct check types:

- **"Full file needed?" = Yes** means the check is skipped (or, for the ENA
  integrity checks, silently not computed) unless the raw fastq was
  downloaded in full — it's gated on `not any_capped` in the code, so a
  `--max-download-mb` spot-check never produces these rows.
- **No** means the check runs identically on a `--max-download-mb`-capped
  sample or the full file — these are the checks worth running by default for
  routine QC, reserving full downloads for when you specifically need
  integrity verification.
- **No file at all** means the check never touches file content — it's pure
  Synapse-annotation/GEO-metadata comparison or Synapse entity/API lookups.

| Tier | Check | What it verifies | Tool / repo | Full file needed? |
|---|---|---|---|---|
| Synapse vs. GEO | `platform_vs_geo` | Synapse `platform` vs GEO `sample_instrument_model` (normalized) | `geo_synapse.geo` (geo_dataset_creation repo) for the GEO fetch; comparison is this pipeline's own code | No file at all |
| Synapse vs. GEO | `reference_vs_geo` | Synapse `referenceSet` vs GEO `assembly` | `geo_synapse.geo` (geo_dataset_creation repo) | No file at all |
| Synapse vs. GEO | `library_version_vs_geo` | Synapse `libraryVersion` vs GEO's reagent-kit description text | `geo_synapse.geo` (geo_dataset_creation repo) | No file at all |
| Synapse vs. GEO | `library_prep_method_vs_geo` | Synapse `libraryPreparationMethod` vs GEO's reagent-kit description text | `geo_synapse.geo` (geo_dataset_creation repo) | No file at all |
| Synapse vs. GEO | `geo_lookup` | Flags a raw run whose GSM isn't found in the given GEO series at all | `geo_synapse.geo` (geo_dataset_creation repo) | No file at all |
| Raw FASTQ vs. ENA | `fastq_resolution` | SRA run accession successfully resolved to real ENA fastq.gz link(s) | `geo_synapse.ena` (geo_dataset_creation repo) | No file at all (API lookup only) |
| Raw FASTQ vs. ENA | `fastq_download` | Download itself succeeded | This pipeline's own code (`requests`) | No — reports whatever was requested, capped or full |
| Raw FASTQ vs. ENA | `fastq_size_vs_ena` | Downloaded byte size vs. ENA's registered `fastq_bytes` | This pipeline's own code vs. ENA registry | **Yes** |
| Raw FASTQ vs. ENA | `fastq_md5_vs_ena` | Downloaded md5 vs. ENA's registered `fastq_md5` | This pipeline's own code (`hashlib`) vs. ENA registry | **Yes** |
| Raw FASTQ vs. ENA | `read_count_vs_ena` | FastQC's total sequence count (summed across files) vs. ENA's registered `read_count` | FastQC + ENA registry | **Yes** |
| Raw FASTQ vs. ENA | `mean_length_vs_registry` | ENA `base_count/read_count`-implied mean length vs. FastQC's observed length(s) | FastQC + ENA registry | **Yes** |
| Raw FASTQ vs. ENA | `paired_fastq_parity` | Split R1/R2 mate files have matching read counts (skipped for a single interleaved file) | FastQC per-file counts | **Yes** |
| FASTQ structure | `fq_lint` | Structural integrity — record completeness, valid alphabet, `+` line, matching seq/quality lengths, well-formed quality string | [`fq`](https://github.com/stjude-rust-labs/fq) (stjude-rust-labs) | No |
| FastQC content | `fastqc_basic_stats` | Total sequences, sequence length, %GC summary | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | No |
| FastQC content | `fastqc_module:<name>` | Any FastQC module that isn't a clean PASS — per-base quality, per-sequence quality, per-base/sequence GC content, per-base N content, sequence length distribution, duplication levels, overrepresented sequences, adapter content, k-mer content | FastQC | No |
| Species/contamination | `kraken2_unclassified` | % of reads not classified against the reference database | [Kraken2](https://github.com/DerrickWood/kraken2) | No |
| Species/contamination | `kraken2_top_species` | Top 5 classified species and their read fractions | Kraken2 | No |
| Species/contamination | `kraken2_species_match` | % of reads matching the GEO-declared organism | Kraken2 + `geo_synapse.geo` for expected species | No |
| Species/contamination | `kraken2_contaminant` | Any other species above 2% of classified reads, flagged as potential contamination | Kraken2 | No |
| Processed matrix | `processed_completeness` | barcodes/features/matrix trio all present for the sample | `synapseclient` (entity check only) | No file at all |
| Processed matrix | `download:<role>` | Per-file download succeeded | `synapseclient` | **Yes** — no capping mechanism exists for matrix files; `syn.get()` always pulls the whole file |
| Processed matrix | `file_integrity` | barcodes/features files are readable; reports counts | This pipeline's own code (`gzip`) | **Yes** |
| Processed matrix | `matrix_integrity` | `matrix.mtx.gz` parses without error | `scipy.io.mmread` | **Yes** |
| Processed matrix | `matrix_dimensions` | Matrix shape matches barcode/feature counts | `scipy`/`numpy` | **Yes** |
| Processed matrix | `matrix_metrics` | Median UMI/cell, median genes/cell, zero-count barcodes, nnz | `numpy` | **Yes** |
| Processed matrix | `empty_barcode_check` | Empty-barcode fraction is sane given the declared raw-vs-filtered `CellrangerOutputClass` | `numpy` | **Yes** |
| Processed matrix | `mito_fraction` | Median mitochondrial read fraction across non-empty barcodes | `numpy` | **Yes** |
