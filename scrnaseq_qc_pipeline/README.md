# scRNA-seq QC Pipeline

An example of what QC for single-cell RNA-seq data can look like when a
dataset originates from a public archive (GEO/SRA) but is also tracked in
[Synapse](https://www.synapse.org) (a data-sharing platform used by several
biomedical research consortia). It was built and tested against one real
public dataset — [GSE293390](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293390),
tracked in Synapse as project `syn73675040` — and is at a preliminary,
exploratory stage: it's a demonstration of an approach, not a finished,
validated gate. See "A note on preliminary status" below.

## Background, for readers new to this kind of data

- **GEO** (Gene Expression Omnibus) is NCBI's public archive for functional
  genomics data — where researchers deposit sequencing datasets alongside
  their metadata (species, tissue, sequencing platform, etc.) so others can
  reuse them. A dataset there is a **GSE** (series) accession; an individual
  sample within it is a **GSM** accession.
- **SRA**/**ENA** (Sequence Read Archive / European Nucleotide Archive) are
  where the actual raw sequencing reads behind a GEO submission live — GEO
  itself usually only stores processed results and metadata, not the raw
  reads. A single sequencing run is an **SRR** accession. ENA mirrors SRA's
  data and — helpfully for this pipeline — makes it available as plain
  `.fastq.gz` files over plain HTTPS, where SRA's own format requires
  specialized tooling to convert.
- A **FASTQ** file holds raw sequencing reads: one entry per read, with the
  sequence and a per-base quality score. Many sequencing runs are
  **paired-end** — each physical DNA/RNA fragment is read from both ends,
  producing two mates per fragment (**R1** and **R2**), either as two
  separate files or interleaved together in one file.
- **CellRanger** is 10x Genomics' standard software for processing single-cell
  RNA-seq reads into a cell-by-gene count matrix (who — which cell barcode —
  expressed how much of which gene). It outputs a **raw** (unfiltered, every
  technically-possible cell barcode included, the overwhelming majority of
  them empty droplets) and/or a **filtered** (CellRanger's own best guess at
  which barcodes are real cells) version of that matrix, stored as three
  files: `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`. Synapse
  records which version a given matrix is via a `CellrangerOutputClass`
  annotation (`raw_feature_bc_matrix` or `filtered_feature_bc_matrix`).
- Synapse's own curation tooling validates that annotations like platform or
  species are legal values from a controlled vocabulary at data-entry time —
  it does not, and cannot, check that a curator picked the value that's
  actually *correct* for this specific sample. That's the gap much of this
  pipeline exists to help catch.

## What this pipeline actually does

It walks a Synapse project/folder, classifies every file as either a raw
sequencing run or a CellRanger-style processed sample, and runs two kinds of
check:

- **External data checked against its own external record** (the main focus
  of this pipeline) — curator-entered Synapse annotations are cross-checked
  against GEO's own record for that sample, and downloaded FASTQ content is
  cross-checked against ENA's record for that same run. Neither side is ever
  checked against the Synapse curation schema (that's the curator tooling's
  job already, described above), and Synapse and GEO/ENA are never
  conflated — each side is independently checked against its own external
  source.
- **Actual data content QC** — does the sequencing data itself look right,
  independent of what any metadata claims: read quality (FastQC), structural
  validity (`fq lint`), species/contamination screening (Kraken2), and
  processed-matrix sanity checks.

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
| `--skip-fastqc` | Skip the whole raw-read tier (`fq lint` + FastQC + Kraken2 + ENA integrity checks) |
| `--skip-matrix` | Skip the processed CellRanger matrix tier |
| `--keep-downloads` | Keep downloaded fastq/matrix files instead of deleting them after QC |

## Dependencies

Python: `synapseclient`, `numpy`, `scipy`, `requests` (see repo-root `requirements.txt`).

External tools, none of which are pip-installable:

| Tool | What it is | Install | Used for |
|---|---|---|---|
| `fastqc` | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) — the standard, widely-used tool for assessing raw sequencing read quality | `brew install fastqc` / `conda install -c bioconda fastqc` | Content quality QC |
| `fq` | A FASTQ file structural validator, by [stjude-rust-labs](https://github.com/stjude-rust-labs/fq) | Download a release binary from [github.com/stjude-rust-labs/fq/releases](https://github.com/stjude-rust-labs/fq/releases) and put it on PATH. No native ARM64 macOS build exists — the `x86_64-apple-darwin` build runs fine under Rosetta 2. | `fq lint` structural validation |
| `kraken2` | [Kraken2](https://github.com/DerrickWood/kraken2) — classifies sequencing reads against a reference database to identify what organism(s)/contaminants they actually came from | `brew install kraken2` | Species/contamination screening |
| Kraken2 database | A pre-built reference index Kraken2 classifies reads against | Pre-built options at [benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2) — "Standard-8" (~5.5GB compressed; a reduced-size index built from NCBI's RefSeq archaea/bacteria/viral/plasmid/human sequences plus UniVec, a database of common vector and adapter sequences) is a reasonable general-purpose default. Extract it and pass the directory via `--kraken2-db`. | Species/contamination screening |
| `geo-synapse` | A companion Python package (from the [`geo_dataset_creation`](https://github.com/sagebio-ada/geo_dataset_creation) repo) that resolves GEO/SRA accessions to real download links and fetches GEO metadata | `pip install git+https://github.com/sagebio-ada/geo_dataset_creation.git` | Resolving SRA run accessions to real ENA fastq download links, and fetching GEO sample metadata |

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
| `geo_metadata_qc.py` | Synapse annotations vs. GEO's own record |
| `fastq_qc.py` | Resolves SRA→ENA fastq links, downloads, checks integrity against ENA's record, calls `fq_lint_qc` and `kraken_qc`, runs FastQC |
| `fq_lint_qc.py` | Wraps `fq lint` for FASTQ structural validation |
| `kraken_qc.py` | Wraps Kraken2 for species/contamination screening |
| `matrix_qc.py` | CellRanger matrix.mtx.gz + barcodes/features sanity checks |
| `qc_common.py` | Shared finding format, CSV report writer |

## A note on preliminary status

This pipeline is a demonstration of an approach, built and tested against one
real dataset, not a finished, broadly-validated tool. Two things follow from
that, both reflected in the report format and worth keeping in mind if you're
adapting this for your own data:

- **This pipeline's own code never asserts a `PASS`/`FAIL` verdict about the
  data.** Only two checks carry a real pass/warn/fail judgment at all —
  `fq lint`'s exit code and FastQC's own built-in per-module thresholds —
  because those are a named tool's own established judgment, not this
  pipeline's opinion. Every comparison this pipeline's own code performs
  (does a checksum match, does a count match, does a Synapse annotation
  match GEO's record) reports `INFO`: both sides of the comparison, with no
  verdict attached, even where the comparison is an exact one. `FAIL`
  appears elsewhere only for a genuine operational failure — a network
  error, a file that won't parse, a required file that's missing — where
  there's no data to report at all, not a judgment about data quality.
  See "Metrics generated" below for exactly which checks fall into which
  category.
- **The report leads with the number, not the verdict.** Columns are ordered
  `sample, check, metric, status, entity_id`, and rows are sorted by
  `sample`/`check` rather than by severity — reading the actual data is the
  point at this stage, not triaging by pass/fail.

## Metrics generated

Every row in `scrnaseq_qc_report.csv` has `sample`, `check`, `metric` (the
actual reported number/fact), `status` (`PASS`/`WARN`/`FAIL`/`INFO`), and
`entity_id` (the specific Synapse file or sample this row is about). 26
distinct check types in the code (the table below lists 29 rows — `fastqc_module:<name>`
is one dynamic check type that names a different FastQC module each time it
fires, split across 4 rows here for the sake of documenting its most common
and most interpretation-worthy names individually), grouped below by
**what kind of check it is**:

- **Metadata cross-check** — compares a declared/registered value against
  another source: a Synapse annotation vs. GEO's record, or a file's own
  measured property vs. ENA's record.
- **File-content check** — pure inspection of the actual file bytes; nothing
  external to compare against.
- **Existence/resolution check** — did we successfully find or fetch the
  thing at all.

Two checks are genuinely hybrid and got a judgment call rather than a clean
fit: `kraken2_species_match` inspects file content *and* compares the result
to GEO's declared species, and `matrix_dimensions` compares two
independently-derived numbers (the matrix file's own header vs. counted
barcode/feature lines). Both are filed under Metadata cross-check since the
comparison is the point of the check, not just the inspection.

**"Full file needed?"** — `Yes` means the check only runs against a fully
downloaded raw fastq, not a `--max-download-mb` spot-check sample (mostly the
checks that need an exact byte-for-byte or read-for-read comparison against
ENA's record). `No` means the check works identically on a small capped
sample — these are the ones worth running by default, reserving full
downloads for when exact integrity verification specifically matters.
`No file at all` means the check never touches file content at all — it's a
metadata comparison or a Synapse/API existence lookup.

**"Tool"** names the actual code that runs the check — a named open-source
tool if one is doing real work, `synapseclient`/`geo-synapse` when the check
is really just an API lookup, or "Custom code" when nothing but this
pipeline's own logic is involved (numpy/scipy noted where they do the
numeric heavy lifting). Nothing from nf-core or `py-dcqc` actually runs here
— both were only researched as references for which standalone tools were
worth adding (that's how Kraken2 and `fq` ended up in this pipeline).
**"Compared against"** names the external source a result is checked against,
where the check involves one at all.

**"Verdict basis"** — only two rows in this table can ever produce a real
`PASS`/`WARN`/`FAIL` judgment: `fq lint`'s exit code and FastQC's own
built-in per-module thresholds, both a named tool's own established
judgment, not this pipeline's opinion. Every other row is `INFO` on success —
reporting both sides of a comparison with no verdict attached, even for an
exact comparison like a checksum — and `FAIL` only for a genuine operational
failure (an exception: a network error, a file that won't parse, a required
file that's missing), never as a graded judgment about the data itself.

The **"Interpreting this"** column calls out where a check's result needs
domain context to read correctly — several of FastQC's and Kraken2's stock
thresholds are calibrated for generic sequencing libraries and routinely
fire on artifacts that are normal, expected byproducts of single-cell 3'
RNA-seq specifically, not signs of bad data. A dash means the check is
straightforward to read as-is.

| Category | Check | What it reports | Tool | Compared against | Full file needed? | Interpreting this |
|---|---|---|---|---|---|---|
| Metadata cross-check | `platform_vs_geo` | Synapse's `platform` annotation vs. GEO's own recorded sequencing instrument, after normalizing formatting differences between the two | `geo-synapse` | GEO's record | No file at all | — |
| Metadata cross-check | `reference_vs_geo` | Synapse's `referenceSet` (reference genome) annotation vs. GEO's recorded genome assembly | `geo-synapse` | GEO's record | No file at all | — |
| Metadata cross-check | `library_version_vs_geo` | Synapse's `libraryVersion` (e.g. 10x Chromium chemistry version) annotation vs. GEO's own free-text description of the sequencing kit used | `geo-synapse` | GEO's record | No file at all | — |
| Metadata cross-check | `library_prep_method_vs_geo` | Synapse's `libraryPreparationMethod` annotation vs. GEO's own free-text kit description | `geo-synapse` | GEO's record | No file at all | — |
| Metadata cross-check | `geo_lookup` | Whether the sample accession was found in the given GEO series at all | `geo-synapse` | GEO's record | No file at all | — |
| Metadata cross-check | `fastq_size_vs_ena` | Downloaded file's byte size vs. the size ENA has registered for it | Custom code | ENA's record | **Yes** | — |
| Metadata cross-check | `fastq_md5_vs_ena` | Downloaded file's checksum vs. the checksum ENA has registered for it | Custom code (MD5) | ENA's record | **Yes** | — |
| Metadata cross-check | `read_count_vs_ena` | FastQC's total-read count vs. the read count ENA has registered for this run | FastQC | ENA's record | **Yes** | — |
| Metadata cross-check | `mean_length_vs_registry` | The read length ENA's own read/base counts imply, alongside FastQC's directly observed read length(s) | FastQC | ENA's record | **Yes** | A single uniform length here, rather than two distinct lengths for a barcode read and a cDNA read, may mean the interleaved-pair structure implied by `pairedEnd`/`libraryVersion` isn't actually present in the public deposit — a data-provenance question for the original submitter, not something this pipeline resolves. |
| Metadata cross-check | `paired_fastq_parity` | For paired-end runs stored as two separate mate files (not interleaved into one), the read count of each mate | FastQC | — (compares FastQC's own per-file counts to each other) | **Yes** | — |
| Metadata cross-check | `matrix_dimensions` | The matrix file's own header dimensions, alongside the barcode/feature counts counted from the other two files | Custom code (`scipy`/`numpy`) | — (compares two internally-derived counts) | **Yes** | — |
| Metadata cross-check | `kraken2_species_match` | % of reads Kraken2 actually classified as the organism GEO says this sample is | Kraken2 | GEO's record | No | A low percentage isn't necessarily concerning for RNA-seq — see `kraken2_unclassified` below. |
| File-content check | `fq_lint` | Structural integrity of the FASTQ file itself — complete records, valid sequence alphabet, correctly formatted separator line, matching sequence/quality-score lengths, well-formed quality scores | [`fq`](https://github.com/stjude-rust-labs/fq) | — | No | — |
| File-content check | `fastqc_basic_stats` | Total read count, read length, overall %GC | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | — | No | — |
| File-content check | `fastqc_module:Per base sequence content` | Whether A/T and G/C percentages diverge too much at any read position | FastQC | — | No | 10x libraries are built via random-hexamer priming and template-switching, which creates a well-documented, expected sequence bias in the first ~20 bases specifically. A WARN/FAIL here often isn't a sign of bad data for this assay type. |
| File-content check | `fastqc_module:Adapter Content` | Whether a known adapter sequence shows up in too many reads at any position | FastQC | — | No | Check *which* adapter triggered it before treating this as a problem. PolyA is very often the transcript's own genuine poly-A tail — this assay captures mRNA via poly-T primers that bind to that same tail — not synthetic library contamination, especially for reads that reach a short transcript's 3' end. |
| File-content check | `fastqc_module:Sequence Duplication Levels` | What fraction of reads are exact duplicates | FastQC | — | No | RNA-seq — especially PCR-amplified single-cell cDNA — naturally has far higher duplication than genomic DNA-seq, since highly-expressed transcripts are legitimately sequenced many times over. A WARN/FAIL here doesn't carry the same weight it would for DNA-seq. |
| File-content check | `fastqc_module:<other names>` | Any other FastQC module that isn't a clean pass — per-sequence quality scores, per-sequence GC content, per-base ambiguous-base (N) content, read length distribution, overrepresented sequences, repetitive k-mer content. Only relays FastQC's verdict word for these; the three modules above get FastQC's underlying position/adapter/percentage extracted too | FastQC | — | No | — |
| File-content check | `kraken2_unclassified` | % of reads Kraken2 couldn't confidently classify against the reference database at all | [Kraken2](https://github.com/DerrickWood/kraken2) | — | No | Expect a much higher rate here than for DNA-seq: the reference database is genome-based, not transcript/splice-aware, so reads spanning splice junctions or from lightly-annotated transcripts often don't classify even when the sample and organism are entirely fine. |
| File-content check | `kraken2_top_species` | The top 5 organisms Kraken2 actually found in the reads, and what fraction of reads matched each | Kraken2 | — | No | — |
| File-content check | `file_integrity` | Whether the barcode and feature list files are readable, with their line counts | Custom code | — | **Yes** | — |
| File-content check | `matrix_integrity` | Whether the count-matrix file parses as a valid matrix, with its parsed shape | Custom code (`scipy`) | — | **Yes** | — |
| File-content check | `matrix_metrics` | Median UMI (unique transcript molecule) count per cell barcode, median genes detected per cell barcode, how many barcodes have zero counts, total nonzero matrix entries | Custom code (`numpy`) | — | **Yes** | — |
| File-content check | `empty_barcode_check` | What fraction of barcodes have zero counts, reported alongside whether Synapse declares this the raw (unfiltered) or filtered CellRanger output | Custom code (`numpy`) | — | **Yes** | For this dataset specifically, both matrices show 0.0% empty barcodes despite being declared `raw_feature_bc_matrix` — genuine unfiltered 10x output should have millions of largely-empty barcodes, so this looks like a mislabeled filtered matrix rather than a clean result. |
| File-content check | `mito_fraction` | Median fraction of each cell's counts coming from mitochondrial genes (a common per-cell quality signal — a very high fraction often indicates a dying/stressed cell), among barcodes with any counts at all | Custom code (`numpy`) | — | **Yes** | — |
| Existence/resolution check | `fastq_resolution` | Whether the archived sequencing run's accession successfully resolved to a real, directly downloadable file, with the actual resolved URL(s) | `geo-synapse` | — | No file at all (API lookup only) | — |
| Existence/resolution check | `fastq_download` | Whether the download itself succeeded, with the filename and size | Custom code | — | No — reports whatever was requested, capped or full | — |
| Existence/resolution check | `processed_completeness` | Whether all three required files (barcodes, features, matrix) are present for this sample, with the specific Synapse file ID behind each one | `synapseclient` | — | No file at all | — |
| Existence/resolution check | `download:<role>` | Whether each of the three processed-data files downloaded successfully, with its Synapse file ID, filename, and size | `synapseclient` | — | **Yes** — no capped-download option exists for these files; a full download always happens | — |
