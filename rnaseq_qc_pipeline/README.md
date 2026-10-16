# RNA-seq QC pipeline

`rnaseq_qc_metrics.py` summarizes QC metrics from an nf-core/rnaseq run into one TSV row per
sample. Every output column is defined in
[rnaseq_qc_metrics_data_dictionary.md](rnaseq_qc_metrics_data_dictionary.md).

```
python rnaseq_qc_metrics.py --rnaseq-dir <nfcore_outdir> --samplesheet <samplesheet.csv> --out qc_metrics.tsv
```

Requirements: Python 3 with `numpy` and `scipy`. No other tools are needed, because the script only
reads files that nf-core/rnaseq already produced. The `namhub_atac_qc` conda environment
(`../atacseq_qc_pipeline/environment.yml`) includes both packages.

`--samplesheet` is the `samplesheet.csv` that was given to nf-core/rnaseq (`--input`) for this run:
the same file. Only its `sample` column is read (the list of samples; replicate rows share a
sample). It is not the replicate sheet below. Without it, the samples are those with a STAR
`Log.final.out`.

Options: `--replicate-sheet` (an optional, separate CSV with `sample, group, donor`; not the nf-core samplesheet. Same donor = isogenic pair, different donors = anisogenic; adds the Spearman replicate-concordance columns and a `_replicate_pairs.tsv` file; see the data dictionary), `--aligner` (default `star_salmon`), `--species` (expected species for the
contamination check, default `Homo sapiens`), `--samples` (restrict to some samples).

The fq-lint, read length, FastQC and contamination sections live in
[`../qc_common/fastq_sections.py`](../qc_common/fastq_sections.py) so other pipelines can reuse them;
the STAR, Qualimap, featureCounts, rRNA and Picard sections stay in the script.

## Known limitation: aligners other than star_salmon

The script has only been tested with the `star_salmon` aligner/quantifier output.

- Passing `--samplesheet` means sample names no longer come from STAR `Log.final.out` files,
  so sample discovery will not fail for other aligners.
- However, the STAR, Qualimap, featureCounts, Picard and Salmon sections still read
  `star_salmon`-style files under `<rnaseq-dir>/<aligner>/`. For other aligners these columns
  may come out blank (a warning is printed per sample and section) or need code changes.
- Other aligners have not been tested yet. Fixes will be made once they are.

## Test data

The script has been tested on one public dataset so far:

- **[GSE195688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195688)** (PubMed
  [36920035](https://pubmed.ncbi.nlm.nih.gov/36920035/)): "Lineage-Specific Differences and
  Inference of Regulatory Networks Governing Human Chondrocyte Development". Human, in-vitro
  chondrocytes derived from H9 cells. Four paired-end RNA-seq samples were used, 2 treatment
  groups x 2 replicates, each with 4 sequencing runs that nf-core merges:

  | `sample` | Group | GEO sample (title) |
  |---|---|---|
  | SRX13980357 | TGF-beta | GSM5848444 (H9_T_44_rep1) |
  | SRX13980358 | TGF-beta | GSM5848445 (H9_T_44_rep2) |
  | SRX13980360 | BMP | GSM5848447 (H9_B_44_rep1) |
  | SRX13980361 | BMP | GSM5848448 (H9_B_44_rep2) |

  The FASTQs were downloaded with nf-core/fetchngs, and nf-core/rnaseq was run with the
  `star_salmon` aligner on Ensembl GRCh38 (release 116). The replicate sheet for this test is
  [`replicate_metadata_GSE195688.csv`](replicate_metadata_GSE195688.csv) (an input, passed with
  `--replicate-sheet`).

**Example output.** [`test_outputs/qc_metrics_GSE195688.tsv`](test_outputs/qc_metrics_GSE195688.tsv)
is the table `rnaseq_qc_metrics.py` wrote for this test (one row per sample). It is a reference copy
for seeing what the output looks like; the script never writes there.

More testing is planned (more samples of this dataset, other aligners); see the limitation above.
