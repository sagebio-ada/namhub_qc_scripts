#!/bin/bash
# Full ATAC-seq QC for one nf-core/atacseq run, as separate Slurm jobs that run at the same time:
#   encode    atacseq_encode_qc_metrics.py (ENCODE metrics of the merged libraries)
#   idr_<g>   atacseq_idr_concordance.py, one job per group with >= 2 replicates
#   tools     atacseq_tools_qc.py, a job array with one task per library (fq lint, kraken2/Bracken)
#   table     atacseq_qc_metrics.py, starts when all of the above have finished successfully
# Run it from any node (it only calls sbatch):
#   bash submit_atacseq_qc.sh                        # output in $ATAC_DIR/qc_metrics
#   QC_DIR=/path/to/out bash submit_atacseq_qc.sh    # output elsewhere
# Edit the paths below. Rerunning resumes: finished tool output and IDR steps are skipped.
# Logs: $QC_DIR/logs/. If one job fails the table job is cancelled (afterok); fix and resubmit.
# Needs the conda env in environment.yml.

set -euo pipefail

ATAC_DIR=/gpfs/data/proteomics/home/snk218/namhub/nfcore/atacseq/GSE195688_reps
SAMPLESHEET=/gpfs/data/proteomics/home/snk218/namhub/nfcore/fetchngs/GSE195688/SRR_Acc_list_AtacSeq/samplesheet/samplesheet_fixed_reps.csv
KRAKEN_DB=/gpfs/data/proteomics/home/snk218/namhub/nfcore/ensembl/hsapiens/index/k2_pluspf_16_GB_20260626
KRAKEN_MAX_READS=0          # reads (pairs) screened per file; 0 = all (as nf-core/rnaseq); 1000000 for a quick run
MITO=MT                     # mitochondrial contig name (UCSC: chrM)
GENOME_FAI=$ATAC_DIR/genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai   # genome index used for the run
KEEP_CHR_REGEX='^([0-9]+|X|Y)$'   # IDR peaks kept on these chromosomes (UCSC: '^chr([0-9]+|X|Y)$')
# IDR/overlap peaks overlapping this are dropped, as ENCODE does; use the blacklist given to nf-core
# (its params.yaml). Empty = no blacklist filter.
BLACKLIST=/gpfs/data/proteomics/home/snk218/namhub/nfcore/blacklists/v3.0/hg38-blacklist.v3.ensembl.bed
SCRIPTS=/gpfs/data/proteomics/home/snk218/namhub/namhub_qc_scripts/atacseq_qc_pipeline
QC_DIR=${QC_DIR:-$ATAC_DIR/qc_metrics}
PARTITION=cpu_short,cpu_medium
CONDA_SH=/gpfs/data/proteomics/home/snk218/software/miniforge3/etc/profile.d/conda.sh

# Per-job resources (CPUs, memory, time limit)
ENCODE_RES="-c 12 --mem=32G -t 12:00:00"
IDR_RES="-c 12 --mem=32G -t 12:00:00"       # per group
TOOLS_RES="-c 8 --mem=32G -t 12:00:00"      # per library; kraken2 holds the 16 GB database in memory
TABLE_RES="-c 2 --mem=8G -t 01:00:00"

source "$CONDA_SH"
conda activate namhub_atac_qc
mkdir -p "$QC_DIR/logs" "$QC_DIR/idr"

BAM_DIR=$ATAC_DIR/bwa/merged_library
# groups with >= 2 replicates: BAM names are <group>_REP<n>.mLb.clN.sorted.bam
IDR_GROUPS=$(for b in "$BAM_DIR"/*_REP*.mLb.clN.sorted.bam; do
                 basename "$b" | sed -E 's/_REP[0-9]+\.mLb\.clN\.sorted\.bam$//'
             done | sort | uniq -c | awk '$1 >= 2 {print $2}')
N_LIBS=$(python "$SCRIPTS/atacseq_tools_qc.py" --samplesheet "$SAMPLESHEET" --list-libraries | wc -l)

# The jobs inherit these variables (sbatch exports the environment), so their scripts are plain text
export ATAC_DIR SAMPLESHEET KRAKEN_DB KRAKEN_MAX_READS MITO GENOME_FAI KEEP_CHR_REGEX BLACKLIST \
       SCRIPTS QC_DIR BAM_DIR CONDA_SH
COMMON=(--parsable --partition="$PARTITION")

ENC=$(sbatch "${COMMON[@]}" --job-name=atac_encode $ENCODE_RES --output="$QC_DIR/logs/encode_%j.log" <<'EOF'
#!/bin/bash
set -euo pipefail
source "$CONDA_SH"; conda activate namhub_atac_qc
python "$SCRIPTS/atacseq_encode_qc_metrics.py" \
    --bam-dir "$BAM_DIR" \
    --out "$QC_DIR/atacseq_encode_qc_metrics.tsv" \
    --mito-name "$MITO" --threads 2 --jobs 6
EOF
)
echo "encode: job $ENC"
DEPS="$ENC"

for g in $IDR_GROUPS; do
    mkdir -p "$QC_DIR/idr/$g"
    J=$(GROUP=$g sbatch "${COMMON[@]}" --job-name="atac_idr_$g" $IDR_RES \
        --output="$QC_DIR/logs/idr_${g}_%j.log" <<'EOF'
#!/bin/bash
set -euo pipefail
source "$CONDA_SH"; conda activate namhub_atac_qc
OUT="$QC_DIR/idr/$GROUP"
cut -f1,2 "$GENOME_FAI" > "$OUT/chrom.sizes"
python "$SCRIPTS/atacseq_idr_concordance.py" \
    --bam-dir "$BAM_DIR" \
    --out-dir "$OUT" --groups "$GROUP" \
    --chrom-sizes "$OUT/chrom.sizes" \
    --keep-chr-regex "$KEEP_CHR_REGEX" \
    ${BLACKLIST:+--blacklist "$BLACKLIST"} \
    --macs2 "$SCRIPTS/macs2_cluster.sh" \
    --threads 4 --jobs 4
EOF
)
    echo "idr $g: job $J"
    DEPS="$DEPS:$J"
done
[ -n "$IDR_GROUPS" ] || echo "No group with >= 2 replicates: no IDR jobs, IDR columns will be blank"

TOOLS=$(sbatch "${COMMON[@]}" --job-name=atac_tools $TOOLS_RES --array=1-"$N_LIBS" \
        --output="$QC_DIR/logs/tools_%A_%a.log" <<'EOF'
#!/bin/bash
set -euo pipefail
source "$CONDA_SH"; conda activate namhub_atac_qc
LIB=$(python "$SCRIPTS/atacseq_tools_qc.py" --samplesheet "$SAMPLESHEET" --list-libraries \
      | sed -n "${SLURM_ARRAY_TASK_ID}p")
python "$SCRIPTS/atacseq_tools_qc.py" \
    --atac-dir "$ATAC_DIR" --samplesheet "$SAMPLESHEET" --qc-dir "$QC_DIR" \
    --samples "$LIB" \
    --kraken-db "$KRAKEN_DB" --kraken-max-reads "$KRAKEN_MAX_READS" --threads 8
EOF
)
echo "tools: array job $TOOLS ($N_LIBS libraries)"
DEPS="$DEPS:$TOOLS"

TABLE=$(sbatch "${COMMON[@]}" --job-name=atac_table $TABLE_RES --dependency=afterok:"$DEPS" \
        --kill-on-invalid-dep=yes --output="$QC_DIR/logs/table_%j.log" <<'EOF'
#!/bin/bash
set -euo pipefail
source "$CONDA_SH"; conda activate namhub_atac_qc
# no group with >= 2 replicates: no IDR tables, and the IDR columns are blank
IDR_ARGS=()
shopt -s nullglob
IDR_FILES=("$QC_DIR"/idr/*/concordance_metrics.tsv)
[ ${#IDR_FILES[@]} -eq 0 ] || IDR_ARGS=(--idr-metrics "${IDR_FILES[@]}")
python "$SCRIPTS/atacseq_qc_metrics.py" \
    --atac-dir "$ATAC_DIR" --samplesheet "$SAMPLESHEET" --qc-dir "$QC_DIR" \
    --encode-qc "$QC_DIR/atacseq_encode_qc_metrics.tsv" \
    ${IDR_ARGS[@]+"${IDR_ARGS[@]}"} \
    --out "$QC_DIR/atacseq_qc_metrics.tsv"
EOF
)
echo "table: job $TABLE (after $DEPS)"
echo "Output: $QC_DIR"
