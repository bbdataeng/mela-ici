#!/bin/sh

set -euo pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base directory containing all melanoma datasets
export BASE_DIR=/data/alice/melanoma

# Dataset name passed as first argument (e.g., Hugo-2016)
export DATASET=${1:?Please provide dataset name, e.g. Hugo-2016}

# Path to the dataset folder
export DATASET_DIR=$BASE_DIR/${DATASET}

# Subdirectory containing FASTQs
export FASTQ_DIR=$DATASET_DIR/fastq

# Output and log directories
export LOG_DIR=$BASE_DIR/logs
export RNASEQ_DIR=$DATASET_DIR/rnaseq
export HLATYPING_DIR=$DATASET_DIR/hlatyping

# Reference data
export REF_DIR=$BASE_DIR/reference


# SET WHICH COMPONENTS TO RUN (0 = skip, 1 = run)
export RUN_RNASEQ_SAMPLESHEET=1
export RUN_RNASEQ=0



# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

if [ ! -d $LOG_DIR ]; then mkdir -p $LOG_DIR;fi



# ==============================================================================
# FUNCTION: Convert seconds to HH:MM:SS
# ==============================================================================
format_time() {
  local total_seconds=$1
  local hours=$((total_seconds / 3600))
  local minutes=$(((total_seconds % 3600) / 60))
  local seconds=$((total_seconds % 60))
  printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds"
}
export -f format_time



# ==============================================================================
# FUNCTION: Create samplesheet for nf-core/rnaseq
# ==============================================================================

RNASEQ_SAMPLESHEET() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting nf-core/rnaseq SAMPLESHEET for $DATASET..."

  mkdir -p "$RNASEQ_DIR"

  local start_time=$(date +%s)
  local samplesheet="$RNASEQ_DIR/samplesheet.csv"

  echo "sample,fastq_1,fastq_2,strandedness" > "$samplesheet"
  for base in $(ls "$FASTQ_DIR"/*_1.fastq.gz | sed 's/_1.fastq.gz//'); do
    sample=$(basename "$base")
    echo "$sample,${base}_1.fastq.gz,${base}_2.fastq.gz,auto" >> "$samplesheet"
  done

  local end_time=$(date +%s)
  local duration=$((end_time - start_time))

  echo "✅ [$(date '+%Y-%m-%d %H:%M:%S')] Completed nf-core/rnaseq SAMPLESHEET for $DATASET in $(format_time $duration)"
}
export -f RNASEQ_SAMPLESHEET


echo "sample,fastq_1,fastq_2,strandedness" > /data/alice/melanoma/Auslander-2018/samplesheet.csv
for base in $(ls /data/alice/melanoma/Auslander-2018/fastq/*_1.fastq.gz | sed  's/_1.fastq.gz//');  do 
    sample=$(basename $base);
    echo -e "$sample,${base}_1.fastq.gz,${base}_2.fastq.gz,auto"; 
done >> /data/alice/melanoma/Auslander-2018/samplesheet.csv



# ==============================================================================
# FUNCTION: Run nf-core/rnaseq
# ==============================================================================

RNASEQ() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting nf-core/rnaseq for $DATASET..."
  
  local start_time=$(date +%s)

  nextflow 24.10.5 run nf-core/rnaseq \
    --input "$RNASEQ_DIR/samplesheet.csv" \
    --outdir "$RNASEQ_DIR/nextflow" \
    --gtf "$REF_DIR/gencode.v47.primary_assembly.basic.annotation.gtf.gz" \
    --fasta "$REF_DIR/GRCh38.primary_assembly.genome.fa.gz" \
    -r 3.18.0 \
    -profile docker 2>&1 | tee "$RNASEQ_DIR/output_rnaseq.log"

  local end_time=$(date +%s)
  local duration=$((end_time - start_time))

  echo "✅ [$(date '+%Y-%m-%d %H:%M:%S')] Completed nf-core/rnaseq for $DATASET in $(format_time $duration)"
}
export -f RNASEQ



# ==============================================================================
# FUNCTION: Log header information
# ==============================================================================

LOG_HEADER() {
  cat <<EOF
################################################################################
# Dataset: $DATASET
# Run samplesheet flag: $RUN_RNASEQ_SAMPLESHEET
# Run nf-core/rnaseq flag: $RUN_RNASEQ
################################################################################
EOF
}
export -f LOG_HEADER



# ==============================================================================
# MAIN PIPELINE FUNCTION
# ==============================================================================

PIPELINE() {
  local start_time=$(date +%s)
  echo "================== PIPELINE START for $DATASET =================="
  echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
  
  local sample_log_dir=$LOG_DIR/${DATASET}
  mkdir -p "$sample_log_dir"
  exec > >(tee -i "$sample_log_dir/${DATASET}_pipeline.log") 2>&1

  LOG_HEADER

  [ $RUN_RNASEQ_SAMPLESHEET -eq 1 ] && RNASEQ_SAMPLESHEET
  [ $RUN_RNASEQ -eq 1 ] && RNASEQ

  local end_time=$(date +%s)
  local duration=$((end_time - start_time))

  echo "================== PIPELINE END for $DATASET =================="
  echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "Total duration: $(format_time $duration)"
}
export -f PIPELINE



# ==============================================================================
# RUN
# ==============================================================================

PIPELINE