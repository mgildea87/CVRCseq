#!/usr/bin/env bash

# Subsample paired FASTQ files to a fixed number of reads and write
# them with generic names suitable for integration testing.
# Uses only zcat/head/gzip — no extra tools required.
# Pairing is preserved because the same number of lines is
# taken from both R1 and R2.

set -euo pipefail

FASTQ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/fastq" && pwd)"
READS=50000
LINES=$(( READS * 4 ))

# original_R1  original_R2  generic_basename
declare -a PAIRS=(
  "293-1_S4_L002_R1_001.fastq.gz  293-1_S4_L002_R2_001.fastq.gz  test_ATAC"
  "FemaleBL1_R1_002.fastq.gz      FemaleBL1_R2_002.fastq.gz      test_RNAseq"
  "m1_-H3K18La_pAb_S3_R1_001.fastq.gz  m1_-H3K18La_pAb_S3_R2_001.fastq.gz  test_Antibody"
  "m1_-_IgG_S2_R1_001.fastq.gz   m1_-_IgG_S2_R2_001.fastq.gz   test_Control"
)

# Rename any previously generated _sub files to the generic names
declare -A SUB_TO_GENERIC=(
  ["293-1_S4_L002_R1_001_sub.fastq.gz"]="test_ATAC_R1.fastq.gz"
  ["293-1_S4_L002_R2_001_sub.fastq.gz"]="test_ATAC_R2.fastq.gz"
  ["FemaleBL1_R1_002_sub.fastq.gz"]="test_RNAseq_R1.fastq.gz"
  ["FemaleBL1_R2_002_sub.fastq.gz"]="test_RNAseq_R2.fastq.gz"
  ["m1_-H3K18La_pAb_S3_R1_001_sub.fastq.gz"]="test_Antibody_R1.fastq.gz"
  ["m1_-H3K18La_pAb_S3_R2_001_sub.fastq.gz"]="test_Antibody_R2.fastq.gz"
  ["m1_-_IgG_S2_R1_001_sub.fastq.gz"]="test_Control_R1.fastq.gz"
  ["m1_-_IgG_S2_R2_001_sub.fastq.gz"]="test_Control_R2.fastq.gz"
)

for old in "${!SUB_TO_GENERIC[@]}"; do
  old_path="${FASTQ_DIR}/${old}"
  new_path="${FASTQ_DIR}/${SUB_TO_GENERIC[$old]}"
  if [[ -f "${old_path}" && ! -f "${new_path}" ]]; then
    echo "Renaming ${old} -> ${SUB_TO_GENERIC[$old]}"
    mv "${old_path}" "${new_path}"
  fi
done

for pair in "${PAIRS[@]}"; do
  read -r r1 r2 base <<< "${pair}"
  r1_out="${FASTQ_DIR}/${base}_R1.fastq.gz"
  r2_out="${FASTQ_DIR}/${base}_R2.fastq.gz"

  if [[ -f "${r1_out}" && -f "${r2_out}" ]]; then
    echo "Skipping ${base} (already exists)"
    continue
  fi

  r1_path="${FASTQ_DIR}/${r1}"
  r2_path="${FASTQ_DIR}/${r2}"

  echo "Subsampling ${r1} -> ${base}_R1.fastq.gz ..."
  { zcat "${r1_path}" | head -n "${LINES}" || true; } | gzip > "${r1_out}"

  echo "Subsampling ${r2} -> ${base}_R2.fastq.gz ..."
  { zcat "${r2_path}" | head -n "${LINES}" || true; } | gzip > "${r2_out}"

  echo "  Done."
done

echo ""
echo "Subsampled FASTQ files written to ${FASTQ_DIR}"
