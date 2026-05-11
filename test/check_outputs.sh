#!/usr/bin/env bash

# Validate integration test outputs for a given workflow run.
# Called automatically by run_dryrun_tests.sh after each integration run.
#
# Usage: check_outputs.sh <workflow_name> <run_dir>

set -euo pipefail

WF_NAME="$1"
RUN_DIR="$2"
RESULTS="${RUN_DIR}/${WF_NAME}/results"

PASS=0
FAIL=0

check_exists_nonempty() {
  local label="$1"
  local file="$2"
  if [[ -f "${file}" && -s "${file}" ]]; then
    echo "  PASS: ${label}"
    PASS=$((PASS + 1))
  else
    echo "  FAIL: ${label} — missing or empty: ${file}"
    FAIL=$((FAIL + 1))
  fi
}

check_exists() {
  local label="$1"
  local file="$2"
  if [[ -f "${file}" ]]; then
    echo "  PASS: ${label}"
    PASS=$((PASS + 1))
  else
    echo "  FAIL: ${label} — file missing: ${file}"
    FAIL=$((FAIL + 1))
  fi
}

check_min_lines() {
  local label="$1"
  local file="$2"
  local min="$3"
  if [[ ! -f "${file}" ]]; then
    echo "  FAIL: ${label} — file missing: ${file}"
    FAIL=$((FAIL + 1))
    return
  fi
  local lines
  lines=$(wc -l < "${file}")
  if [[ "${lines}" -ge "${min}" ]]; then
    echo "  PASS: ${label} (${lines} lines)"
    PASS=$((PASS + 1))
  else
    echo "  FAIL: ${label} — only ${lines} lines (expected >= ${min}): ${file}"
    FAIL=$((FAIL + 1))
  fi
}

check_bam() {
  local label="$1"
  local file="$2"
  if [[ ! -f "${file}" ]]; then
    echo "  FAIL: ${label} — file missing: ${file}"
    FAIL=$((FAIL + 1))
    return
  fi
  if samtools quickcheck "${file}" 2>/dev/null; then
    echo "  PASS: ${label} (valid BAM)"
    PASS=$((PASS + 1))
  else
    echo "  FAIL: ${label} — samtools quickcheck failed: ${file}"
    FAIL=$((FAIL + 1))
  fi
}

echo "Checking outputs for ${WF_NAME} in ${RUN_DIR}"

case "${WF_NAME}" in

  ATACseq_PE)
    for bam in "${RESULTS}"/alignment/*_filtered_sorted.bam; do
      check_bam "filtered sorted BAM: $(basename "${bam}")" "${bam}"
    done
    for bw in "${RESULTS}"/alignment/*_dedup_filtered_sorted.bw; do
      check_exists_nonempty "bigwig: $(basename "${bw}")" "${bw}"
    done
    for peak in "${RESULTS}"/peaks/MACS2/*_peaks.narrowPeak; do
      check_min_lines "narrowPeak: $(basename "${peak}")" "${peak}" 1
    done
    check_min_lines "FRiP summary" "${RESULTS}/peaks/MACS2/qc/frip_summary_detailed.tsv" 2
    ;;

  ChIPseq_PE)
    for bam in "${RESULTS}"/alignment/*_sorted.bam; do
      check_bam "sorted BAM: $(basename "${bam}")" "${bam}"
    done
    for bw in "${RESULTS}"/alignment/*.bw; do
      check_exists_nonempty "bigwig: $(basename "${bw}")" "${bw}"
    done
    for peak in "${RESULTS}"/peaks/MACS2/*_peaks.narrowPeak; do
      check_min_lines "narrowPeak: $(basename "${peak}")" "${peak}" 1
    done
    check_min_lines "FRiP summary" "${RESULTS}/peaks/MACS2/qc/frip_summary_detailed.tsv" 2
    ;;

  CUT-RUN_PE)
    for bam in "${RESULTS}"/alignment/*_sorted.bam; do
      check_bam "sorted BAM: $(basename "${bam}")" "${bam}"
    done
    for peak in "${RESULTS}"/peaks/MACS2/*_peaks.broadPeak; do
      check_min_lines "broadPeak: $(basename "${peak}")" "${peak}" 1
    done
    check_min_lines "MACS2 FRiP summary" "${RESULTS}/peaks/MACS2/qc/frip_summary_detailed.tsv" 2
    ;;

  RNAseq_PE|RNAseq_SE|sRNAseq_SE)
    check_min_lines "count table" "${RESULTS}/feature_counts/count_table.txt" 2
    for bam in "${RESULTS}"/alignment/*.bam; do
      check_bam "alignment BAM: $(basename "${bam}")" "${bam}"
    done
    ;;

  RNAseqTE_PE)
    check_min_lines "TE count table" "${RESULTS}/TEcount/count_table_all.csv" 2
    for bam in "${RESULTS}"/alignment/*.bam; do
      check_bam "alignment BAM: $(basename "${bam}")" "${bam}"
    done
    ;;

  RNAseq_PE_HISAT2_stringtie)
    check_min_lines "gene count matrix" "${RESULTS}/stringtie/gene_count_matrix.csv" 2
    check_min_lines "transcript count matrix" "${RESULTS}/stringtie/transcript_count_matrix.csv" 2
    for bam in "${RESULTS}"/alignment/*.bam; do
      check_bam "alignment BAM: $(basename "${bam}")" "${bam}"
    done
    ;;

  RNAseq_PE_HISAT2_stringtie_nvltrx)
    check_min_lines "gene count matrix" "${RESULTS}/stringtie/merged/gene_count_matrix.csv" 2
    check_min_lines "transcript count matrix" "${RESULTS}/stringtie/merged/transcript_count_matrix.csv" 2
    for bam in "${RESULTS}"/alignment/*.bam; do
      check_bam "alignment BAM: $(basename "${bam}")" "${bam}"
    done
    ;;

  *)
    echo "  WARNING: no output checks defined for workflow '${WF_NAME}'"
    ;;
esac

echo ""
echo "Output checks: ${PASS} passed, ${FAIL} failed"
[[ ${FAIL} -eq 0 ]]
