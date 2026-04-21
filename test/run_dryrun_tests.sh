#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="${ROOT_DIR}/.test-work"
SNAKEMAKE_CMD="${SNAKEMAKE_CMD:-snakemake}"
KEEP=false
INTEGRATION=false
WORKFLOW_FILTER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --keep)
      KEEP=true
      shift
      ;;
    --integration)
      INTEGRATION=true
      shift
      ;;
    --workflow=*)
      WORKFLOW_FILTER="${1#--workflow=}"
      shift
      ;;
    --workflow|-w)
      if [[ $# -lt 2 ]]; then
        echo "Error: $1 requires a workflow name"
        echo "Usage: $0 [--integration] [--keep] [--workflow NAME]"
        exit 1
      fi
      WORKFLOW_FILTER="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      echo "Usage: $0 [--integration] [--keep] [--workflow NAME]"
      exit 1
      ;;
  esac
done

declare -a CONFIGS=(
  "test/config/config_ATAC.yaml"
  "test/config/config_CUT-RUN.yaml"
  "test/config/config_ChIPseq.yaml"
  "test/config/config_RNAseq_PE.yaml"
  "test/config/config_RNAseq_SE.yaml"
  "test/config/config_RNAseqTE_PE.yaml"
  "test/config/config_RNAseq_PE_HISAT2_stringtie.yaml"
  "test/config/config_RNAseq_PE_HISAT2_stringtie_nvltrx.yaml"
  "test/config/config_sRNAseq_SE.yaml"
)

mkdir -p "${WORK_DIR}"

if ! command -v "${SNAKEMAKE_CMD%% *}" >/dev/null 2>&1; then
  echo "Error: could not find Snakemake command: ${SNAKEMAKE_CMD}"
  echo "Set SNAKEMAKE_CMD if needed, for example:"
  echo "  SNAKEMAKE_CMD='conda run -n CVRCseq snakemake' bash test/run_dryrun_tests.sh"
  exit 1
fi

if ! ${SNAKEMAKE_CMD} --version >/dev/null 2>&1; then
  echo "Error: Snakemake command is present but not runnable: ${SNAKEMAKE_CMD}"
  echo "Activate the proper environment, or set SNAKEMAKE_CMD explicitly."
  echo "Example:"
  echo "  SNAKEMAKE_CMD='conda run -n CVRCseq snakemake' bash test/run_dryrun_tests.sh"
  exit 1
fi

if [[ "${INTEGRATION}" == true ]]; then
  echo "Running Snakemake INTEGRATION tests (jobs will be submitted via Slurm)"
else
  echo "Running Snakemake dry-run tests"
fi
echo "Root: ${ROOT_DIR}"
echo "Work: ${WORK_DIR}"

declare -a FAILED_WORKFLOWS=()

prepare_rule_directories() {
  local run_dir="$1"
  local wf_name="$2"

  case "${wf_name}" in
    "ATACseq_PE")
      mkdir -p "${run_dir}/${wf_name}/results/peaks/MACS2/qc"
      mkdir -p "${run_dir}/${wf_name}/results/alignment/frag_len"
      mkdir -p "${run_dir}/${wf_name}/results/alignment/idxstat"
      mkdir -p "${run_dir}/${wf_name}/results/fastqc_post_trim"
      mkdir -p "${run_dir}/${wf_name}/results/logs/trim_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/alignment_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/MACS2"
      ;;
    "CUT-RUN_PE")
      mkdir -p "${run_dir}/${wf_name}/results/peaks/seacr/qc"
      mkdir -p "${run_dir}/${wf_name}/results/peaks/MACS2/qc"
      mkdir -p "${run_dir}/${wf_name}/results/alignment/bed"
      mkdir -p "${run_dir}/${wf_name}/results/alignment/frag_len"
      mkdir -p "${run_dir}/${wf_name}/results/fastqc_post_trim"
      mkdir -p "${run_dir}/${wf_name}/results/logs/trim_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/alignment_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/MACS2"
      ;;
    "ChIPseq_PE")
      mkdir -p "${run_dir}/${wf_name}/results/peaks/MACS2/qc"
      mkdir -p "${run_dir}/${wf_name}/results/alignment/frag_len"
      mkdir -p "${run_dir}/${wf_name}/results/fastqc_post_trim"
      mkdir -p "${run_dir}/${wf_name}/results/logs/trim_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/alignment_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/MACS2"
      ;;
    "RNAseq_PE")
      mkdir -p "${run_dir}/${wf_name}/results/feature_counts"
      mkdir -p "${run_dir}/${wf_name}/results/fastqc_post_trim"
      mkdir -p "${run_dir}/${wf_name}/results/logs/trim_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/alignment_reports"
      ;;
    "RNAseq_SE")
      mkdir -p "${run_dir}/${wf_name}/results/feature_counts"
      mkdir -p "${run_dir}/${wf_name}/results/fastqc_post_trim"
      mkdir -p "${run_dir}/${wf_name}/results/logs/trim_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/alignment_reports"
      ;;
    "RNAseqTE_PE")
      mkdir -p "${run_dir}/${wf_name}/results/TEcount"
      mkdir -p "${run_dir}/${wf_name}/results/fastqc_post_trim"
      mkdir -p "${run_dir}/${wf_name}/results/logs/trim_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/alignment_reports"
      ;;
    "RNAseq_PE_HISAT2_stringtie")
      mkdir -p "${run_dir}/${wf_name}/results/stringtie"
      mkdir -p "${run_dir}/${wf_name}/results/logs/trim_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/alignment_reports"
      ;;
    "RNAseq_PE_HISAT2_stringtie_nvltrx")
      mkdir -p "${run_dir}/${wf_name}/results/stringtie/merged"
      mkdir -p "${run_dir}/${wf_name}/results/logs/trim_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/alignment_reports"
      ;;
    "sRNAseq_SE")
      mkdir -p "${run_dir}/${wf_name}/results/feature_counts"
      mkdir -p "${run_dir}/${wf_name}/results/umi_tools_trim"
      mkdir -p "${run_dir}/${wf_name}/results/fastqc_post_trim"
      mkdir -p "${run_dir}/${wf_name}/results/logs/umi_tools_trim_reports"
      mkdir -p "${run_dir}/${wf_name}/results/logs/alignment_reports"
      ;;
  esac
}

prepare_fastq_links() {
  local run_dir="$1"
  local wf_name="$2"
  local sample_file_abs="$3"
  local wf_fastq_dir="${run_dir}/${wf_name}/inputs/fastq"
  local src_fastq_dir="${ROOT_DIR}/test/fastq"

  mkdir -p "${wf_fastq_dir}"

  for src in "${src_fastq_dir}"/*.fastq.gz; do
    ln -sfn "${src}" "${wf_fastq_dir}/$(basename "${src}")"
  done

  resolve_fastq_source() {
    local filename="$1"
    local exact="${src_fastq_dir}/${filename}"
    local alias
    local lane_pattern
    local lane_match

    if [[ -f "${exact}" ]]; then
      echo "${exact}"
      return 0
    fi

    # Backward-compatible aliases for legacy sample tables.
    case "${filename}" in
      293-1_S4_R1_001.fastq.gz|293-1_S4_L002_R1_001.fastq.gz|293-1_S4_L002_R1_001_sub.fastq.gz)
        alias="test_ATAC_R1.fastq.gz"
        ;;
      293-1_S4_R2_001.fastq.gz|293-1_S4_L002_R2_001.fastq.gz|293-1_S4_L002_R2_001_sub.fastq.gz)
        alias="test_ATAC_R2.fastq.gz"
        ;;
      FemaleBL1_R1_002.fastq.gz|FemaleBL1_R1_002_sub.fastq.gz)
        alias="test_RNAseq_R1.fastq.gz"
        ;;
      FemaleBL1_R2_002.fastq.gz|FemaleBL1_R2_002_sub.fastq.gz)
        alias="test_RNAseq_R2.fastq.gz"
        ;;
      m1_-H3K18La_pAb_S3_R1_001.fastq.gz|m1_-H3K18La_pAb_S3_R1_001_sub.fastq.gz)
        alias="test_Antibody_R1.fastq.gz"
        ;;
      m1_-H3K18La_pAb_S3_R2_001.fastq.gz|m1_-H3K18La_pAb_S3_R2_001_sub.fastq.gz)
        alias="test_Antibody_R2.fastq.gz"
        ;;
      m1_-_IgG_S2_R1_001.fastq.gz|m1_-_IgG_S2_R1_001_sub.fastq.gz)
        alias="test_Control_R1.fastq.gz"
        ;;
      m1_-_IgG_S2_R2_001.fastq.gz|m1_-_IgG_S2_R2_001_sub.fastq.gz)
        alias="test_Control_R2.fastq.gz"
        ;;
      *)
        alias=""
        ;;
    esac

    if [[ -n "${alias}" && -f "${src_fastq_dir}/${alias}" ]]; then
      echo "${src_fastq_dir}/${alias}"
      return 0
    fi

    lane_pattern="$(echo "${filename}" | sed -E 's/_R([12])_/_L*_R\1_/')"
    lane_match="$(compgen -G "${src_fastq_dir}/${lane_pattern}" | head -n 1 || true)"

    if [[ -n "${lane_match}" ]]; then
      echo "${lane_match}"
      return 0
    fi

    return 1
  }

  while IFS=$'\t' read -r file_r1 file_r2 sample condition replicate antibody sample_name _ || [[ -n "${file_r1:-}" ]]; do
    if [[ "${file_r1}" == "File_Name_R1" ]]; then
      continue
    fi

    if [[ "${wf_name}" == "ChIPseq_PE" || "${wf_name}" == "CUT-RUN_PE" ]]; then
      sample_alias="${sample}_${condition}_${replicate}_${antibody}"
    else
      sample_alias="${sample}_${condition}_${replicate}"
    fi

    if ! src_r1="$(resolve_fastq_source "${file_r1}")"; then
      echo "Error: could not resolve FASTQ source for ${file_r1}" >&2
      exit 1
    fi
    ln -sfn "${src_r1}" "${wf_fastq_dir}/${sample_alias}_R1.fastq.gz"

    if [[ -n "${file_r2:-}" ]]; then
      if ! src_r2="$(resolve_fastq_source "${file_r2}")"; then
        echo "Error: could not resolve FASTQ source for ${file_r2}" >&2
        exit 1
      fi
      ln -sfn "${src_r2}" "${wf_fastq_dir}/${sample_alias}_R2.fastq.gz"
    fi
  done < "${sample_file_abs}"
}

ACTIVE_CONFIGS=("${CONFIGS[@]}")

if [[ -n "${WORKFLOW_FILTER}" ]]; then
  _found=false
  for _c in "${ACTIVE_CONFIGS[@]}"; do
    _wf="$(awk -F '"' '/^workflow:/ {print $2}' "${ROOT_DIR}/${_c}")"
    [[ "${_wf}" == "${WORKFLOW_FILTER}" ]] && _found=true && break
  done
  if [[ "${_found}" == false ]]; then
    echo "Error: workflow '${WORKFLOW_FILTER}' not found."
    echo "Available workflows:"
    for _c in "${ACTIVE_CONFIGS[@]}"; do
      awk -F '"' '/^workflow:/ {print "  - " $2}' "${ROOT_DIR}/${_c}"
    done
    exit 1
  fi
fi

for cfg in "${ACTIVE_CONFIGS[@]}"; do
  cfg_path="${ROOT_DIR}/${cfg}"
  wf_name="$(awk -F '"' '/^workflow:/ {print $2}' "${cfg_path}")"

  if [[ -n "${WORKFLOW_FILTER}" && "${wf_name}" != "${WORKFLOW_FILTER}" ]]; then
    continue
  fi

  sample_file_rel="$(awk -F '"' '/^sample_file:/ {print $2}' "${cfg_path}")"
  sample_file_abs="${ROOT_DIR}/${sample_file_rel}"
  run_dir="${WORK_DIR}/${wf_name}"

  rm -rf "${run_dir}"
  mkdir -p "${run_dir}/inputs"
  mkdir -p "${run_dir}/${wf_name}"
  prepare_rule_directories "${run_dir}" "${wf_name}"
  prepare_fastq_links "${run_dir}" "${wf_name}" "${sample_file_abs}"
  ln -sfn "${run_dir}/${wf_name}/inputs/fastq" "${run_dir}/inputs/fastq"

  echo ""
  if [[ "${INTEGRATION}" == true ]]; then
    echo "--- [${wf_name}] integration run with ${cfg} ---"
    mkdir -p "${run_dir}/slurm_logs"
    if ${SNAKEMAKE_CMD} \
      --snakefile "${ROOT_DIR}/workflow/Snakefile" \
      --configfile "${cfg_path}" \
      --config "sample_file=${sample_file_abs}" \
      --directory "${run_dir}" \
      --profile "${ROOT_DIR}/config/profile" \
      --rerun-incomplete; then
      echo "Pipeline complete. Checking outputs..."
      if bash "${ROOT_DIR}/test/check_outputs.sh" "${wf_name}" "${run_dir}"; then
        echo "PASS: ${wf_name}"
      else
        echo "FAIL: ${wf_name} (output validation failed)"
        FAILED_WORKFLOWS+=("${wf_name}")
      fi
    else
      echo "FAIL: ${wf_name} (pipeline error)"
      FAILED_WORKFLOWS+=("${wf_name}")
    fi
  else
    echo "--- [${wf_name}] dry-run with ${cfg} ---"
    if ${SNAKEMAKE_CMD} \
      --snakefile "${ROOT_DIR}/workflow/Snakefile" \
      --configfile "${cfg_path}" \
      --config "sample_file=${sample_file_abs}" \
      --directory "${run_dir}" \
      --cores 1 \
      --dry-run; then
      echo "PASS: ${wf_name}"
    else
      echo "FAIL: ${wf_name}"
      FAILED_WORKFLOWS+=("${wf_name}")
    fi
  fi
done

echo ""
if [[ ${#FAILED_WORKFLOWS[@]} -eq 0 ]]; then
  echo "All configured dry-run tests completed successfully."
  if [[ "${KEEP}" == false ]]; then
    echo "Cleaning up ${WORK_DIR}"
    rm -rf "${WORK_DIR}"
  fi
else
  echo "Dry-run failures (${#FAILED_WORKFLOWS[@]}): ${FAILED_WORKFLOWS[*]}"
  exit 1
fi