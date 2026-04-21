#!/bin/bash -l

#help function
Help()
{
   # Display Help
   echo
   echo
   echo "CVRCseq"
   echo
   echo "   This script exectues scripts/cat_rename.py to concatenate .fastq files from multiple lanes and renames them based on config/samples_info.tab. It then loads the environment. Finally, it launches the appropriate pipeline"
   echo "   For more detail see https://github.com/mgildea87/CVRCseq"
   echo
   echo "   Syntax: scriptTemplate [-h|-s|-w|-d|-c|-i]"
   echo "       options:"
   echo "           -h     help"
   echo "           -d     .fastq directory. location of .fastq files"
   echo "           -s     additional arguments to pass to snakemake (quote multiple flags: -s \"--dryrun --quiet\")"
   echo "           -c     Skip cat_rename.py. Use to skip copying, concatenating, and renaming of .fastq files to local directory." 
   echo "           -i     Singularity .sif image path override. If omitted, default image is used when available"
   echo "           -w     workflow. Can be 1 of:"
   echo "                             'RNAseq_SE' - single end reads, fastqc, fastp, STAR, featurecounts"
   echo "                             'RNAseq_PE' - paired end reads, fastqc, fastp, STAR, featurecounts"
   echo "                             'RNAseq_PE_HISAT2_stringtie' - paired end reads, fastqc, fastp, HISAT2, stringtie"
   echo "                             'RNAseq_PE_HISAT2_stringtie_nvltrx' - paired end reads, fastqc, fastp, HISAT2, stringtie (novel transcript assembly)"
   echo "                             'RNAseqTE_PE' - paired end reads, fastqc, fastp, STAR, TEcount"
   echo "                             'sRNAseq_SE' - single end reads, fastqc, umi-tools, STAR, featurecounts"
   echo "                             'ChIPseq_PE' - paired end reads, fastqc, fastp, bowtie2, macs2"
   echo "                             'CUT-RUN_PE' - paired end reads, fastqc, fastp, bowtie2, seacr"
   echo "                             'ATACseq_PE' - paired end reads, fastqc, fastp, bowtie2, macs2"
   echo
   echo "   If snakemake reports a locked directory (e.g. after a hard crash), load the environment and unlock with:"
   echo "       snakemake --unlock --profile config/profile"
   echo
   echo "   Container defaults:"
   echo "       Default image path: /gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq.sif"
   echo "       Override default with env var: CVRCSEQ_SIF=/path/to/image.sif"
   echo
}

#parse arguments
sif_path_set="no"
while getopts ":w:s:ci:d:h" arg; do
    case $arg in
        w) workflow=$OPTARG;;
        s) snakemake_arg=$OPTARG;;
        c) skip_cat_rename='skip';;
        i) sif_path=$OPTARG; sif_path_set="yes";;
        d) fastq_directory=$OPTARG;;
        h) # display help 
            Help
            exit;; 
    esac
done

#Check required arguments
if [[ -z "$workflow" || -z "$fastq_directory" ]]; then
    echo "Error: -w (workflow) and -d (fastq directory) are required."
    Help
    exit 1
fi

#Check if workflow (-w) exists in available workflows. If not, exit.
workflow_options=( "RNAseq_SE" "sRNAseq_SE" "RNAseq_PE" "RNAseq_PE_HISAT2_stringtie" "RNAseq_PE_HISAT2_stringtie_nvltrx" "ChIPseq_PE" "CUT-RUN_PE" "ATACseq_PE" "RNAseqTE_PE" )

if printf '%s\n' "${workflow_options[@]}" | grep -Fxq -- "$workflow"; then
    echo "$workflow"
else
    echo "Workflow does not exist. Select one from the list in -h"
    exit 1
fi

# Configure execution mode: default to singularity image, fallback to host conda.
default_sif_path="${CVRCSEQ_SIF:-/gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq.sif}"

if [[ "$sif_path_set" = "yes" ]]; then
  selected_sif_path="$sif_path"
else
  selected_sif_path="$default_sif_path"
fi

if [[ -n "$selected_sif_path" && -f "$selected_sif_path" ]]; then
  sif_path="$selected_sif_path"

  # Auto-load singularity module when needed so -i works without manual setup.
  if ! command -v singularity >/dev/null 2>&1; then
    if type module >/dev/null 2>&1; then
      module load singularity/3.11.5 >/dev/null 2>&1 || true
    fi

    if ! command -v singularity >/dev/null 2>&1; then
      echo "Error: 'singularity' not found. Load singularity module or check PATH."
      exit 1
    fi
  fi

  echo
  echo "==================================================="
  echo " Execution mode : CONTAINER (Singularity)"
  echo " Image           : $sif_path"
  echo "==================================================="
  echo
  use_container="yes"
else
  if [[ "$sif_path_set" = "yes" ]]; then
    echo "Error: Singularity image not found: $selected_sif_path"
    exit 1
  fi

  echo
  echo "==================================================="
  echo " Execution mode : HOST (conda environment)"
  echo " Reason         : No Singularity image found at:"
  echo "                  $selected_sif_path"
  echo "==================================================="
  echo
  use_container="no"
fi

# Build snakemake container flags — snakemake runs natively but wraps each
# submitted job inside the container via --use-singularity.
if [[ "$use_container" = "yes" ]]; then
  container_snakemake_args=(--use-singularity --singularity-args "--bind /gpfs")
  sed -i "s|^singularity_image:.*|singularity_image: \"$sif_path\"|" config/config.yaml
else
  container_snakemake_args=()
  sed -i 's|^singularity_image:.*|singularity_image: ""|' config/config.yaml
fi

run_tool() {
  if [[ "$use_container" = "yes" ]]; then
    singularity exec --bind /gpfs "$sif_path" "$@"
  else
    "$@"
  fi
}

# load conda environment on host only
if [[ "$use_container" = "no" ]]; then
  source workflow/scripts/condaload_CVRCseq.sh
fi

# output environment info
run_tool conda list > conda_env.txt

# Write workflow into config so it doesn't need to be set there manually
sed -i "s/^workflow:.*/workflow: \"$workflow\"/" config/config.yaml

mkdir -p "$workflow"/inputs/fastq
mkdir -p slurm_logs

printf "fastq_directory: %s\nworkflow: %s\n" "$fastq_directory" "$workflow" > snakemake_init_commands.txt

skip_cat_rename=${skip_cat_rename:-'dont_skip'}


if [[ $skip_cat_rename = "skip" ]] ; then
  #launch snakemake without running cat_rename.py first
  snakemake $snakemake_arg "${container_snakemake_args[@]}" --profile config/profile --config workflow=$workflow --rerun-incomplete || exit 1
  snakemake --report workflow/snake_make_report.html || exit 1
  run_tool multiqc . --force || exit 1
else
  if ! run_tool python workflow/scripts/cat_rename.py "$fastq_directory" "$workflow"; then
    echo "Exiting..."
    exit 1
  fi
  #launch snakemake
  snakemake $snakemake_arg "${container_snakemake_args[@]}" --profile config/profile --config workflow=$workflow --rerun-incomplete || exit 1
  snakemake --report workflow/snake_make_report.html || exit 1
  run_tool multiqc . --force --interactive || exit 1
fi

