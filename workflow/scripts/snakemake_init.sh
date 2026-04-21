#!/bin/bash -l

#help function
Help()
{
   # Display Help
   echo
   echo
   echo "CVRCseq"
   echo
   echo "   This script exectues scripts/cat_rename.py to concatenate .fastq files from multiple lanes and renames them based on config/samples_info.tab. It then loads the conda environment via condaload_CVRCseq.sh. Finally, it launches the appropriate pipeline"
   echo "   For more detail see https://github.com/mgildea87/CVRCseq"
   echo
   echo "   Syntax: scriptTemplate [-h|-s|-w|-d|-c|-i]"
   echo "       options:"
   echo "           -h     help"
   echo "           -d     .fastq directory. location of .fastq files"
   echo "           -s     additional arguments to pass to snakemake (quote multiple flags: -s \"--dryrun --quiet\")"
   echo "           -c     Skip cat_rename.py. Use to skip copying, concatenating, and renaming of .fastq files to local directory.
           -i     Path to Singularity .sif image. When provided, all pipeline commands run inside the container.
                  Build the image once with: singularity build --fakeroot CVRCseq.sif CVRCseq.def
                  Example: -i /gpfs/data/cvrcbioinfolab/gildem01/CVRCseq.sif" 
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
   echo "   If snakemake reports a locked directory (e.g. after a hard crash), unlock with:"
   echo "       snakemake --unlock --profile config/profile"
   echo "       singularity exec --bind /gpfs CVRCseq.sif snakemake --unlock --profile config/profile  # if using -i"
   echo
}

#parse arguments
while getopts ":w:s:cd:i:h" arg; do
    case $arg in
        w) workflow=$OPTARG;;
        s) snakemake_arg=$OPTARG;;
        c) skip_cat_rename='skip';;
        d) fastq_directory=$OPTARG;;
        i) sif_path=$OPTARG;;
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

# Set up execution environment — singularity or conda
if [[ -n "$sif_path" ]]; then
    if [[ ! -f "$sif_path" ]]; then
        echo "Error: Singularity image not found: $sif_path"
        exit 1
    fi
    SIF_EXEC="singularity exec --bind /gpfs $sif_path"
    echo "Using Singularity image: $sif_path"
    $SIF_EXEC conda list > conda_env.txt
else
    SIF_EXEC=""
    source workflow/scripts/condaload_CVRCseq.sh
    conda list > conda_env.txt
fi

# Write workflow into config so it doesn't need to be set there manually
sed -i "s/^workflow:.*/workflow: \"$workflow\"/" config/config.yaml

mkdir -p "$workflow"/inputs/fastq
mkdir -p slurm_logs

printf "fastq_directory: %s\nworkflow: %s\n" "$fastq_directory" "$workflow" > snakemake_init_commands.txt

skip_cat_rename=${skip_cat_rename:-'dont_skip'}


if [[ $skip_cat_rename = "skip" ]] ; then
  #launch snakemake without running cat_rename.py first
  $SIF_EXEC snakemake $snakemake_arg --profile config/profile --config workflow=$workflow --rerun-incomplete || exit 1
  $SIF_EXEC snakemake --report workflow/snake_make_report.html
  $SIF_EXEC multiqc . --force
else
  if ! $SIF_EXEC python workflow/scripts/cat_rename.py "$fastq_directory" "$workflow"; then
    echo "Exiting..."
    exit 1
  fi
  #launch snakemake
  $SIF_EXEC snakemake $snakemake_arg --profile config/profile --config workflow=$workflow --rerun-incomplete || exit 1
  $SIF_EXEC snakemake --report workflow/snake_make_report.html
  $SIF_EXEC multiqc . --force --interactive
fi

