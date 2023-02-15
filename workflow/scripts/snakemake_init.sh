#!/bin/bash -l

#help function
Help()
{
   # Display Help
   echo
   echo
   echo "CVRCseq"
   echo
   echo "   This script exectues scripts/cat_rename.py to concatenate .fastq files from multiple lanes and renames them based on config/samples_info.tab. It then launches the appropriate pipeline"
   echo "   For more detail see github/...."
   echo
   echo "   Syntax: scriptTemplate [-h|-c|-w|-d]"
   echo "       options:"
   echo "           -h     help"
   echo "           -d     .fastq directory"
   echo "           -c     Absolute path to conda environment"
   echo "           -w     workflow. Can be 1 of:"
   echo "                             'RNAseq_SE' - single end reads, fastqc, fastp, STAR, featurecounts"
   echo "                             'RNAseq_PE' - paired end reads, fastqc, fastp, STAR, featurecounts"
   echo "                             'RNAseq_PE_HISAT2_stringtie' - paired end reads, fastqc, fastp, HISAT2, stringtie"
   echo "                             'RNAseq_PE_HISAT2_stringtie_nvltrx' - paired end reads, fastqc, fastp, HISAT2, stringtie (novel transcript assembly)"
   echo "                             'ChIPseq' - paired end reads, fastqc, fastp, bowtie2, macs2"
   echo "                             'CUT-RUN' - paired end reads, fastqc, fastp, bowtie2, seacr"
   echo
}

#parse arguments
while getopts ":c:w:d:h" arg; do
    case $arg in
        c) conda_env=$OPTARG;;
        w) workflow=$OPTARG;;
        d) fastq_directory=$OPTARG;;
        h) # display help 
            Help
            exit;; 
    esac
done

#Check if workflow (-w) exists in available workflows. If not, exit.
workflow_options=( "RNAseq_SE" "RNAseq_PE" "RNAseq_PE_HISAT2_stringtie" "RNAseq_PE_HISAT2_stringtie_nvltrx" "ChIPseq" "CUT-RUN" )

if printf '%s\n' "${workflow_options[@]}" | grep -Fxq -- $workflow; then
    echo $workflow
else
    echo "Workflow does not exist. Select one from the list in -h"
    exit
fi

module load miniconda3/cpu/4.9.2
conda activate $conda_env
mkdir -p "$workflow"/inputs/fastq
mkdir slurm_logs

#exit if this script errors
if ! python workflow/scripts/cat_rename.py $fastq_directory $workflow ${1}; then
    echo "Exiting..."
    exit
fi

#launch snakemake
snakemake --profile config/profile --config workflow=$workflow
snakemake --report workflow/snake_make_report.html
multiqc .