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
   echo "   Syntax: scriptTemplate [-h|-s|-w|-d]"
   echo "       options:"
   echo "           -h     help"
   echo "           -d     .fastq directory"
   echo "           -s     arguments to pass to snakemake"
   echo "           -w     workflow. Can be 1 of:"
   echo "                             'RNAseq_SE' - single end reads, fastqc, fastp, STAR, featurecounts"
   echo "                             'RNAseq_PE' - paired end reads, fastqc, fastp, STAR, featurecounts"
   echo "                             'RNAseq_PE_HISAT2_stringtie' - paired end reads, fastqc, fastp, HISAT2, stringtie"
   echo "                             'RNAseq_PE_HISAT2_stringtie_nvltrx' - paired end reads, fastqc, fastp, HISAT2, stringtie (novel transcript assembly)"
   echo "                             'ChIPseq_PE' - paired end reads, fastqc, fastp, bowtie2, macs2"
   echo "                             'CUT-RUN_PE' - paired end reads, fastqc, fastp, bowtie2, seacr"
   echo "                             'ATACseq_PE' - paired end reads, fastqc, fastp, bowtie2, macs2"  
   echo
}

#parse arguments
while getopts ":w:s:d:h" arg; do
    case $arg in
        w) workflow=$OPTARG;;
        s) snakemake_arg=$OPTARG;;
        d) fastq_directory=$OPTARG;;
        h) # display help 
            Help
            exit;; 
    esac
done

#Check if workflow (-w) exists in available workflows. If not, exit.
workflow_options=( "RNAseq_SE" "RNAseq_PE" "RNAseq_PE_HISAT2_stringtie" "RNAseq_PE_HISAT2_stringtie_nvltrx" "ChIPseq_PE" "CUT-RUN_PE" "ATACseq_PE" )

if printf '%s\n' "${workflow_options[@]}" | grep -Fxq -- $workflow; then
    echo $workflow
else
    echo "Workflow does not exist. Select one from the list in -h"
    exit
fi

# load conda environment
source workflow/scripts/condaload_CVRCseq.sh
mkdir -p "$workflow"/inputs/fastq
mkdir slurm_logs

printf "fastq_directory: %s\nworkflow: %s\n" "$fastq_directory" "$workflow" > snakemake_init_commands.txt

#exit if this script errors
if ! python workflow/scripts/cat_rename.py $fastq_directory $workflow ${1}; then
    echo "Exiting..."
    exit
fi

#launch snakemake
snakemake "$snakemake_arg" --profile config/profile --config workflow=$workflow
snakemake --report workflow/snake_make_report.html
multiqc . --force