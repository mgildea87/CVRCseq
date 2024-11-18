#!/bin/bash
#SBATCH -J ATACseq_launch
#SBATCH --mem=10000
#SBATCH --cpus-per-task=1
#SBATCH -p fn_medium
#SBATCH --export=ALL
#SBATCH --time=48:00:00

 bash workflow/scripts/snakemake_init.sh -d ../fastq/ -w ATACseq_PE
