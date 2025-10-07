# CVRCseq 
This is a collection of commonly used pipelines integrated into a single workflow via snakemake. Previously, I had all of these as individual snakemake workflows. This workflow is designed to run on NYU's UltraViolet HPC, which utilizes Slurm and has a variety of different node types.

## RNA-seq
There are currently 4 RNA-seq analysis pipelines available

	1.	RNAseq_PE
			paired-end data fastqc > fastp > STAR > featurecounts
	2.	RNAseq_SE
			single-end data fastqc > fastp > STAR > featurecounts
	3.	RNAseq_HISAT2_stringtie
			paired-end data fastqc > fastp > HISAT2 > stringtie
	4.	RNAseq_HISAT2_stringtie_nvltrx
			paired-end data fastqc > fastp > HISAT2 > stringtie novel transcript identification

## small RNA-seq
There is currently 1 small RNA-seq analysis pipeline available. This is designed to work with QIAseq miRNA Library Kit from Qiagen.

	1.	sRNAseq_SE
			single-end data fastqc > umi-tools > STAR > featurecounts

## DNA Binding/enrichment 
There are currently 3 analysis pipelines available

	1.	ChIPseq_PE
			paired-end data fastqc > fastp > bowtie2 > macs2
	2.	CUT-RUN_PE
			paired-end data fastqc > fastp > bowtie2 > seacr & macs2
	3.	ATACseq_PE
			paired-end data fastqc > fastp > bowtie2 > macs2


# Description of files:

## Snakefiles
workflow/Snakefile launches the individual pipelines in workflow/rules

## config/samples_info.tab

	This file contains a tab deliminated table with:
		1. The names of R1 and R2 of each fastq file as received from the sequencing center. If sample was split over multiple lanes, remove the lane number ('L00X') from the fastq file name. cat_rename.py removes this when it concatenates .fastq files split over multiple lanes.
		2. Simple sample names
		3. Condition (e.g. diabetic vs non_diabetic)
		4. Replicate #
		5. If using ChIPseq or CUT-RUN a column titles 'antibody' is required. antibody specifies if the sample is the ChIP antibody or a control (input or IgG etc...)
		6. Sample name is the concatenated final sample_id. This is a concatenation of the sample name, condition, replicate, and antibody (if present) columns  
		7. Additional metadata can be added to this table for downstream analysis
		8. For ChIPseq and CUT-RUN, sample name, condition, and replicate should be identical for each pair of antibody and control fastq files. The antibody column specifies which of the pair is antibody and which is control.

## config/config.yaml
This file contains required general and workflow specific configuaration info.
	
	Generic requirements
		sample_file: Where to locate the samples_info.tab file (default config/samples_info.tab)
		workflow: name of workflow being used
		genome: location of indexed genome. 
			1. For RNAseq_PE, RNAseq_SE, or sRNAseq_SE - star 2.7.7a index
			2. For HISAT2 workflows - HISAT2 index
			3. For ChIPseq/CUT-RUN/ATACseq - bowtie2 index
		GTF: location of .gtf file
	CUT-RUN_PE
		spike_genome: Location of spike-in genome index. This is only implemented in CUT-RUN. bowtie2 index
		chromosome_lengths: location of chromosome lengths file. required for spike-in normalization in CUT-RUN
		effective_genome_size: Effective genome size for MACS2
	ChIPseq_PE
		effective_genome_size: Effective genome size for MACS2
	ATACseq_PE
		effective_genome_size: Effective genome size for MACS2
	RNAseq_HISAT2_stringtie or RNAseq_HISAT2_stringtie_nvltrx
		prepDE_length: Average fragment length for stringtie prepDE script

## config/profile/config.yaml
This file contains the default slurm resources for each rule

## workflow/scripts/cat_rename.py

	This script:
		1. Concatenates fastq files for samples that were split over multiple sequencing lanes
		2. Renames the fastq files from the generally verbose ids given by the sequencing center to those supplied in Samples_info.tab.
		3. The sample name, condition, and replicate columns are concatenated and form the new sample_id_Rx.fastq.gz files
		4. This script is executed via snakemake_init.sh prior to launching the appropriate snakemake pipeline
	Skip this script with the -c option when launching pipeline with snakemake_init.sh

## workflow/scripts/snakemake_init.sh
	
	This bash script:
		1. Executes cat_rename.py
		2. executes conda_load script
		3. Executes snakemake
		4. Runs multiqc

## workflow/scripts/launch_sbatch.sh
This script will launch the pipeline from a compute node vs a login node. We should always do this. Edit the snakemake_init.sh command in the script with desired parameters and launch via sbatch. 

## workflow/scripts/condaload_CVRCseq.sh
This script sets some environment variable and loads the conda environment.

## workflow/scripts/FRP.py
This file computes the fraction of reads in peaks (FRP) and outputs a table with FRP, total fragments, and fragments within peaks.

## workflow/envs/CVRCseq.yml
This file contains the info for the conda environment used by this pipeline.
 
## Usage
	
	When starting a new project:
		1. Clone the git repo using 'git clone https://github.com/mgildea87/CVRCseq.git'
		2. Update the samples_info.tab file with fastq.gz file names and desired sample, condition, replicate names, and Antibody/IgG control status (if using)
		3. Update config.yaml
		4. Modify parameters in the appropriate worklow/rules .smk file if desired. e.g. alignment parameters. 
		5. Run 'bash workflow/scripts/snakemake_init.sh'
			Description of parameters
				-h	help"
				-d	.fastq directory"
				-s	parameters to pass to snakemake (e.g. --unlock)
				-w	workflow name (e.g. 'RNAseq_PE')
				-c	Skip cat_rename.py. Use to skip copying, concatenating, and renaming of .fastq files to local directory

# To-do

	* Add testing data and tests
	* Add better error reporting. It's difficult sometimes to ID why the pipeline has failed.
	* Enrichment pipelines
		* Add parameter for specification of MACS2 narrow vs broad. Add to config file. Would need to edit snakefiles and FRP script
		* Add irreproducible discovery rate (IDR) for identifying robust peak sets between replicates. See ENCODE pipeline
		* Add deduplication as option for ChIPseq analysis.
		* enable more efficient handling of experimental designs where the same input is used for multiple pull-down/antibody samples. e.g. ChIRPseq.
		* Incorporate spike-in normalization to CUT&RUN macs2 peak calling
	* Simplify cat_rename.py to take sample prefixes (text upstream of the lane number '_L00X') supplied via samples_info.tab. see: https://github.com/macs3-project/MACS/issues/356
	* Add parameter to specify output directory name. Right now its given the pipeline name.
	* Add salmon pipeline for RNAseq
	* Add rules to snakefiles containing R scripts for some downstream QC and plotting. e.g. for RNAseq: PCA, replicate scatter plots, count statistics or for ATACseq: fragment length distributions, FRP plots, replicate comparisons. 

