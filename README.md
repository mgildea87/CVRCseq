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

## DNA Binding/enrichment 
There are currently 2 analysis pipelines available

	1.	ChIPseq
			paired-end data fastqc > fastp > bowtie2 > macs2
	2.	CUT-RUN
			paired-end data fastqc > fastp > bowtie2 > seacr


# Description of files:

## Snakefiles
workflow/Snakefile contains the master workflow which launches the individual pipelines in workflow/rules

## config/samples_info.tab
This file contains a tab deliminated table with:

		1. The names of R1 and R2 of each fastq file as received from the sequencing center. If sample was not split over multiple lanes, remove the lane number (L001) from the fastq file name. cat_rename.py removes this when it concatenates .fastq files split over multiple lanes.
		2. Simple sample names
		3. Condition (e.g. diabetic vs non_diabetic)
		4. Replicate #
		5. If using ChIPseq or CUT-RUN a column titles 'antibody' is required. antibody specifies if the sample is the ChIP antibody or a control (input or IgG etc...)
		6. Sample name is the concatenated final sample_id. This is a concatenation of the sample name, condition, replicate, and antibody (if present) columns  
		7. Additional metadata can be added to this table for downstream analysis
		8. For ChIPseq and CUT-RUN, sample name, condition, and replicate should be identical for each pair of antibody and control fastq files. The antibody column specifies which of the pair is antibody and which is control.

## config/config.yaml
This file contains required general and workflow specific configuaration info.
	
	Generic
		sample_file: Where to locate the samples_info.tab file (should be config/samples_info.tab)
		genome: location of indexed genome. 
			1. For RNAseq_PE or RNAseq_SE - star 2.7.7a index
			2. For HISAT2 workflows - HISAT2 index
			3. For ChIPseq or CUT-RUN - bowtie2 index
		GTF: location of .gtf file
	CUT-RUN
		spike_genome: Location of spike-in genome index. This is only implemented in CUT-RUN. bowtie2
		chromosome_lengths: location of chromosome lengths file. required for spike-in normalization in CUT-RUN
	ChIPseq
		effective_genome_size: Effective genome size for MACS2
	RNAseq_HISAT2_stringtie or RNAseq_HISAT2_stringtie_nvltrx
		prepDE_length: Average fragment length for stringtie prepDE script

## config/profile/config.yaml
This file contains the default slurm parameters for each rule

## workflow/scripts/cat_rename.py
This script:

		1. Concatenates fastq files for samples that were split over multiple sequencing lanes
		2. Renames the fastq files from the generally verbose ids given by the sequencing center to those supplied in the Samples_info.tab file.
		3. The sample name, condition, and replicate columns are concatenated and form the new sample_id_Rx.fastq.gz files
		4. This script is executed snakemake_init.sh prior to snakemake execution
## workflow/scripts/snakemake_init.sh
This bash script:

		1. loads the miniconda3/cpu/4.9.2 module
		3. Executes snakemake

## workflow/scripts/FRP.py
This file computes the fraction of reads in peaks (FRP) and outputs a table with FRP, total fragments, and fragments within peaks.

## workflow/envs/CVRCseq.yml
This file contains the conda environment info used by this pipeline.
 
## Usage
When starting a new project:

		1. Clone the git repo using 'git clone https://github.com/mgildea87/CVRCseq.git'
		2. Update the samples_info.tab file with fastq.gz file names and desired sample, condition, replicate names, and Antibody/IgG control status (if using)
		3. Update config.yaml
		4. Run 'bash test/snakemake_init.sh'
			Description of parameters
				-h	help"
				-d	.fastq directory"
				-c	Absolute path to conda environment"
				-w	workflow
