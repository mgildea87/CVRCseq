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
			1. For RNAseq_PE or RNAseq_SE - star 2.7.7a index
			2. For HISAT2 workflows - HISAT2 index
			3. For ChIPseq or CUT-RUN - bowtie2 index
		GTF: location of .gtf file
	CUT-RUN
		spike_genome: Location of spike-in genome index. This is only implemented in CUT-RUN. bowtie2 index
		chromosome_lengths: location of chromosome lengths file. required for spike-in normalization in CUT-RUN
	ChIPseq
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
## workflow/scripts/snakemake_init.sh
	
	This bash script:
		1. Executes cat_rename.py
		2. loads the miniconda3/cpu/4.9.2 module on Ultraviolet
		3. Executes snakemake
		4. Runs multiqc


## workflow/scripts/FRP.py
This file computes the fraction of reads in peaks (FRP) and outputs a table with FRP, total fragments, and fragments within peaks.

## workflow/envs/CVRCseq.yml
This file contains the conda environment info used by this pipeline.
 
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
				-c	Absolute path to conda environment"
				-w	workflow

# To-do

	* Add testing data and tests
	* Enrichment pipelines
		* Add ATACseq pipeline
		* Add irreproducible discovery rate (IDR) for identifying robust peak sets between replicates. See ENCODE pipeline
		* Add deduplication by default. Likely via Picard prior to peak calling. 
		* enable more efficient handling of experimental designs where the same input is used for multiple pull-down/antibody samples. e.g. ChIRPseq.
	* Simplify cat_rename.py to take sample prefixes (text upstream of the lane number '_L00X') supplied via samples_info.tab.
	* Add parameter to specify output directory name. Right now its given the pipeline name.
	* Add rules to snakefiles containing R scripts for some downstream QC and plotting. e.g. for RNAseq: PCA, replicate scatter plots, count statistics or for ATACseq: fragment length distributions, FRP plots, replicate comparisons. 

