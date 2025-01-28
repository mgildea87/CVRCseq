import pandas as pd
import os


for directory in ['sRNAseq_SE/results', 'sRNAseq_SE/results/fastqc','sRNAseq_SE/results/fastqc_post_trim', 'sRNAseq_SE/results/umi_tools_trim', 'sRNAseq_SE/results/logs', 'sRNAseq_SE/results/logs/umi_tools_trim_reports', 'sRNAseq_SE/results/alignment', 'sRNAseq_SE/results/logs/alignment_reports', 'sRNAseq_SE/results/feature_counts']:
	if not os.path.isdir(directory):
		os.mkdir(directory)

sample_file = config["sample_file"]
GTF = config["GTF"]
table = pd.read_table(sample_file)
sample = table['Sample']
replicate = table['Replicate']
condition = table['Condition']
File_R1 = table['File_Name_R1']
genome = config["genome"]

sample_ids = []
for i in range(len(sample)):
	sample_ids.append('%s_%s_%s' % (sample[i], condition[i], replicate[i]))

read = ['_R1']


rule all:
	input:
		'results/feature_counts/count_table.txt',
		expand('results/fastqc/{sample}{read}_fastqc.html', sample = sample_ids, read = read),
		expand('results/fastqc_post_trim/{sample_file}_trimmed{read}_fastqc.html', sample_file = sample_ids, read = read)

rule fastqc:
	input: 
		fastq = "inputs/fastq/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc/{sample}{read}_fastqc.html",
		"results/fastqc/{sample}{read}_fastqc.zip"
	params:
		'sRNAseq_SE/results/fastqc/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule umi_tools_trim:
	input:
		R1='inputs/fastq/{sample}_R1.fastq.gz'
	output:
		R1='results/umi_tools_trim/{sample}_trimmed_R1.fastq.gz'
	resources: 
		time_min=120, mem_mb=10000, cpus=1
	log:
		'results/logs/umi_tools_trim_reports/{sample}.log'
	shell:
		'umi_tools extract --extract-method=regex --bc-pattern=".+(?P<discard_1>AACTGTAGGCACCATCAAT){{s<=2}}(?P<umi_1>.{{12}})(?P<discard_2>.+)" -I {input.R1} -S {output.R1} -L {log}'

rule fastqc_post_trim:
	input: 
		fastq = "results/umi_tools_trim/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc_post_trim/{sample}{read}_fastqc.html",
	params:
		'sRNAseq_SE/results/fastqc_post_trim/'
	shell: 
		'fastqc {input.fastq} -o {params}' 

rule align:
	input:
		R1='results/umi_tools_trim/{sample}_trimmed_R1.fastq.gz'
	output:
		bam = 'results/alignment/{sample}.bam'
	threads: 10
	resources: 
		time_min=240, mem_mb=60000, cpus=20
	params:
		'--readFilesCommand zcat --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate --alignEndsType EndToEnd --outFilterMismatchNmax 1'
		'--outFilterMultimapScoreRange 0 --outFilterMultimapNmax 10 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0'
		'--outFilterMatchNmin 16 --alignSJDBoverhangMin 1000 --alignIntronMax 1'
	shell:
		'STAR {params} --genomeDir %s --runThreadN {threads} --readFilesIn {input.R1} --outFileNamePrefix sRNAseq_SE/results/alignment/{wildcards.sample}_ | samtools view -bh > sRNAseq_SE/results/alignment/{wildcards.sample}.bam' % (genome)
rule count:
	input:
		bam = expand('results/alignment/{sample}.bam', sample = sample_ids)
	output:
		counts = 'results/feature_counts/count_table.txt'
	threads: 20
	resources: 
		time_min=480, mem_mb=30000, cpus=20
	params:
		'-g gene_id -s 1 -Q 5 -F GTF --extraAttributes gene_type,gene_name'
	shell:
		'featureCounts {params} -T {threads} -a %s -o {output.counts} {input.bam}' % (GTF)
