import pandas as pd
import os


for directory in ['RNAseq_SE/results', 'RNAseq_SE/results/fastqc','RNAseq_SE/results/fastqc_post_trim', 'RNAseq_SE/results/trim', 'RNAseq_SE/results/logs', 'RNAseq_SE/results/logs/trim_reports', 'RNAseq_SE/results/alignment', 'RNAseq_SE/results/logs/alignment_reports', 'RNAseq_SE/results/feature_counts']:
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
	threads: 1
	params:
		'RNAseq_SE/results/fastqc/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule trim:
	input:
		R1='inputs/fastq/{sample}_R1.fastq.gz',
	output:
		R1='results/trim/{sample}_trimmed_R1.fastq.gz',
		html='results/logs/trim_reports/{sample}.html',
		json='results/logs/trim_reports/{sample}.json'
	threads: 16
	resources: 
		time_min=240, mem_mb=20000
	log:
		'results/logs/trim_reports/{sample}.log'
	shell:
		'fastp -w {threads} {params} -i {input.R1} -o {output.R1} --html {output.html} --json {output.json} 2> {log}'

rule fastqc_post_trim:
	input: 
		fastq = "results/trim/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc_post_trim/{sample}{read}_fastqc.html"
	threads: 1
	params:
		'RNAseq_SE/results/fastqc_post_trim/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule align:
	input:
		R1='results/trim/{sample}_trimmed_R1.fastq.gz'
	output:
		bam = 'results/alignment/{sample}.bam'
	threads: 16
	resources: 
		time_min=240, mem_mb=60000
	params:
		'--readFilesCommand zcat --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate --alignMatesGapMax 1000000 --outFilterMismatchNmax 999 --alignIntronMax 1000000 ' 
		'--alignSplicedMateMapLmin 3 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04 --outSAMattributes All --alignIntronMin 20 '
		'--outFilterIntronMotifs RemoveNoncanonicalUnannotated --scoreGapNoncan -14 --outSJfilterReads Unique --outFilterMultimapNmax 10'
	shell:
		'STAR {params} --genomeDir %s --runThreadN {threads} --readFilesIn {input.R1} --outFileNamePrefix RNAseq_SE/results/alignment/{wildcards.sample}_ | samtools view -bh > RNAseq_SE/results/alignment/{wildcards.sample}.bam' % (genome)
rule count:
	input:
		bam = expand('results/alignment/{sample}.bam', sample = sample_ids)
	output:
		counts = 'results/feature_counts/count_table.txt'
	threads: 16
	resources: 
		time_min=480, mem_mb=30000
	params:
		'-g gene_id -s 2 -Q 5 --extraAttributes gene_type,gene_name'
	shell:
		'featureCounts {params} -T {threads} -a %s -o {output.counts} {input.bam}' % (GTF)