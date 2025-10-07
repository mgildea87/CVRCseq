import pandas as pd
import os

for directory in ['RNAseq_PE/results', 'RNAseq_PE/results/fastqc','RNAseq_PE/results/fastqc_post_trim', 'RNAseq_PE/results/trim', 'RNAseq_PE/results/logs', 'RNAseq_PE/results/logs/trim_reports', 'RNAseq_PE/results/alignment', 'RNAseq_PE/results/logs/alignment_reports', 'RNAseq_PE/results/feature_counts']:
	if not os.path.isdir(directory):
		os.mkdir(directory)

sample_file = config["sample_file"]
GTF = config["GTF"]
table = pd.read_table(sample_file)
sample = table['Sample']
replicate = table['Replicate']
condition = table['Condition']
File_R1 = table['File_Name_R1']
File_R2 = table['File_Name_R2']
File_names = File_R1.append(File_R2)
genome = config["genome"]

sample_ids = []
for i in range(len(sample)):
	sample_ids.append('%s_%s_%s' % (sample[i], condition[i], replicate[i]))

read = ['_R1', '_R2']


rule all:
	input:
		'results/feature_counts/count_table.txt',
		expand('results/fastqc/{sample}{read}_fastqc.html', sample = sample_ids, read = read),
		expand('results/fastqc_post_trim/{sample}_trimmed{read}_fastqc.html', sample = sample_ids, read = read)

rule fastqc:
	input: 
		fastq = "inputs/fastq/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc/{sample}{read}_fastqc.html",
		"results/fastqc/{sample}{read}_fastqc.zip"
	params:
		'RNAseq_PE/results/fastqc/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule trim:
	input:
		R1='inputs/fastq/{sample}_R1.fastq.gz',
		R2='inputs/fastq/{sample}_R2.fastq.gz'
	output:
		R1='results/trim/{sample}_trimmed_R1.fastq.gz',
		R2='results/trim/{sample}_trimmed_R2.fastq.gz',
		html='results/logs/trim_reports/{sample}.html',
		json='results/logs/trim_reports/{sample}.json'
	threads: 16
	resources: 
		time_min=240, mem_mb=20000, cpus=16
	log:
		'results/logs/trim_reports/{sample}.log'
	params:
		'--detect_adapter_for_pe'
#		'--adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT'
	shell:
		'fastp -w {threads} {params} -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} 2> {log}'

rule fastqc_post_trim:
	input: 
		fastq = "results/trim/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc_post_trim/{sample}{read}_fastqc.html",
	params:
		'RNAseq_PE/results/fastqc_post_trim/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule align:
	input:
		R1='results/trim/{sample}_trimmed_R1.fastq.gz',
		R2='results/trim/{sample}_trimmed_R2.fastq.gz'
	output:
		bam = 'results/alignment/{sample}.bam'
	threads: 16
	resources: 
		time_min=240, mem_mb=60000, cpus=16
	params:
		'--readFilesCommand zcat --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate --alignMatesGapMax 1000000 --outFilterMismatchNmax 999 --alignIntronMax 1000000 ' 
		'--alignSplicedMateMapLmin 3 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04 --outSAMunmapped Within KeepPairs --outSAMattributes All --alignIntronMin 20 '
		'--outFilterIntronMotifs RemoveNoncanonicalUnannotated --scoreGapNoncan -14 --outSJfilterReads Unique --outFilterMultimapNmax 10'
	shell:
		'STAR {params} --genomeDir %s --runThreadN {threads} --readFilesIn {input.R1} {input.R2} --outFileNamePrefix RNAseq_PE/results/alignment/{wildcards.sample}_ | samtools view -bh > RNAseq_PE/results/alignment/{wildcards.sample}.bam' % (genome)

rule count:
	input:
		bam = expand('results/alignment/{sample}.bam', sample = sample_ids)
	output:
		counts = 'results/feature_counts/count_table.txt'
	threads: 16
	resources: 
		time_min=480, mem_mb=30000, cpus=16
	params:
		'-p --countReadPairs -g gene_id -s 2 -Q 5 --extraAttributes gene_type,gene_name'
	shell:
		'featureCounts {params} -T {threads} -a %s -o {output.counts} {input.bam}' % (GTF)

