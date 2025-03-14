import pandas as pd
import os

for directory in ['RNAseq_PE_HISAT2_stringtie_nvltrx/results', 'RNAseq_PE_HISAT2_stringtie_nvltrx/results/fastqc', 'RNAseq_PE_HISAT2_stringtie_nvltrx/results/trim', 'RNAseq_PE_HISAT2_stringtie_nvltrx/results/logs', 'RNAseq_PE_HISAT2_stringtie_nvltrx/results/logs/trim_reports', 'RNAseq_PE_HISAT2_stringtie_nvltrx/results/alignment', 'RNAseq_PE_HISAT2_stringtie_nvltrx/results/logs/alignment_reports', 'RNAseq_PE_HISAT2_stringtie_nvltrx/results/stringtie', 'RNAseq_PE_HISAT2_stringtie_nvltrx/results/stringtie/merged']:
	if not os.path.isdir(directory):
		os.mkdir(directory)

configfile: "config.yaml"
sample_file = config["sample_file"]
GTF = config["GTF"]
prepDE_length = config["prepDE_length"]

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
		'results/stringtie/merged/gene_count_matrix.csv',
		'results/stringtie/merged/transcript_count_matrix.csv',
		expand('results/fastqc/{sample}{read}_fastqc.html', sample = sample_ids, read = read)

rule fastqc:
	input: 
		fastq = "inputs/fastq/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc/{sample}{read}_fastqc.html",
		"results/fastqc/{sample}{read}_fastqc.zip"
	params:
		'RNAseq_PE_HISAT2_stringtie_nvltrx/results/fastqc/'
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
	shell:
		'fastp -w {threads} {params} -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} 2> {log}'

rule align:
	input:
		R1='results/trim/{sample}_trimmed_R1.fastq.gz',
		R2='results/trim/{sample}_trimmed_R2.fastq.gz'
	output:
		bam = 'results/alignment/{sample}.bam'
	threads: 16
	resources: 
		time_min=240, mem_mb=60000, cpus=16
	log:
		'results/logs/alignment_reports/{sample}.log'
	params:
		'--phred33 --rna-strandness RF --dta'
	shell:
		'hisat2 {params} -p {threads} -x %s -1 {input.R1} -2 {input.R2} 2> {log} | samtools sort - -o RNAseq_PE_HISAT2_stringtie_nvltrx/results/alignment/{wildcards.sample}.bam -@ {threads}' % (genome)

rule count:
	input:
		bam = 'results/alignment/{sample}.bam'
	output:
		trans_counts = 'results/stringtie/{sample}.gtf',
		gene_counts = 'results/stringtie/{sample}.tab'
	threads: 16
	resources: 
		time_min=240, mem_mb=40000, cpus=16
	params:
		'--rf --conservative'
	shell:
		'stringtie -p {threads} {params} -G %s -o {output.trans_counts} -l {wildcards.sample} -A {output.gene_counts} {input.bam}' % (GTF)

rule combine_gtf:
	input:
		trans_counts = expand('results/stringtie/{sample}.gtf' , sample = sample_ids)
	output:
		gtf_list = 'results/stringtie/merged.txt'
	shell:
		'ls RNAseq_PE_HISAT2_stringtie_nvltrx/results/stringtie/*.gtf > {output.gtf_list}'

rule merge:
	input:
		gtf_list = 'results/stringtie/merged.txt'
	output:
		merged_gtf = 'results/stringtie/merged.gtf'
	threads: 16
	resources: 
		time_min=240, mem_mb=40000, cpus=16
	shell:
		'stringtie --merge -p {threads} -G %s -o {output.merged_gtf} {input.gtf_list}' % (GTF)

rule count_2:
	input:
		bam = 'results/alignment/{sample}.bam',
		merged_gtf = 'results/stringtie/merged.gtf'
	output:
		trans_counts = 'results/stringtie/merged/{sample}/{sample}_merged.gtf'
	threads: 16
	resources: 
		time_min=240, mem_mb=40000, cpus=16
	params:
		'--rf -e -B'
	shell:
		'stringtie -p {threads} {params} -G  {input.merged_gtf} -o {output.trans_counts} {input.bam}'

rule deseq_prep:
	input:
		expand('results/stringtie/merged/{sample}/{sample}_merged.gtf', sample = sample_ids),
		str_dir = 'results/stringtie/merged/' 
	output:
		gene_counts = 'results/stringtie/merged/gene_count_matrix.csv',
		transcript_counts = 'results/stringtie/merged/transcript_count_matrix.csv'
	shell:
		'prepDE.py -l %s -i {input.str_dir} -g {output.gene_counts} -t {output.transcript_counts}' % (prepDE_length)