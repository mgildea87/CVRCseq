import pandas as pd
import os

for directory in ['ChIPseq_PE/results', 'ChIPseq_PE/results/fastqc', 'ChIPseq_PE/results/fastqc_post_trim', 'ChIPseq_PE/results/trim', 'ChIPseq_PE/results/logs', 'ChIPseq_PE/results/logs/trim_reports', 'ChIPseq_PE/results/alignment', 'ChIPseq_PE/results/alignment/frag_len', 'ChIPseq_PE/results/logs/alignment_reports', 'ChIPseq_PE/results/peaks', 'ChIPseq_PE/results/logs/MACS2']:
	if not os.path.isdir(directory):
		os.mkdir(directory)

sample_file = config["sample_file"]
genome = config["genome"]
effective_genome_size = config["effective_genome_size"]

table = pd.read_table(sample_file)
sample = table['Sample']
replicate = table['Replicate']
condition = table['Condition']
Antibody = table['Antibody']
File_R1 = table['File_Name_R1']
File_R2 = table['File_Name_R2']
File_names = File_R1.append(File_R2)

sample_ids = []
for i in range(len(sample)):
	sample_ids.append('%s_%s_%s' % (sample[i], condition[i], replicate[i]))
sample_ids = pd.unique(sample_ids).tolist()

sample_ids_file = []
for i in range(len(sample)):
	sample_ids_file.append('%s_%s_%s_%s' % (sample[i], condition[i], replicate[i], Antibody[i]))

read = ['_R1', '_R2']

rule all:
	input:
		expand('results/fastqc/{sample_file}{read}_fastqc.html', sample_file = sample_ids_file, read = read),
		expand('results/fastqc_post_trim/{sample_file}_trimmed{read}_fastqc.html', sample_file = sample_ids_file, read = read),
		expand('results/peaks/{sample}_peaks.narrowPeak', sample = sample_ids),
		'results/FRP.txt',
		expand('results/alignment/{sample}_sorted.bam', sample = sample_ids_file),
		expand('results/alignment/{sample}_sorted.bam.bai', sample = sample_ids_file),
		expand('results/alignment/frag_len/{sample}.txt', sample = sample_ids_file)

rule fastqc:
	input: 
		fastq = "inputs/fastq/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc/{sample}{read}_fastqc.html"
	threads: 1
	params:
		'ChIPseq_PE/results/fastqc/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule fastqc_post_trim:
	input: 
		fastq = "results/trim/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc_post_trim/{sample}{read}_fastqc.html"
	threads: 1
	params:
		'ChIPseq_PE/results/fastqc_post_trim/'
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
		time_min=240, mem_mb=20000
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
		'results/alignment/{sample}.bam'
	threads: 16
	resources: 
		time_min=240, mem_mb=60000
	log:
		'results/logs/alignment_reports/{sample}.log'
	params:
		'--end-to-end --very-sensitive --no-mixed --no-unal --no-discordant --phred33'
	shell:
		'bowtie2 {params} -x %s --threads {threads} -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -bh -q 3 > ChIPseq_PE/results/alignment/{wildcards.sample}.bam' % (genome)

rule sort:
	input:
		'results/alignment/{sample}.bam'
	output:
		'results/alignment/{sample}_sorted.bam'	
	threads: 16
	resources: 
		time_min=240, mem_mb=20000
	shell:
		'samtools sort -@ {threads} {input} > {output}'

rule index:
	input:
		'results/alignment/{sample}_sorted.bam'
	output:
		'results/alignment/{sample}_sorted.bam.bai'	
	threads: 16
	resources: 
		time_min=240, mem_mb=20000
	shell:
		'samtools index -@ {threads} {input} > {output}'		

rule MACS2:
  	input:
  		exp='results/alignment/{sample}_Antibody.bam',
  		con='results/alignment/{sample}_Control.bam'
  	output:
  		'results/peaks/{sample}_peaks.narrowPeak'
  	threads: 1
  	log:
  		'results/logs/MACS2/{sample}.log'
  	params:
  		'-B --outdir ChIPseq_PE/results/peaks/ -g %s -q 0.05 --keep-dup auto -f BAMPE' % (effective_genome_size)
  	shell:
  		'macs2 callpeak -t {input.exp} -c {input.con} {params} -n {wildcards.sample} 2> {log}'

rule fragment_size:
 	input:
 		'results/alignment/{sample}.bam'
 	output:
 		'results/alignment/frag_len/{sample}.txt'
 	threads: 1
 	shell:
 		"""
 		samtools view {input} | awk -F'\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{print abs($9)}}' | sort | uniq -c | awk -v OFS="\t" '{{print $2, $1/2}}' > {output}
 		"""

rule FRP:
 	input:
 		expand('results/peaks/{sample}_peaks.narrowPeak', sample = sample_ids)
 	output:
 		'results/FRP.txt'
 	threads: 1
 	script:
 		'../scripts/FRP.py'
