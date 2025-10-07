import pandas as pd
import os

for directory in ['CUT-RUN_PE/results', 'CUT-RUN_PE/results/fastqc', 'CUT-RUN_PE/results/fastqc_post_trim', 'CUT-RUN_PE/results/trim', 'CUT-RUN_PE/results/logs', 'CUT-RUN_PE/results/logs/MACS2', 'CUT-RUN_PE/results/logs/trim_reports', 'CUT-RUN_PE/results/alignment','CUT-RUN_PE/results/alignment/bed', 'CUT-RUN_PE/results/alignment/frag_len', 'CUT-RUN_PE/results/logs/alignment_reports', 'CUT-RUN_PE/results/peaks', 'CUT-RUN_PE/results/peaks/seacr', 'CUT-RUN_PE/results/peaks/MACS2']:
	if not os.path.isdir(directory):
		os.mkdir(directory)

sample_file = config["sample_file"]
genome = config["genome"]
spike_genome = config["spike_genome"]
chr_lens = config["chromosome_lengths"]
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
		expand('results/peaks/seacr/{sample}.stringent.bed', sample = sample_ids),
		'results/FRP_seacr.txt',
		expand('results/peaks/MACS2/{sample}_peaks.broadPeak', sample = sample_ids),
		'results/FRP_MACS2.txt',
		expand('results/alignment/frag_len/{sample}.txt', sample = sample_ids_file),
		expand('results/alignment/{sample}_sorted.bam.bai', sample = sample_ids_file),
		expand('results/alignment/{sample}_sorted.bam', sample = sample_ids_file)

rule fastqc:
	input: 
		fastq = "inputs/fastq/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc/{sample}{read}_fastqc.html",
	params:
		'CUT-RUN_PE/results/fastqc/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule fastqc_post_trim:
	input: 
		fastq = "results/trim/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc_post_trim/{sample}{read}_fastqc.html",
	params:
		'CUT-RUN_PE/results/fastqc_post_trim/'
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
	threads: 20
	resources: 
		time_min=240, mem_mb=10000, cpus=20
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
	threads: 30
	resources: 
		time_min=240, mem_mb=60000, cpus=30
	log:
		'results/logs/alignment_reports/{sample}.log'
	params:
		'--end-to-end --very-sensitive --no-mixed --no-unal --no-discordant --phred33 -I 10 -X 700'
	shell:
		'bowtie2 {params} -x %s --threads {threads} -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -bh -q 3 > CUT-RUN_PE/results/alignment/{wildcards.sample}.bam' % (genome)

rule align_spike:
	input:
		R1='results/trim/{sample}_trimmed_R1.fastq.gz',
		R2='results/trim/{sample}_trimmed_R2.fastq.gz'
	output:
		'results/alignment/{sample}_ecoli.bam'
	threads: 30
	resources: 
		time_min=240, mem_mb=60000, cpus=30
	log:
		'results/logs/alignment_reports/{sample}_ecoli.log'
	params:
		'--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-unal --no-discordant --phred33 -I 10 -X 700'
	shell:
		'bowtie2 {params} -x %s --threads {threads} -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -bh -q 3 > CUT-RUN_PE/results/alignment/{wildcards.sample}_ecoli.bam' % (spike_genome)

rule sort:
	input:
		'results/alignment/{sample}.bam'
	output:
		'results/alignment/{sample}_sorted.bam'	
	threads: 20
	resources: 
		time_min=240, mem_mb=20000, cpus=20
	shell:
		'samtools sort -@ {threads} {input} > {output}'

rule index:
	input:
		'results/alignment/{sample}_sorted.bam'
	output:
		'results/alignment/{sample}_sorted.bam.bai'	
	threads: 20
	resources: 
		time_min=240, mem_mb=20000, cpus=20
	shell:
		'samtools index -@ {threads} {input} > {output}'

rule spike_in_norm:
	input:
		sample_bam='results/alignment/{sample}.bam',
		spike_bam='results/alignment/{sample}_ecoli.bam'
	output:
		'results/alignment/bed/{sample}.bed',
		'results/alignment/bed/{sample}.bedgraph'
	shell:
		"""
		depth=`samtools view CUT-RUN_PE/results/alignment/{wildcards.sample}_ecoli.bam | wc -l`
		depth=$((depth/2))
		echo $depth
		scale_fac=`echo "10000 / $depth" | bc -l`
		echo $scale_fac
		bedtools bamtobed -bedpe -i CUT-RUN_PE/results/alignment/{wildcards.sample}.bam | cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n > CUT-RUN_PE/results/alignment/bed/{wildcards.sample}.bed
		"""
		'bedtools genomecov -bg -i CUT-RUN_PE/results/alignment/bed/{wildcards.sample}.bed -scale $scale_fac -g %s > CUT-RUN_PE/results/alignment/bed/{wildcards.sample}.bedgraph' % (chr_lens)

rule SEACR:
	input:
		exp='results/alignment/bed/{sample}_Antibody.bedgraph',
		con='results/alignment/bed/{sample}_Control.bedgraph'
	output:
		'results/peaks/seacr/{sample}.stringent.bed'
	resources: 
		time_min=120, mem_mb=40000
	params:
		'non stringent'
	shell:
		'bash SEACR_1.3.sh {input.exp} {input.con} {params} CUT-RUN_PE/results/peaks/seacr/{wildcards.sample}'

rule MACS2:
	input:
		exp='results/alignment/{sample}_Antibody.bam',
		con='results/alignment/{sample}_Control.bam'
	output:
		'results/peaks/MACS2/{sample}_peaks.broadPeak'
	resources: 
		time_min=120, mem_mb=40000
	log:
		'results/logs/MACS2/{sample}.log'
	params:
		'-B --outdir CUT-RUN_PE/results/peaks/MACS2/ -g %s -q 0.05 -f BAMPE --broad' % (effective_genome_size)
	shell:
		'macs2 callpeak -t {input.exp} -c {input.con} {params} -n {wildcards.sample} 2> {log}'

rule fragment_size:
	input:
		'results/alignment/{sample}.bam'
	output:
		'results/alignment/frag_len/{sample}.txt'
	shell:
		"""
		samtools view {input} | awk -F'\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{print abs($9)}}' | sort | uniq -c | awk -v OFS="\t" '{{print $2, $1/2}}' > {output}
		"""

rule FRP_seacr:
	input:
		expand('results/peaks/seacr/{sample}.stringent.bed', sample = sample_ids)
	output:
		'results/FRP_seacr.txt'
	resources: 
		time_min=300, mem_mb=40000
	params:
		peak_caller='seacr'
	script:
		'../scripts/FRP.py'

rule FRP_MACS2:
	input:
		expand('results/peaks/MACS2/{sample}_peaks.broadPeak', sample = sample_ids)
	output:
		'results/FRP_MACS2.txt'
	resources: 
		time_min=300, mem_mb=40000
	params:
		peak_caller='MACS2'
	script:
		'../scripts/FRP.py'

