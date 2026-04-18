import pandas as pd
import os

for directory in ['ATACseq_PE/results', 'ATACseq_PE/results/fastqc', 'ATACseq_PE/results/fastqc_post_trim', 'ATACseq_PE/results/trim', 'ATACseq_PE/results/logs', 'ATACseq_PE/results/logs/trim_reports', 'ATACseq_PE/results/alignment', 'ATACseq_PE/results/alignment/frag_len', 'ATACseq_PE/results/alignment/idxstat', 'ATACseq_PE/results/logs/alignment_reports', 'ATACseq_PE/results/peaks/MACS2','ATACseq_PE/results/peaks/MACS2/qc', 'ATACseq_PE/results/logs/MACS2']:
	if not os.path.isdir(directory):
		os.mkdir(directory)

sample_file = config["sample_file"]
genome = config["genome"]
effective_genome_size = config["effective_genome_size"]

table = pd.read_table(sample_file)
sample = table['Sample']
replicate = table['Replicate']
condition = table['Condition']
File_R1 = table['File_Name_R1']
File_R2 = table['File_Name_R2']
File_names = File_R1.append(File_R2)

sample_ids = []
for i in range(len(sample)):
	sample_ids.append('%s_%s_%s' % (sample[i], condition[i], replicate[i]))
sample_ids = pd.unique(sample_ids).tolist()

read = ['_R1', '_R2']

rule all:
	input:
		expand('results/fastqc/{sample_file}{read}_fastqc.html', sample_file = sample_ids, read = read),
		expand('results/fastqc_post_trim/{sample_file}_trimmed{read}_fastqc.html', sample_file = sample_ids, read = read),
		expand('results/peaks/MACS2/{sample}_peaks.narrowPeak', sample = sample_ids),
		expand('results/peaks/MACS2/qc/{sample}_frip_stats.txt', sample= sample_ids),
		expand('results/alignment/{sample}.bam', sample = sample_ids),
		expand('results/alignment/{sample}.bam.bai', sample = sample_ids),
		expand('results/alignment/{sample}_dedup.bam', sample = sample_ids),
		expand('results/alignment/{sample}_dedup.bam.bai', sample = sample_ids),
		expand('results/alignment/{sample}_filtered_sorted.bam', sample = sample_ids),
		expand('results/alignment/{sample}_filtered_sorted.bam.bai', sample = sample_ids),
		expand('results/alignment/{sample}_dedup_filtered_sorted.bam', sample = sample_ids),
		expand('results/alignment/{sample}_dedup_filtered_sorted.bam.bai', sample = sample_ids),
		expand('results/alignment/{sample}_dedup_filtered_sorted.bw', sample = sample_ids),
		expand('results/alignment/frag_len/{sample}.txt', sample = sample_ids),
		expand('results/alignment/idxstat/{sample}_idxstat.tab', sample = sample_ids),
		expand('results/alignment/idxstat/{sample}_dedup_idxstat.tab', sample = sample_ids),
		"results/peaks/MACS2/qc/frip_summary_detailed.tsv"

rule fastqc:
	input: 
		fastq = "inputs/fastq/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc/{sample}{read}_fastqc.html"
	threads: 1
	params:
		'ATACseq_PE/results/fastqc/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule fastqc_post_trim:
	input: 
		fastq = "results/trim/{sample}{read}.fastq.gz"
	output:  
		"results/fastqc_post_trim/{sample}{read}_fastqc.html"
	threads: 1
	params:
		'ATACseq_PE/results/fastqc_post_trim/'
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
	threads: 24
	resources: 
		time_min=719, mem_mb=60000
	log:
		'results/logs/alignment_reports/{sample}.log'
	params:
		'--end-to-end --very-sensitive --no-mixed --no-unal --no-discordant --phred33'
	shell:
		'bowtie2 {params} -x %s --threads {threads} -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -bh -q 3 -f 3 | samtools sort -@ {threads} > ATACseq_PE/results/alignment/{wildcards.sample}.bam' % (genome)

rule idxstat:
	input:
		bam='results/alignment/{sample}.bam'
	output:
		idx='results/alignment/idxstat/{sample}_idxstat.tab'
	threads: 1
	resources: 
		time_min=10, mem_mb=5000		
	shell:
		'samtools idxstat {input.bam} > {output.idx}'

rule dedup_bam:
	input:
		bam='results/alignment/{sample}.bam'
	output:
		'results/alignment/{sample}_dedup.bam'
	threads: 16
	resources: 
		time_min=120, mem_mb=30000
	shell:
		'samtools collate {input} -O -@ {threads} | samtools fixmate -m -@ {threads} - - | samtools sort -@ {threads} | samtools markdup - ATACseq_PE/results/alignment/{wildcards.sample}_dedup.bam -@ {threads} -rsS'


rule filter_bam:
	input:
		'results/alignment/{sample}.bam'
	output:
		'results/alignment/{sample}_filtered_sorted.bam'
	threads: 16
	resources: 
		time_min=240, mem_mb=30000
	shell:
		'samtools view -h {input} | grep -v chrM | samtools view -bh > {output}'


rule index:
	input:
		'results/alignment/{sample}.bam'
	output:
		'results/alignment/{sample}.bam.bai'	
	threads: 16
	resources: 
		time_min=240, mem_mb=30000
	shell:
		'samtools index -@ {threads} {input} > {output}'


rule bam2bed:
	input:
		'results/alignment/{sample}_dedup_filtered_sorted.bam'
	output:
		'results/alignment/{sample}_dedup_filtered_sorted.bed'
	threads: 1
	shell:
		"""
		bedtools bamtobed -i {input} | awk -F$'\t' 'BEGIN {{OFS = FS}}{{ if ($6 == "+") {{$2 = $2 + 4}} else if ($6 == "-") {{$3 = $3 - 5}} print $0}}' > {output}
		"""

rule MACS2:
	input:
		exp='results/alignment/{sample}_dedup_filtered_sorted.bed'
	output:
		'results/peaks/MACS2/{sample}_peaks.narrowPeak'
	threads: 1
	resources:
		time_min=240, mem_mb=30000
	log:
		'results/logs/MACS2/{sample}.log'
	params:
		'-B --outdir ATACseq_PE/results/peaks/MACS2/ -g %s -p 0.01 --keep-dup all -f BED --nomodel --shift -75 --extsize 150 --call-summits' % (effective_genome_size)
	shell:
		'macs2 callpeak -t {input.exp} {params} -n {wildcards.sample} 2> {log}'

rule bam2bw:
	input:
		BAM='results/alignment/{sample}_dedup_filtered_sorted.bam',
		BAI='results/alignment/{sample}_dedup_filtered_sorted.bam.bai'
	output:
		'results/alignment/{sample}_dedup_filtered_sorted.bw'
	threads: 16
	params:
	  '--normalizeUsing RPKM --outFileFormat bigwig --binSize 1'
	shell:
		"""
		bamCoverage -b {input.BAM} -o {output} {params} --numberOfProcessors {threads}
		"""

rule fragment_size:
	input:
		'results/alignment/{sample}_dedup_filtered_sorted.bam'
	output:
		'results/alignment/frag_len/{sample}.txt'
	threads: 1
	shell:
		"""
		samtools view {input} | awk -F'\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{print abs($9)}}' | sort | uniq -c | awk -v OFS="\t" '{{print $2, $1/2}}' > {output}
		"""

rule FRP:
	input:
		bam = "results/alignment/{sample}_dedup_filtered_sorted.bam",
		peaks = "results/peaks/MACS2/{sample}_peaks.narrowPeak"
	output:
		stats = "results/peaks/MACS2/qc/{sample}_frip_stats.txt"
	threads:8
	resources: 
		mem_mb=50000
	shell:
		"""
		# Count total mapped reads
		total_reads=$(samtools view -c -F 260 {input.bam})
    	total_fragments=$(( total_reads / 2 ))
		# Count reads overlapping peaks
		reads_in_peaks=$(samtools sort -n -@ {threads} -m 3G {input.bam} | bedtools bamtobed -bedpe -i stdin | bedtools intersect -a stdin -b {input.peaks} -u | wc -l)
        
		# Count total number of peaks called
		num_peaks=$(wc -l < {input.peaks})

		# Calculate FRiP
		frip=$(awk -v a="$reads_in_peaks" -v b="$total_fragments" \
		    'BEGIN {{ if (b>0) printf "%.4f", a/b; else print "0" }}')

		# Save all values to a single line
		echo -e "{wildcards.sample}\\t$total_fragments\\t$num_peaks\\t$reads_in_peaks\\t$frip" > {output.stats}
		"""

rule aggregate_qc_summary:
	input:
		# Collects all stats files from the previous step
		stats_files = expand("results/peaks/MACS2/qc/{sample}_frip_stats.txt", sample=sample_ids)
	output:
		summary = "results/peaks/MACS2/qc/frip_summary_detailed.tsv"
	script:
		"../scripts/aggregate_peak_qc.py"
