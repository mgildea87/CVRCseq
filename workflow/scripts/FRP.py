import subprocess, sys, pandas as pd, yaml
with open('config/config.yaml', 'r') as file:
	config = yaml.safe_load(file)
workflow = config['workflow']
sample_file = config['sample_file']

table = pd.read_table(sample_file)
sample = table['Sample']
replicate = table['Replicate']
condition = table['Condition']

#peak files get named different things with macs vs seacr will need to update if peak caller is changed or other workflows are added
if workflow == 'ChIPseq' or workflow == 'ATACseq_PE':
	peak_file_suffix = '_peaks.narrowPeak'
if workflow == 'CUT-RUN':
	peak_file_suffix = '.stringent.bed'

if workflow == 'ATACseq_PE':
	sample_ids = []
	for i in range(len(sample)):
		sample_ids.append('%s_%s_%s' % (sample[i], condition[i], replicate[i]))
	sample_ids = pd.unique(sample_ids).tolist()
else:
	Antibody = table['Antibody']
	sample_ids = []
	for i in range(len(sample)):
		sample_ids.append('%s_%s_%s_%s' % (sample[i], condition[i], replicate[i], Antibody[i]))

total_fragments = []
peak_fragments = []


def FRP(sample, workflow, peak_file_suffix):
	peak_id = sample.split('_')[0:4]
	peak_id = "_".join(peak_id)
	print(sample)
	cmd_1 = ['samtools', 'sort', '-n', '-@', '8', '-m', '3G', '%s/results/alignment/%s.bam' % (workflow, sample)] 
	cmd_2 = ['bedtools', 'bamtobed', '-bedpe', '-i', 'stdin'] 
	cmd_3 = ['bedtools', 'intersect', '-a', 'stdin', '-b', '%s/results/peaks/%s%s' % (workflow, peak_id, peak_file_suffix), '-u']
	cmd_4 = ['wc', '-l']
	step_1 = subprocess.Popen(cmd_1, stdout=subprocess.PIPE)
	step_2 = subprocess.Popen(cmd_2, stdin = step_1.stdout, stdout=subprocess.PIPE)
	step_3 = subprocess.Popen(cmd_3, stdin = step_2.stdout, stdout=subprocess.PIPE)
	step_4 = subprocess.Popen(cmd_4, stdin = step_3.stdout, stdout=subprocess.PIPE)
	frag_peaks = step_4.stdout.read()
	frag_peaks.strip()
	peak_fragments.append(float(frag_peaks))

	cmd_5 = ['samtools', 'view', '%s/results/alignment/%s.bam' % (workflow, sample)]
	step_1 = subprocess.Popen(cmd_5, stdout=subprocess.PIPE)
	step_2 = subprocess.Popen(cmd_4, stdin = step_1.stdout, stdout=subprocess.PIPE)

	frag = step_2.stdout.read()
	frag.strip()
	frag = float(frag)
	total_fragments.append(frag/2)

for samples in sample_ids:
	FRP(samples, workflow, peak_file_suffix)

with open('%s/results/FRP.txt' % (workflow), 'w') as out:
	out.write('Sample\tTotal_fragments\tFragments_in_peaks\tFRP\n')
	for i in range(len(sample_ids)):
		out.write('%s\t%s\t%s\t%s\n' % (sample_ids[i], total_fragments[i], peak_fragments[i], peak_fragments[i]/total_fragments[i]))
	out.close()	
