import subprocess, sys, pandas as pd, yaml
with open('config/config.yaml', 'r') as file:
	config = yaml.safe_load(file)
workflow = config['workflow']
sample_file = config['sample_file']
sample = pd.read_table(sample_file)['Sample']
replicate = pd.read_table(sample_file)['Replicate']
condition = pd.read_table(sample_file)['Condition']
Antibody = pd.read_table(sample_file)['Antibody']

#peak files get named different things with macs vs seacr will need to update if peak caller is changed or other workflows are added
if workflow == 'ChIPseq':
	peak_file_suffix = '_peaks.narrowPeak'
if workflow == 'CUT-RUN':
	peak_file_suffix = '.stringent.bed'

print(peak_file_suffix)

sample_ids = []
for i in range(len(sample)):
	sample_ids.append('%s_%s_%s' % (sample[i], condition[i], replicate[i]))
sample_ids = pd.unique(sample_ids).tolist()

sample_ids_file = []
for i in range(len(sample)):
	sample_ids_file.append('%s_%s_%s_%s' % (sample[i], condition[i], replicate[i], Antibody[i]))

total_fragments = []
peak_fragments = []

def FRP(sample, workflow, peak_file_suffix):
	peak_id = sample.split('_')[0:-1]
	peak_id = "_".join(peak_id)
	print(sample)
	cmd_1 = ['bedtools', 'bamtobed', '-bedpe', '-i', '%s/results/alignment/%s.bam' % (workflow, sample)] 
	cmd_2 = ['bedtools', 'intersect', '-a', 'stdin', '-b', '%s/results/peaks/%s%s' % (workflow, peak_id, peak_file_suffix)]
	cmd_3 = ['wc', '-l']
	step_1 = subprocess.Popen(cmd_1, stdout=subprocess.PIPE)
	step_2 = subprocess.Popen(cmd_2, stdin = step_1.stdout, stdout=subprocess.PIPE)
	step_3 = subprocess.Popen(cmd_3, stdin = step_2.stdout, stdout=subprocess.PIPE)
	frag_peaks = step_3.stdout.read()
	frag_peaks.strip()
	peak_fragments.append(float(frag_peaks))

	cmd_4 = ['samtools', 'view', '%s/results/alignment/%s.bam' % (workflow, sample)]
	step_1 = subprocess.Popen(cmd_4, stdout=subprocess.PIPE)
	step_2 = subprocess.Popen(cmd_3, stdin = step_1.stdout, stdout=subprocess.PIPE)

	frag = step_2.stdout.read()
	frag.strip()
	frag = float(frag)
	total_fragments.append(frag/2)

for samples in sample_ids_file:
	FRP(samples, workflow, peak_file_suffix)

with open('%s/results/FRP.txt' % (workflow), 'w') as out:
	out.write('Sample\tTotal_fragments\tFragments_in_peaks\tFRP\n')
	for i in range(len(sample_ids_file)):
		out.write('%s\t%s\t%s\t%s\n' % (sample_ids_file[i], total_fragments[i], peak_fragments[i], peak_fragments[i]/total_fragments[i]))
	out.close()	
