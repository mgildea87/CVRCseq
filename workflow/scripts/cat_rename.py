import os, subprocess, sys, yaml
import pandas as pd

#get config file sample_file
with open('config/config.yaml', 'r') as file:
	config = yaml.safe_load(file)
sample_file = config['sample_file']

#Combine fastq files from multiple lanes
def concat():
	cur_dir = sys.argv[2]+'/inputs/fastq/'
	fast_dir = sys.argv[1]
	files = os.listdir(fast_dir)
	samples = {}
	for file in files:
		file = '%s%s' % (fast_dir, file)
		if 'L00' in file:
			if os.path.isfile(file) and file.endswith('fastq.gz'):
					a = file.split('L00')
					f = '%s%s' % (a[0], a[1][2:])
					if f not in samples:
						samples[f] = [file]
					else:
						samples[f].append(file)
		elif os.path.isfile(file) and file.endswith('fastq.gz'):
			samples[file] = [file]
	for sample in samples:
		com = samples[sample]
		com.sort()
		com.insert(0, 'cat')
		with open('%s%s' % (cur_dir, sample.split('/')[-1]), 'w') as R1:			
			subprocess.run(com, stdout=R1)

#rename samples based on sample ids in samples_info.tab. This differs with different workflows
def rename_RNA_SE(sample_file):
	cur_dir = sys.argv[2]+'/inputs/fastq'
	table = pd.read_table(sample_file)
	sample = table['Sample']
	replicate = table['Replicate']
	condition = table['Condition']
	File_R1 = table['File_Name_R1']

	for i in range(len(File_R1)):
		if 'L00' in File_R1[i]:
			file = File_R1[i].split('L00')
			file = '%s%s' % (file[0], file[1][2:])
		else:
			file = File_R1[i]
		if os.path.exists('%s/%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])):
			continue
		if os.path.exists('%s/%s' % (cur_dir,file)):
			os.system('mv %s/%s %s/%s_%s_%s_R1.fastq.gz' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i]))
		elif os.path.exists('%s/%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])) == False:
			print('%s/%s or %s/%s_%s_%s_R1.fastq.gz do not exist!' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i]))
			sys.exit(1)

def rename_RNA_PE(sample_file):
	cur_dir = sys.argv[2]+'/inputs/fastq'
	table = pd.read_table(sample_file)
	sample = table['Sample']
	replicate = table['Replicate']
	condition = table['Condition']
	File_R1 = table['File_Name_R1']
	File_R2 = table['File_Name_R2']

	for i in range(len(File_R1)):
		if 'L00' in File_R1[i]:
			file = File_R1[i].split('L00')
			file = '%s%s' % (file[0], file[1][2:])
		else:
			file = File_R1[i]
		if os.path.exists('%s/%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])):
			continue
		if os.path.exists('%s/%s' % (cur_dir,file)):
			os.system('mv %s/%s %s/%s_%s_%s_R1.fastq.gz' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i]))
		elif os.path.exists('%s/%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])) == False:
			print('%s/%s or %s/%s_%s_%s_R1.fastq.gz do not exist!' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i]))
			sys.exit(1)


	for i in range(len(File_R2)):
		if 'L00' in File_R2[i]:
			file = File_R2[i].split('L00')
			file = '%s%s' % (file[0], file[1][2:])
		else:
			file = File_R2[i]
		if os.path.exists('%s/%s_%s_%s_R2.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])):
			continue
		if os.path.exists('%s/%s' % (cur_dir,file)):
			os.system('mv %s/%s %s/%s_%s_%s_R2.fastq.gz' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i]))
		elif os.path.exists('%s/%s_%s_%s_R2.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])) == False:
			print('%s/%s or %s/%s_%s_%s_R2.fastq.gz do not exist!' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i]))
			sys.exit(1)

def rename_ChIP(sample_file):
	cur_dir = sys.argv[2]+'/inputs/fastq'
	table = pd.read_table(sample_file)
	sample = table['Sample']
	replicate = table['Replicate']
	condition = table['Condition']
	Antibody = table['Antibody']
	File_R1 = table['File_Name_R1']
	File_R2 = table['File_Name_R2']

	for i in range(len(File_R1)):
		if 'L00' in File_R1[i]:
			file = File_R1[i].split('L00')
			file = '%s%s' % (file[0], file[1][2:])
		else:
			file = File_R1[i]
		if os.path.exists('%s/%s_%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i], Antibody[i])):
			continue
		if os.path.exists('%s/%s' % (cur_dir,file)):
			os.system('mv %s/%s %s/%s_%s_%s_%s_R1.fastq.gz' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i], Antibody[i]))
		elif os.path.exists('%s/%s_%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i], Antibody[i])) == False:
			print('%s/%s or %s/%s_%s_%s_%s_R1.fastq.gz do not exist!' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i], Antibody[i]))
			sys.exit(1)

	for i in range(len(File_R2)):
		if 'L00' in File_R2[i]:
			file = File_R2[i].split('L00')
			file = '%s%s' % (file[0], file[1][2:])
		else:
			file = File_R2[i]
		if os.path.exists('%s/%s_%s_%s_%s_R2.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i], Antibody[i])):
			continue
		if os.path.exists('%s/%s' % (cur_dir,file)):
			os.system('mv %s/%s %s/%s_%s_%s_%s_R2.fastq.gz' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i], Antibody[i]))
		elif os.path.exists('%s/%s_%s_%s_%s_R2.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i], Antibody[i])) == False:
			print('%s/%s or %s/%s_%s_%s_%s_R2.fastq.gz do not exist!' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i], Antibody[i]))
			sys.exit(1)

concat()

if sys.argv[2] == 'CUT-RUN_PE' or sys.argv[2] == 'ChIPseq_PE':
	rename_ChIP(sample_file)
if sys.argv[2] == 'RNAseq_PE' or sys.argv[2] == 'RNAseq_PE_HISAT2_stringtie' or sys.argv[2] == 'RNAseq_PE_HISAT2_stringtie_nvltrx' or sys.argv[2] == 'ATACseq_PE':
	rename_RNA_PE(sample_file)
if sys.argv[2] == 'RNAseq_SE':
	rename_RNA_SE(sample_file)