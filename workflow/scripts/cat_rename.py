import os, shutil, subprocess, sys, yaml
import pandas as pd

def copy_file(src, dst):
	shutil.copy2(src, dst)

def remove_file(path):
	os.remove(path)

def main():
	if len(sys.argv) < 3:
		print('Usage: python cat_rename.py <fastq_dir> <workflow>')
		sys.exit(1)
	print('copying and concatenating .fastq files')
	#get config file sample_file
	with open('config/config.yaml', 'r') as file:
		config = yaml.safe_load(file)
	sample_file = config['sample_file']
	sample_table = pd.read_table(sample_file)
	
	concat_file_names = concat(sample_table)
	workflow = sys.argv[2]
	
	if workflow == 'CUT-RUN_PE' or workflow == 'ChIPseq_PE':
		rename_ChIP(sample_table, concat_file_names)
	elif workflow == 'RNAseq_PE' or workflow == 'RNAseq_PE_HISAT2_stringtie' or workflow == 'RNAseq_PE_HISAT2_stringtie_nvltrx' or workflow == 'ATACseq_PE' or workflow == 'RNAseqTE_PE':
		rename_RNA_PE(sample_table, concat_file_names)
	elif workflow == 'RNAseq_SE' or workflow == 'sRNAseq_SE':
		rename_RNA_SE(sample_table, concat_file_names)
	else:
		print('Unsupported workflow: %s' % workflow)
		sys.exit(1)


#Combine fastq files from multiple lanes
def concat(sample_table):
	cur_dir = sys.argv[2]+'/inputs/fastq/'
	fast_dir = sys.argv[1]
	files = os.listdir(fast_dir)
	samples = {}
	concat_file_names = []
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
			subprocess.run(com, stdout=R1, check=True)
			concat_file_names.append('%s%s' % (cur_dir, sample.split('/')[-1]))
	return(concat_file_names)

def rename_RNA_SE(sample_table, concat_file_names):
	cur_dir = sys.argv[2]+'/inputs/fastq'
	sample = sample_table['Sample']
	replicate = sample_table['Replicate']
	condition = sample_table['Condition']
	File_R1 = sample_table['File_Name_R1']
	for i in range(len(File_R1)):
		if 'L00' in File_R1[i]:
			file = File_R1[i].split('L00')
			file = '%s%s' % (file[0], file[1][2:])
		else:
			file = File_R1[i]
		if os.path.exists('%s/%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])):
			continue
		if os.path.exists('%s/%s' % (cur_dir,file)):
			copy_file('%s/%s' % (cur_dir,file), '%s/%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i]))
		elif os.path.exists('%s/%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])) == False:
			print('%s/%s or %s/%s_%s_%s_R1.fastq.gz do not exist!' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i]))
			sys.exit(1)
	# Remove unused files
	for file in concat_file_names:
		if os.path.exists(file):
			remove_file(file)

def rename_RNA_PE(sample_table, concat_file_names):
	cur_dir = sys.argv[2]+'/inputs/fastq'
	sample = sample_table['Sample']
	replicate = sample_table['Replicate']
	condition = sample_table['Condition']
	File_R1 = sample_table['File_Name_R1']
	File_R2 = sample_table['File_Name_R2']
	for i in range(len(File_R1)):
		if 'L00' in File_R1[i]:
			file = File_R1[i].split('L00')
			file = '%s%s' % (file[0], file[1][2:])
		else:
			file = File_R1[i]
		if os.path.exists('%s/%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])):
			continue
		if os.path.exists('%s/%s' % (cur_dir,file)):
			copy_file('%s/%s' % (cur_dir,file), '%s/%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i]))
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
			copy_file('%s/%s' % (cur_dir,file), '%s/%s_%s_%s_R2.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i]))
		elif os.path.exists('%s/%s_%s_%s_R2.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i])) == False:
			print('%s/%s or %s/%s_%s_%s_R2.fastq.gz do not exist!' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i]))
			sys.exit(1)
	# Remove unused files
	for file in concat_file_names:
		if os.path.exists(file):
			remove_file(file)

def rename_ChIP(sample_table, concat_file_names):
	cur_dir = sys.argv[2]+'/inputs/fastq'
	sample = sample_table['Sample']
	replicate = sample_table['Replicate']
	condition = sample_table['Condition']
	Antibody = sample_table['Antibody']
	File_R1 = sample_table['File_Name_R1']
	File_R2 = sample_table['File_Name_R2']
	for i in range(len(File_R1)):
		if 'L00' in File_R1[i]:
			file = File_R1[i].split('L00')
			file = '%s%s' % (file[0], file[1][2:])
		else:
			file = File_R1[i]
		if os.path.exists('%s/%s_%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i], Antibody[i])):
			continue
		if os.path.exists('%s/%s' % (cur_dir,file)):
			copy_file('%s/%s' % (cur_dir,file), '%s/%s_%s_%s_%s_R1.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i], Antibody[i]))
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
			copy_file('%s/%s' % (cur_dir,file), '%s/%s_%s_%s_%s_R2.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i], Antibody[i]))
		elif os.path.exists('%s/%s_%s_%s_%s_R2.fastq.gz' % (cur_dir,sample[i],condition[i],replicate[i], Antibody[i])) == False:
			print('%s/%s or %s/%s_%s_%s_%s_R2.fastq.gz do not exist!' % (cur_dir,file,cur_dir,sample[i],condition[i],replicate[i], Antibody[i]))
			sys.exit(1)
	# Remove unused files
	for file in concat_file_names:
		if os.path.exists(file):
			remove_file(file)

if __name__ == "__main__":
    main()
