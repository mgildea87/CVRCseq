# 🐍 **CVRCseq** 🐍

CVRCseq is a collection of commonly used pipelines integrated into a single workflow via **Snakemake**. Previously, these existed as individual Snakemake workflows. This unified workflow is designed to run on **NYU's UltraViolet HPC**, which utilizes **Slurm** and offers a variety of different node types.

---

## **Available Pipelines**

### 🧬 **RNA-seq Analysis**
Currently, there are **5 RNA-seq analysis pipelines** available:

1. **RNAseq_PE**  
   *Paired-end data*: `fastqc → fastp → STAR → featurecounts`

2. **RNAseq_SE**  
   *Single-end data*: `fastqc → fastp → STAR → featurecounts`

3. **RNAseq_HISAT2_stringtie**  
   *Paired-end data*: `fastqc → fastp → HISAT2 → stringtie`

4. **RNAseq_HISAT2_stringtie_nvltrx**  
   *Paired-end data*: `fastqc → fastp → HISAT2 → stringtie → novel transcript identification`

5. **RNAseqTE_PE**  
   *Paired-end data*: `fastqc → fastp → STAR → TEcount`

---

### 🔬🧬 **Small RNA-seq Analysis**
Currently, there is **1 small RNA-seq analysis pipeline** available, designed to work with the **QIAseq miRNA Library Kit** from Qiagen:

1. **sRNAseq_SE**  
   *Single-end data*: `fastqc → umi-tools → STAR → featurecounts`

---

### 🧬 **DNA Binding/Enrichment Analysis**
Currently, there are **3 DNA binding/enrichment pipelines** available:

1. **ChIPseq_PE**  
   *Paired-end data*: `fastqc → fastp → bowtie2 → macs2`

2. **CUT-RUN_PE**  
   *Paired-end data*: `fastqc → fastp → bowtie2 → seacr & macs2`

3. **ATACseq_PE**  
   *Paired-end data*: `fastqc → fastp → bowtie2 → macs2`


---

## **📁 File Descriptions**

### **Snakefiles**
- **`workflow/Snakefile`** - Launches individual pipelines located in `workflow/rules`

---

### **Configuration Files**

#### **`config/samples_info.tab`**
Tab-delimited file containing sample metadata:

1. **R1 and R2 fastq file names** - As received from sequencing center (remove lane numbers `L00X` for multi-lane samples)
2. **Simple sample names** - User-defined identifiers
3. **Condition** - Experimental condition (e.g., diabetic vs non_diabetic)
4. **Replicate number** - Biological replicate identifier
5. **Antibody column** - Required for ChIPseq/CUT-RUN (specifies antibody vs control samples)
6. **Final sample ID** - Concatenation of sample name, condition, replicate, and antibody columns
7. **Additional metadata** - Can be added for downstream analysis

> **Notes:** 

1. `cat_rename.py` handles concatenation of multi-lane fastq files and renaming based on this table.
2. **Paired samples** - For ChIPseq/CUT-RUN, sample name/condition/replicate should be identical between antibody and control pairs

---

#### **`config/config.yaml`**
Contains general and workflow-specific configuration parameters:

##### **Generic Requirements:**
- **`sample_file`** - Location of `samples_info.tab` (default: `config/samples_info.tab`)
- **`workflow`** - Name of workflow being used
- **`genome`** - Location of indexed genome:
  - RNAseq_PE/RNAseq_SE/sRNAseq_SE: STAR 2.7.7a index
  - HISAT2 workflows: HISAT2 index
  - ChIPseq/CUT-RUN/ATACseq: bowtie2 index
- **`GTF`** - Location of annotation file

##### **Workflow-Specific config settings:**

**CUT-RUN_PE:**
- `spike_genome` - Spike-in genome index (bowtie2)
- `chromosome_lengths` - Required for spike-in normalization. This file can be found in the STAR genome index folder (chrLength.txt)
- `effective_genome_size` - For MACS2

**ChIPseq_PE & ATACseq_PE:**
- `effective_genome_size` - For MACS2

**RNAseq_HISAT2_stringtie variants:**
- `prepDE_length` - Average fragment length for stringtie prepDE script

**RNAseqTE_PE:**
- `TE_GTF` - GTF file with TE annotations (available from [MGH lab](https://www.dropbox.com/scl/fo/jdpgn6fl8ngd3th3zebap/ACdZkShDC1au-OckIipI5kM/TEtranscripts/TE_GTF?rlkey=41oz6ppggy82uha5i3yo1rnlx&e=1&subfolder_nav_tracking=1&dl=0))

---

#### **`config/profile/config.yaml`**
Defines default **Slurm resources** for each rule.

---

### 📝 **Scripts**

#### **`workflow/scripts/cat_rename.py`**
Preprocessing script that:
1. Concatenates fastq files split across multiple sequencing lanes
2. Renames fastq files from verbose sequencing center IDs to user-defined names
3. Creates new files as `sample_id_Rx.fastq.gz`
4. Executed automatically via `snakemake_init.sh`

> **Skip option:** Use `-c` flag with `snakemake_init.sh` to bypass this step.

---

#### **`workflow/scripts/snakemake_init.sh`**
Main execution script that:
1. Executes `cat_rename.py`
2. Loads conda environment
3. Launches Snakemake pipeline
4. Runs MultiQC for quality control

---

#### **`workflow/scripts/launch_sbatch.sh`**
Launches pipeline from compute node (recommended over login node). Edit the `snakemake_init.sh` command with desired parameters and submit via `sbatch`.

---

#### **`workflow/scripts/condaload_CVRCseq.sh`**
Sets environment variables and loads the conda environment.

---

#### **`workflow/scripts/FRP.py`**
Computes **Fraction of Reads in Peaks (FRP)** and outputs a summary table with:
- FRP values
- Total fragments
- Fragments within peaks

---

#### **`workflow/scripts/combine_TE_counts.py`**
combines counts from TEcount into a single .csv file.

---

### **Environment**
#### **`workflow/envs/CVRCseq.yml`**
Contains conda environment specifications for the pipeline.
 
---

## **🚀 Usage Instructions**

### **Getting Started:**
1. **Clone repository:**
  git clone https://github.com/mgildea87/CVRCseq.git
2. **Update sample information**
  Edit config/samples_info.tab with fastq.gz file names and desired sample, condition, replicate names, and Antibody/IgG control status (if using)
3. **Configure workflow**
  Update config.yaml with project-specific settings
4. **Customize parameters (optional)** 
  Set workflow specific parameters in the appropriate worklow/rules .smk file if desired. e.g. alignment parameters. 
5. **Launch pipeline**
  bash workflow/scripts/snakemake_init.sh
    Description of parameters:
      -h	help"
      -d	.fastq directory"
      -s	parameters to pass to snakemake (e.g. --unlock)
      -w	workflow name (e.g. 'RNAseq_PE')
      -c	Skip cat_rename.py. Use to skip copying, concatenating, and renaming of .fastq files to the *workflow*/inputs/fastq/ local directory

## Software links

[snakemake](https://snakemake.github.io/), 
[STAR](https://github.com/alexdobin/STAR), 
[fastqc](https://github.com/s-andrews/FastQC), 
[fastp](https://github.com/OpenGene/fastp), 
[subread - featurecounts](https://github.com/ShiLab-Bioinformatics/subread), 
[HISAT2](https://daehwankimlab.github.io/hisat2/), 
[stringtie](https://ccb.jhu.edu/software/stringtie/), 
[TEcount](https://github.com/mhammell-laboratory/TEtranscripts), 
[umi-tools](https://github.com/CGATOxford/UMI-tools), 
[bowtie2](https://github.com/BenLangmead/bowtie2), 
[macs2](https://pypi.org/project/MACS2/), 
[seacr](https://github.com/FredHutch/SEACR)

















