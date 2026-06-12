# CVRCseq

CVRCseq is a unified Snakemake workflow collection for common NGS analyses on Slurm-based HPC systems (developed for NYU UltraViolet).

## Available Workflows

### RNA-seq

- `RNAseq_PE`: paired-end, `fastqc -> fastp -> STAR -> featureCounts`
- `RNAseq_SE`: single-end, `fastqc -> fastp -> STAR -> featureCounts`
- `RNAseq_PE_HISAT2_stringtie`: paired-end, `fastqc -> fastp -> HISAT2 -> StringTie`
- `RNAseq_PE_HISAT2_stringtie_nvltrx`: paired-end, `fastqc -> fastp -> HISAT2 -> StringTie -> novel transcript workflow`
- `RNAseqTE_PE`: paired-end, `fastqc -> fastp -> STAR -> TEcount`

### Small RNA-seq

- `sRNAseq_SE`: single-end, `fastqc -> umi-tools -> STAR -> featureCounts`

### DNA Binding / Enrichment

- `ChIPseq_PE`: paired-end, `fastqc -> fastp -> bowtie2 -> MACS2`
- `CUT-RUN_PE`: paired-end, `fastqc -> fastp -> bowtie2 -> MACS2`
- `ATACseq_PE`: paired-end, `fastqc -> fastp -> bowtie2 -> MACS2`

## Repository Structure

- `workflow/Snakefile`: top-level workflow entry point; loads one rules file based on `workflow` in config.
- `workflow/rules/*.smk`: per-workflow rule definitions.
- `workflow/scripts/snakemake_init.sh`: main launcher script.
- `workflow/scripts/cat_rename.py`: optional preprocessing step for lane concatenation and FASTQ renaming.
- `config/config.yaml`: global and workflow-specific parameters.
- `config/samples_info.tab`: sample metadata table.
- `config/profile/config.yaml`: Snakemake profile and Slurm defaults.
- `workflow/envs/CVRCseq.yml`: conda environment definition.

## Configuration

### Sample Metadata (`config/samples_info.tab`)

Expected columns include:

1. FASTQ file names (R1/R2)
2. User-friendly sample name
3. Condition
4. Replicate
5. Antibody/control label (required for ChIP-seq and CUT-RUN)
6. Final sample ID (used for renamed FASTQ output)
7. Optional additional metadata

Notes:

- `cat_rename.py` concatenates multi-lane FASTQs and renames files from this table.
- For ChIP-seq and CUT-RUN pairs, keep sample name/condition/replicate consistent between IP and control rows.

### Main Config (`config/config.yaml`)

Common keys:

- `sample_file`: path to sample table (default `config/samples_info.tab`)
- `workflow`: active workflow name (set automatically by `snakemake_init.sh`)
- `genome`: index path (STAR, HISAT2, or bowtie2 depending on workflow)
- `GTF`: annotation file path

Workflow-specific keys:

- `CUT-RUN_PE`:
  - `spike_genome`
  - `chromosome_lengths`
  - `effective_genome_size`
- `ChIPseq_PE`, `ATACseq_PE`:
  - `effective_genome_size`
- `RNAseq_PE_HISAT2_stringtie`, `RNAseq_PE_HISAT2_stringtie_nvltrx`:
  - `prepDE_length`
  - `stringtie_strandedness` (example: `"--rf"`)
- `RNAseqTE_PE`:
  - `TE_GTF`
  - `TE_strandedness` (example: `"reverse"`)
- `RNAseq_PE`, `RNAseq_SE`, `sRNAseq_SE`:
  - `featurecounts_strandedness` (`0`, `1`, or `2`)

## Running the Pipeline

### 1) Clone

```bash
git clone https://github.com/mgildea87/CVRCseq.git
cd CVRCseq
```

### 2) Prepare inputs

- Update `config/samples_info.tab`.
- Update `config/config.yaml` for your references and workflow settings.

### 3) Launch

```bash
bash workflow/scripts/snakemake_init.sh -d /path/to/fastq -w RNAseq_PE
```

Options:

- `-h`: help
- `-d`: FASTQ directory (required)
- `-w`: workflow name (required)
- `-s`: extra Snakemake args (quote multiple flags, for example `-s "--dryrun --quiet"`)
- `-c`: skip `cat_rename.py`
- `-i`: override Singularity image path

If needed, unlock a stale Snakemake directory:

```bash
snakemake --unlock --profile config/profile
```

This requires loading the container or conda evironment where snakemake is installed

## Execution Mode (Container vs Host)

Default behavior:

- Uses Singularity image at `/gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq.sif` if available.
- Falls back to host conda environment (`/gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq`) if the image is absent and `-i` is not provided.

Pull the image manually:

```bash
module load singularity/3.11.5
singularity pull --dir /gpfs/data/cvrcbioinfolab/shared_conda_envs/ docker://mgildea87/cvrcsseq:latest
```

For additional container details, see [container/README.md](container/README.md).

## Host Orchestrator Requirements (Option 1)

In the default architecture, Snakemake runs on the host as the workflow
orchestrator, and each rule executes in the Singularity container via
`--use-singularity`.

This means one host-side Snakemake installation is still required even in
container mode.

### How `snakemake_init.sh` selects Snakemake

`workflow/scripts/snakemake_init.sh` resolves the Snakemake executable in this
order:

1. `CVRCSEQ_SNAKEMAKE_BIN` (explicit override).
2. `${CVRCSEQ_HOST_ENV:-/gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq}/bin/snakemake` (preferred default).
3. A runnable `snakemake` found on `PATH`.
4. Exit with an error if none of the above are valid.

### Recommended host versions

To minimize orchestration drift, keep host Snakemake and Python aligned with
the CVRCseq environment definition in `workflow/envs/CVRCseq.yml`:

- Python: `3.10.2`
- Snakemake: `7.21.0`

Note: exact host/container version matching is not strictly required for tools
executed inside rules, but Snakemake major/minor compatibility on the host is
strongly recommended.

## Running on a Different System

If you are running outside NYU UltraViolet, do not rely on default `/gpfs/...`
paths. Set explicit host and container paths before launching.

### Minimum requirements

- A runnable host Snakemake installation (recommended: Snakemake `7.21.0` with Python `3.10.x`)
- Singularity/Apptainer available on the host
- A CVRCseq container image (`.sif`) accessible on your filesystem
- A Slurm environment, or a compatible profile if adapting to a different scheduler

### Recommended portable launch pattern

Set explicit paths so `snakemake_init.sh` does not depend on site-specific defaults:

```bash
export CVRCSEQ_SNAKEMAKE_BIN=/path/to/host/env/bin/snakemake
export CVRCSEQ_HOST_ENV=/path/to/host/env
export CVRCSEQ_SIF=/path/to/CVRCseq.sif

bash workflow/scripts/snakemake_init.sh -d /path/to/fastq -w RNAseq_PE
```

Optional: you can still pass `-i /path/to/CVRCseq.sif` on the command line to
override image selection for a single run.

### How host Snakemake is discovered

With the launcher logic in `workflow/scripts/snakemake_init.sh`, executable
resolution is:

1. `CVRCSEQ_SNAKEMAKE_BIN`
2. `${CVRCSEQ_HOST_ENV}/bin/snakemake` (or default host env path)
3. `snakemake` from `PATH`

### Portability note

Current container wrapping binds `/gpfs` into the container. On non-GPFS
systems, you may need to adjust bind paths in the launcher/profile to match
your local filesystem layout.

## Running on a Compute Node

Launching from a compute node is recommended. Update `workflow/scripts/launch_sbatch.sh` and submit:

```bash
sbatch workflow/scripts/launch_sbatch.sh
```

## Tool Links

- [Snakemake](https://snakemake.github.io/)
- [STAR](https://github.com/alexdobin/STAR)
- [FastQC](https://github.com/s-andrews/FastQC)
- [fastp](https://github.com/OpenGene/fastp)
- [subread/featureCounts](https://github.com/ShiLab-Bioinformatics/subread)
- [HISAT2](https://daehwankimlab.github.io/hisat2/)
- [StringTie](https://ccb.jhu.edu/software/stringtie/)
- [TEtranscripts/TEcount](https://github.com/mhammell-laboratory/TEtranscripts)
- [UMI-tools](https://github.com/CGATOxford/UMI-tools)
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [MACS2](https://pypi.org/project/MACS2/)
















