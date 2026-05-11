# CVRCseq Container Guide (Docker + Singularity)

This project uses a Docker-built software environment that is converted to Singularity on HPC.

## Why this setup

- The CVRCseq software stack is large and difficult to rebuild on shared HPC systems without elevated permissions.
- Direct `singularity build --fakeroot` is not always available to all users.
- Building with GitHub Actions + Docker Hub makes the container reproducible and easy to update.

## What is in the container

- The full CVRCseq conda software environment from `workflow/envs/CVRCseq.yml`.
- Tools such as Snakemake, Python, MultiQC, and workflow dependencies.

## What is NOT in the container

- Your project data.
- Your active run directory.
- Workflow execution state.

The workflow code and data stay on GPFS and are mounted at runtime with `--bind /gpfs`.

## Build architecture

1. You update and commit `workflow/envs/CVRCseq.yml` and/or `Dockerfile`.
2. GitHub Actions workflow `.github/workflows/build-docker.yml` builds and pushes `docker://mgildea87/cvrcsseq:latest`.
3. On HPC, pull and convert to Singularity `.sif` using `singularity pull`.

## One-time setup

1. Docker Hub account and repository `mgildea87/cvrcsseq`.
2. GitHub repository secrets:
   - `DOCKERHUB_USERNAME`
   - `DOCKERHUB_TOKEN`

## Automatic rebuild triggers

The workflow rebuilds the Docker image on pushes to `main` or `development` when either of these files changes:

- `workflow/envs/CVRCseq.yml`
- `Dockerfile`

## Update process

1. Update software in the shared conda environment.
2. Export and clean the env file:

```bash
conda env export -p /gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq > workflow/envs/CVRCseq.yml
sed -i '/^prefix:/d' workflow/envs/CVRCseq.yml
```

3. Commit and push:

```bash
git add workflow/envs/CVRCseq.yml
git commit -m "Update CVRCseq container environment"
git push origin development
```

4. Wait for GitHub Actions build completion.
5. Pull the updated container on HPC:

```bash
module load singularity/3.11.5
singularity pull --dir /gpfs/data/cvrcbioinfolab/shared_conda_envs/ docker://mgildea87/cvrcsseq:latest
mv /gpfs/data/cvrcbioinfolab/shared_conda_envs/cvrcsseq_latest.sif /gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq.sif
```

## Validate pulled container

```bash
singularity exec --bind /gpfs /gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq.sif snakemake --version
singularity exec --bind /gpfs /gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq.sif python --version
singularity exec --bind /gpfs /gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq.sif multiqc --version
```

## How runtime execution works

- `singularity exec` runs tools from inside the container.
- `--bind /gpfs` exposes your repo, configs, and data inside the container.
- The pipeline reads and writes to your GPFS paths as usual.

`snakemake_init.sh` now defaults to container execution if a default image exists at `/gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq.sif`.

Use `-i` to override the image path explicitly:

```bash
bash workflow/scripts/snakemake_init.sh -w RNAseq_PE -d /path/to/fastq -i /gpfs/data/cvrcbioinfolab/shared_conda_envs/CVRCseq.sif
```

You can also override the default image path with an environment variable:

```bash
export CVRCSEQ_SIF=/path/to/other_image.sif
bash workflow/scripts/snakemake_init.sh -w RNAseq_PE -d /path/to/fastq
```

If no default image exists and `-i` is not set, the workflow falls back to host conda execution.

## Singularity module behavior

- The Singularity module is required only when running with `-i` (container mode).
- `workflow/scripts/snakemake_init.sh` automatically checks for `singularity` when `-i` is provided.
- If `singularity` is not already on PATH, the script attempts `module load singularity/3.11.5`.
- If `singularity` is still unavailable after that, the script exits with an error.

If you prefer manual setup, load the module yourself before running:

```bash
module load singularity/3.11.5
```

## Troubleshooting

- Build error: `invalid reference format`
  - Cause: Docker image names cannot contain spaces.
  - Fix: use `mgildea87/cvrcsseq:latest`.

- Build does not reflect latest packages
  - Cause: updated conda env was not exported/committed.
  - Fix: re-export `workflow/envs/CVRCseq.yml`, commit, and push again.

- HPC cannot use `--fakeroot`
  - Expected in this setup. Use Docker Hub build + `singularity pull`.

- `singularity pull` creates `cvrcsseq_latest.sif`
  - This is normal. Rename to `CVRCseq.sif` for consistent usage.
