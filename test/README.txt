Testing in this repository is driven by scripts in the test folder.

Scripts:
- test/run_dryrun_tests.sh
	- Main test harness.
	- Dry-run mode (default): builds DAGs and validates rule wiring without running tools.
	- Integration mode (--integration): runs full workflows via the Snakemake Slurm profile, then validates outputs with test/check_outputs.sh.
	- Uses isolated working directories under .test-work/<workflow>.
	- Cleans .test-work after a successful run unless --keep is specified.

- test/subsample_fastq.sh
	- Creates small, generic FASTQ test inputs (50,000 reads) in test/fastq.
	- Output file names:
		- test_ATAC_R1.fastq.gz / test_ATAC_R2.fastq.gz
		- test_RNAseq_R1.fastq.gz / test_RNAseq_R2.fastq.gz
		- test_Antibody_R1.fastq.gz / test_Antibody_R2.fastq.gz
		- test_Control_R1.fastq.gz / test_Control_R2.fastq.gz
	- If old *_sub.fastq.gz files exist, it renames them to the generic names.
	- If generic files already exist, it skips regenerating them.

- test/check_outputs.sh
	- Output validator used by integration tests.
	- Confirms expected outputs exist and are non-empty.
	- Uses basic sanity checks like line-count thresholds and samtools quickcheck for BAM files.
	- Can also be run manually: bash test/check_outputs.sh <workflow_name> <run_dir>

run_dryrun_tests.sh parameters:
- --integration
	- Run real workflows (not dry-run) using config/profile resources.

- --keep
	- Keep .test-work after successful completion.
	- Default behavior is to delete .test-work on success.

- --workflow NAME
- --workflow=NAME
- -w NAME
	- Run only one workflow from the configured test set.
	- Supported workflows:
		- ATACseq_PE
		- CUT-RUN_PE
		- ChIPseq_PE
		- RNAseq_PE
		- RNAseq_PE_HISAT2_stringtie
		- RNAseq_PE_HISAT2_stringtie_nvltrx

Environment override:
- SNAKEMAKE_CMD
	- Optional override for Snakemake invocation.
	- Example:
		- SNAKEMAKE_CMD="conda run -n CVRCseq snakemake" bash test/run_dryrun_tests.sh

How to run:

1) Prepare small test FASTQs (recommended before integration runs):
source /gpfs/data/cvrcbioinfolab/shared_conda_envs/condaload_CVRCseq.sh
bash test/subsample_fastq.sh

2) Run all dry-run tests:
source /gpfs/data/cvrcbioinfolab/shared_conda_envs/condaload_CVRCseq.sh
bash test/run_dryrun_tests.sh

3) Run one dry-run workflow:
bash test/run_dryrun_tests.sh -w ATACseq_PE

4) Run full integration tests (on a compute node):
source /gpfs/data/cvrcbioinfolab/shared_conda_envs/condaload_CVRCseq.sh
bash test/run_dryrun_tests.sh --integration

5) Run one integration workflow and keep outputs:
bash test/run_dryrun_tests.sh --integration --workflow CUT-RUN_PE --keep

Notes:
- Dry-run mode does not produce pipeline result files.
- Integration mode produces outputs inside .test-work/<workflow>/<workflow>/results.
- On failure, .test-work is preserved for debugging even without --keep.
