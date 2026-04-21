CVRCseq Test Suite
==================

Overview
--------
Three levels of testing are available:

  Unit tests      Fast, no dependencies. Test individual scripts in isolation.
  Dry-run tests   Validate Snakemake DAG/rule wiring without running any tools.
  Integration     Run full pipelines end-to-end on a compute node via SLURM.


Quick Start
-----------
# Activate the conda environment (required for dry-run and integration tests)
source /gpfs/data/cvrcbioinfolab/shared_conda_envs/condaload_CVRCseq.sh

# Run all unit tests (no conda required)
bash test/test_snakemake_init.sh
python -m pytest test/test_cat_rename.py -v

# Run all dry-run tests
bash test/run_dryrun_tests.sh

# Run full integration tests (on a compute node)
bash test/run_dryrun_tests.sh --integration


Unit Tests
----------
test/test_snakemake_init.sh
  Tests argument validation in workflow/scripts/snakemake_init.sh.
  Stubs out conda, snakemake, and multiqc — no pipeline dependencies required.
  Run with: bash test/test_snakemake_init.sh

  Tests (10):
    - exits with code 1 when no arguments provided
    - exits with code 1 when -d (fastq directory) is missing
    - exits with code 1 when -w (workflow) is missing
    - exits with code 1 for an unsupported workflow name
    - passes validation for each supported workflow:
        RNAseq_PE, ATACseq_PE, CUT-RUN_PE, ChIPseq_PE, RNAseqTE_PE
    - -c flag is boolean and does not consume the -d value

test/test_cat_rename.py
  Tests all functions in workflow/scripts/cat_rename.py using pytest.
  Uses tempfile.TemporaryDirectory and unittest.mock.patch — no real FASTQ data required.
  Run with: python -m pytest test/test_cat_rename.py -v

  Tests (12):
    TestConcat
      - test_single_lane_passthrough:  file with no L00 tag is passed through as-is
      - test_multi_lane_merging:       L001 + L002 files are merged into one output
      - test_non_fastq_files_ignored:  non-.fastq.gz files are skipped
    TestRenameRNASE
      - test_renames_correctly:        file copied to Sample_Condition_Replicate_R1.fastq.gz
      - test_skips_if_output_exists:   copy skipped if renamed file already exists
      - test_exits_on_missing_file:    sys.exit(1) when neither source nor renamed file found
      - test_removes_concat_files:     intermediate concat files removed after renaming
    TestRenameRNAPE
      - test_renames_r1_and_r2:        both R1 and R2 renamed correctly
      - test_exits_on_missing_r2:      sys.exit(1) when R2 is absent
    TestRenameChIP
      - test_renames_with_antibody:    Antibody column included in output filename
    TestMainArgValidation
      - test_exits_with_no_args:       main() exits with code 1 when no args provided
      - test_exits_with_one_arg:       main() exits with code 1 when only one arg provided


Dry-Run and Integration Tests
------------------------------
test/run_dryrun_tests.sh
  Main test harness for pipeline-level testing.
  Dry-run mode (default): builds Snakemake DAGs and validates rule wiring without running tools.
  Integration mode:       runs full workflows on SLURM, then validates outputs with check_outputs.sh.
  Uses isolated working directories under .test-work/<workflow>.
  Deletes .test-work after a successful run unless --keep is specified.

  Parameters:
    --integration         Run real workflows instead of dry-run.
    --keep                Keep .test-work after successful completion.
    --workflow NAME        Run only one workflow. Also accepted as --workflow=NAME or -w NAME.
    SNAKEMAKE_CMD=...     Environment variable to override the snakemake invocation.
                          Example: SNAKEMAKE_CMD="conda run -n CVRCseq snakemake" bash test/run_dryrun_tests.sh

  Supported workflows:
    ATACseq_PE
    ChIPseq_PE
    CUT-RUN_PE
    RNAseq_PE
    RNAseq_SE
    RNAseqTE_PE
    RNAseq_PE_HISAT2_stringtie
    RNAseq_PE_HISAT2_stringtie_nvltrx
    sRNAseq_SE

  Examples:
    bash test/run_dryrun_tests.sh                                     # all dry-runs
    bash test/run_dryrun_tests.sh -w ATACseq_PE                       # one dry-run
    bash test/run_dryrun_tests.sh --integration                       # all workflows
    bash test/run_dryrun_tests.sh --integration --workflow CUT-RUN_PE --keep

test/subsample_fastq.sh
  Generates small FASTQ test inputs (50,000 reads) in test/fastq/ for integration runs.
  Skips regeneration if output files already exist.
  Output files:
    test_ATAC_R1.fastq.gz / test_ATAC_R2.fastq.gz
    test_RNAseq_R1.fastq.gz / test_RNAseq_R2.fastq.gz
    test_Antibody_R1.fastq.gz / test_Antibody_R2.fastq.gz
    test_Control_R1.fastq.gz / test_Control_R2.fastq.gz

test/check_outputs.sh
  Output validator called automatically by integration tests.
  Can also be run manually: bash test/check_outputs.sh <workflow_name> <run_dir>

  Internal checks:
    check_exists_nonempty   File exists and has non-zero size (bigwigs, SEACR beds).
    check_min_lines         File exists and has at least N lines (count tables, peak files, FRiP summaries).
    check_bam               samtools quickcheck — valid BAM header and EOF block.

  Per-workflow output checks:
    ATACseq_PE                        filtered_sorted BAMs valid; dedup_filtered_sorted bigwigs non-empty;
                                      narrowPeak >= 1 line; FRiP summary >= 2 lines
    ChIPseq_PE                        sorted BAMs valid; bigwigs non-empty;
                                      narrowPeak >= 1 line; FRiP summary >= 2 lines
    CUT-RUN_PE                        sorted BAMs valid; SEACR stringent.bed non-empty;
                                      MACS2 broadPeak >= 1 line; MACS2 and SEACR FRiP summaries >= 2 lines
    RNAseq_PE / RNAseq_SE / sRNAseq_SE   count_table.txt >= 2 lines; alignment BAMs valid
    RNAseqTE_PE                       TEcount/count_table_all.csv >= 2 lines; alignment BAMs valid
    RNAseq_PE_HISAT2_stringtie        gene and transcript count matrices >= 2 lines; alignment BAMs valid
    RNAseq_PE_HISAT2_stringtie_nvltrx same as above but paths are under stringtie/merged/


Notes
-----
- Dry-run mode does not produce pipeline result files.
- Integration mode produces outputs inside .test-work/<workflow>/<workflow>/results.
- On failure, .test-work is preserved for debugging even without --keep.
- RNAseqTE_PE integration requires a valid TE annotation in test/config/config_RNAseqTE_PE.yaml (TE_GTF).
