#!/bin/bash
# Unit tests for workflow/scripts/snakemake_init.sh argument validation.
# Tests only the argument parsing and validation layer — does not invoke
# conda, snakemake, or multiqc. Safe to run on a login node without any
# pipeline dependencies.

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$REPO_ROOT/workflow/scripts/snakemake_init.sh"
PASS=0
FAIL=0

# ---------------------------------------------------------------------------
# Setup: create a temp working directory with stubs so the script can pass
# argument validation without needing conda, snakemake, or multiqc.
# ---------------------------------------------------------------------------
STUB_DIR=$(mktemp -d)
WORK_DIR=$(mktemp -d)

cleanup() { rm -rf "$STUB_DIR" "$WORK_DIR"; }
trap cleanup EXIT

# Stub executables for post-validation commands
for cmd in conda snakemake multiqc; do
    printf '#!/bin/bash\nexit 0\n' > "$STUB_DIR/$cmd"
    chmod +x "$STUB_DIR/$cmd"
done

# Stub singularity to execute container commands against local stubs
cat > "$STUB_DIR/singularity" << 'EOF'
#!/bin/bash
if [[ "$1" == "exec" ]]; then
    shift
    while [[ $# -gt 0 ]]; do
        if [[ "$1" == "--bind" ]]; then
            shift 2
            continue
        fi
        if [[ "$1" == *.sif ]]; then
            shift
            continue
        fi
        break
    done
    "$@"
    exit $?
fi
exit 0
EOF
chmod +x "$STUB_DIR/singularity"

# Stub condaload script (sourced via relative path from working dir)
mkdir -p "$WORK_DIR/workflow/scripts"
printf '#!/bin/bash\n' > "$WORK_DIR/workflow/scripts/condaload_CVRCseq.sh"

# Stub snakemake profile config so --profile doesn't error
mkdir -p "$WORK_DIR/config/profile"
printf 'jobs: 1\n' > "$WORK_DIR/config/profile/config.yaml"

# Stub config so workflow and singularity_image can be written via sed
mkdir -p "$WORK_DIR/config"
printf 'workflow: "RNAseq_PE"\nsingularity_image: ""\n' > "$WORK_DIR/config/config.yaml"

# Dummy sif image file for container mode tests
touch "$WORK_DIR/CVRCseq.sif"

# ---------------------------------------------------------------------------
# Test runner
# ---------------------------------------------------------------------------
run_test() {
    local description="$1"
    local expected_exit="$2"
    shift 2

    cd "$WORK_DIR"
    PATH="$STUB_DIR:$PATH" bash "$SCRIPT" "$@" > /dev/null 2>&1
    actual_exit=$?
    cd "$REPO_ROOT"

    if [[ $actual_exit -eq $expected_exit ]]; then
        echo "PASS: $description"
        ((PASS++))
    else
        echo "FAIL: $description (expected exit $expected_exit, got $actual_exit)"
        ((FAIL++))
    fi
}

run_test_default_container() {
    local description="$1"
    local expected_exit="$2"
    shift 2

    cd "$WORK_DIR"
    CVRCSEQ_SIF="$WORK_DIR/CVRCseq.sif" PATH="$STUB_DIR:$PATH" bash "$SCRIPT" "$@" > /dev/null 2>&1
    actual_exit=$?
    cd "$REPO_ROOT"

    if [[ $actual_exit -eq $expected_exit ]]; then
        echo "PASS: $description"
        ((PASS++))
    else
        echo "FAIL: $description (expected exit $expected_exit, got $actual_exit)"
        ((FAIL++))
    fi
}

echo
echo "Running snakemake_init.sh argument validation tests..."
echo "-------------------------------------------------------"

# Failure cases — script exits before source/snakemake calls
run_test "exits with code 1 when no arguments provided"            1
run_test "exits with code 1 when -d is missing"                    1  -w RNAseq_PE
run_test "exits with code 1 when -w is missing"                    1  -d /some/fastq/dir
run_test "exits with code 1 for unsupported workflow name"         1  -w InvalidWorkflow -d /some/fastq/dir

# Passing validation cases — use -c to skip cat_rename.py and reach snakemake stub
run_test "passes validation for RNAseq_PE"                         0  -w RNAseq_PE   -d /some/fastq/dir -c
run_test "passes validation for ATACseq_PE"                        0  -w ATACseq_PE  -d /some/fastq/dir -c
run_test "passes validation for CUT-RUN_PE"                        0  -w CUT-RUN_PE  -d /some/fastq/dir -c
run_test "passes validation for ChIPseq_PE"                        0  -w ChIPseq_PE  -d /some/fastq/dir -c
run_test "passes validation for RNAseqTE_PE"                       0  -w RNAseqTE_PE -d /some/fastq/dir -c
run_test "-c flag is boolean and does not consume -d value"         0  -w RNAseq_PE -c -d /some/fastq/dir
run_test_default_container "uses container by default when CVRCSEQ_SIF exists" 0 -w RNAseq_PE -d /some/fastq/dir -c
run_test "passes validation in container mode with valid -i image"  0  -w RNAseq_PE -d /some/fastq/dir -c -i "$WORK_DIR/CVRCseq.sif"
run_test "exits with code 1 for missing -i image path"              1  -w RNAseq_PE -d /some/fastq/dir -c -i "$WORK_DIR/missing.sif"

echo "-------------------------------------------------------"
echo "Results: $PASS passed, $FAIL failed"
echo

[[ $FAIL -eq 0 ]]
