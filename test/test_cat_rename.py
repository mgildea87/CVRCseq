"""Unit tests for workflow/scripts/cat_rename.py"""
import os
import sys
import gzip
import tempfile
import pandas as pd
import pytest
from unittest.mock import patch

# Add the scripts directory to the path so cat_rename can be imported
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'workflow', 'scripts'))
import cat_rename


def make_sample_table_SE(fastq_dir):
    """Minimal single-end sample table pointing at files in fastq_dir."""
    return pd.DataFrame({
        'Sample':       ['Sample1'],
        'Condition':    ['BL'],
        'Replicate':    ['1'],
        'File_Name_R1': ['sample1_R1.fastq.gz'],
    })


def make_sample_table_PE(fastq_dir):
    """Minimal paired-end sample table pointing at files in fastq_dir."""
    return pd.DataFrame({
        'Sample':       ['Sample1'],
        'Condition':    ['BL'],
        'Replicate':    ['1'],
        'File_Name_R1': ['sample1_R1.fastq.gz'],
        'File_Name_R2': ['sample1_R2.fastq.gz'],
    })


def make_sample_table_ChIP(fastq_dir):
    """Minimal ChIP/CUT-RUN sample table pointing at files in fastq_dir."""
    return pd.DataFrame({
        'Sample':       ['Sample1'],
        'Condition':    ['BL'],
        'Replicate':    ['1'],
        'Antibody':     ['H3K27ac'],
        'File_Name_R1': ['sample1_R1.fastq.gz'],
        'File_Name_R2': ['sample1_R2.fastq.gz'],
    })


def create_empty_fastq(path):
    """Create an empty gzipped file to simulate a FASTQ."""
    with gzip.open(path, 'wb') as f:
        f.write(b'')


# ---------------------------------------------------------------------------
# concat()
# ---------------------------------------------------------------------------

class TestConcat:
    def test_single_lane_passthrough(self):
        """Single-lane file (no L00 tag) is passed through unchanged."""
        with tempfile.TemporaryDirectory() as base:
            src = os.path.join(base, 'raw') + '/'
            os.makedirs(src)
            fname = 'sample1_R1.fastq.gz'
            create_empty_fastq(os.path.join(src, fname))
            out_dir = os.path.join(base, 'RNAseq_PE', 'inputs', 'fastq')
            os.makedirs(out_dir)
            with patch.object(sys, 'argv', ['cat_rename.py', src, os.path.join(base, 'RNAseq_PE')]):
                result = cat_rename.concat(make_sample_table_SE(src))
            assert len(result) == 1
            assert result[0].endswith(fname)
            assert os.path.exists(result[0])

    def test_multi_lane_merging(self):
        """Two L001/L002 lanes for the same file are merged into one output."""
        with tempfile.TemporaryDirectory() as base:
            src = os.path.join(base, 'raw') + '/'
            os.makedirs(src)
            create_empty_fastq(os.path.join(src, 'sample1_L001_R1.fastq.gz'))
            create_empty_fastq(os.path.join(src, 'sample1_L002_R1.fastq.gz'))
            out_dir = os.path.join(base, 'RNAseq_PE', 'inputs', 'fastq')
            os.makedirs(out_dir)
            with patch.object(sys, 'argv', ['cat_rename.py', src, os.path.join(base, 'RNAseq_PE')]):
                result = cat_rename.concat(make_sample_table_SE(src))
            # Both lanes should be merged into one output file
            assert len(result) == 1

    def test_non_fastq_files_ignored(self):
        """Non-fastq.gz files in the source directory are ignored."""
        with tempfile.TemporaryDirectory() as base:
            src = os.path.join(base, 'raw') + '/'
            os.makedirs(src)
            create_empty_fastq(os.path.join(src, 'sample1_R1.fastq.gz'))
            open(os.path.join(src, 'README.txt'), 'w').close()
            open(os.path.join(src, 'sample1.bam'), 'w').close()
            out_dir = os.path.join(base, 'RNAseq_PE', 'inputs', 'fastq')
            os.makedirs(out_dir)
            with patch.object(sys, 'argv', ['cat_rename.py', src, os.path.join(base, 'RNAseq_PE')]):
                result = cat_rename.concat(make_sample_table_SE(src))
            assert len(result) == 1


# ---------------------------------------------------------------------------
# rename_RNA_SE()
# ---------------------------------------------------------------------------

class TestRenameRNASE:
    def _argv(self, base):
        return ['cat_rename.py', base + '/', os.path.join(base, 'RNAseq_SE')]

    def test_renames_correctly(self):
        """File is copied to Sample_Condition_Replicate_R1.fastq.gz."""
        with tempfile.TemporaryDirectory() as base:
            fastq_dir = os.path.join(base, 'RNAseq_SE', 'inputs', 'fastq')
            os.makedirs(fastq_dir)
            create_empty_fastq(os.path.join(fastq_dir, 'sample1_R1.fastq.gz'))
            expected = os.path.join(fastq_dir, 'Sample1_BL_1_R1.fastq.gz')
            with patch.object(sys, 'argv', self._argv(base)):
                cat_rename.rename_RNA_SE(make_sample_table_SE(fastq_dir), [])
            assert os.path.exists(expected)

    def test_skips_if_output_exists(self):
        """Copy is skipped when the target file already exists."""
        with tempfile.TemporaryDirectory() as base:
            fastq_dir = os.path.join(base, 'RNAseq_SE', 'inputs', 'fastq')
            os.makedirs(fastq_dir)
            create_empty_fastq(os.path.join(fastq_dir, 'sample1_R1.fastq.gz'))
            expected = os.path.join(fastq_dir, 'Sample1_BL_1_R1.fastq.gz')
            create_empty_fastq(expected)
            mtime_before = os.path.getmtime(expected)
            with patch.object(sys, 'argv', self._argv(base)):
                cat_rename.rename_RNA_SE(make_sample_table_SE(fastq_dir), [])
            assert os.path.getmtime(expected) == mtime_before  # not overwritten

    def test_exits_on_missing_file(self):
        """sys.exit(1) is called when neither source nor target file exists."""
        with tempfile.TemporaryDirectory() as base:
            fastq_dir = os.path.join(base, 'RNAseq_SE', 'inputs', 'fastq')
            os.makedirs(fastq_dir)
            # No FASTQ file created
            with patch.object(sys, 'argv', self._argv(base)):
                with pytest.raises(SystemExit) as exc_info:
                    cat_rename.rename_RNA_SE(make_sample_table_SE(fastq_dir), [])
            assert exc_info.value.code == 1

    def test_removes_concat_files(self):
        """Temporary concat files are removed after renaming."""
        with tempfile.TemporaryDirectory() as base:
            fastq_dir = os.path.join(base, 'RNAseq_SE', 'inputs', 'fastq')
            os.makedirs(fastq_dir)
            create_empty_fastq(os.path.join(fastq_dir, 'sample1_R1.fastq.gz'))
            tmp_file = os.path.join(fastq_dir, 'tmp_concat.fastq.gz')
            create_empty_fastq(tmp_file)
            with patch.object(sys, 'argv', self._argv(base)):
                cat_rename.rename_RNA_SE(make_sample_table_SE(fastq_dir), [tmp_file])
            assert not os.path.exists(tmp_file)


# ---------------------------------------------------------------------------
# rename_RNA_PE()
# ---------------------------------------------------------------------------

class TestRenameRNAPE:
    def _argv(self, base):
        return ['cat_rename.py', base + '/', os.path.join(base, 'RNAseq_PE')]

    def test_renames_r1_and_r2(self):
        """Both R1 and R2 files are renamed correctly."""
        with tempfile.TemporaryDirectory() as base:
            fastq_dir = os.path.join(base, 'RNAseq_PE', 'inputs', 'fastq')
            os.makedirs(fastq_dir)
            create_empty_fastq(os.path.join(fastq_dir, 'sample1_R1.fastq.gz'))
            create_empty_fastq(os.path.join(fastq_dir, 'sample1_R2.fastq.gz'))
            with patch.object(sys, 'argv', self._argv(base)):
                cat_rename.rename_RNA_PE(make_sample_table_PE(fastq_dir), [])
            assert os.path.exists(os.path.join(fastq_dir, 'Sample1_BL_1_R1.fastq.gz'))
            assert os.path.exists(os.path.join(fastq_dir, 'Sample1_BL_1_R2.fastq.gz'))

    def test_exits_on_missing_r2(self):
        """sys.exit(1) when R2 file is missing."""
        with tempfile.TemporaryDirectory() as base:
            fastq_dir = os.path.join(base, 'RNAseq_PE', 'inputs', 'fastq')
            os.makedirs(fastq_dir)
            create_empty_fastq(os.path.join(fastq_dir, 'sample1_R1.fastq.gz'))
            # R2 intentionally not created
            with patch.object(sys, 'argv', self._argv(base)):
                with pytest.raises(SystemExit) as exc_info:
                    cat_rename.rename_RNA_PE(make_sample_table_PE(fastq_dir), [])
            assert exc_info.value.code == 1


# ---------------------------------------------------------------------------
# rename_ChIP()
# ---------------------------------------------------------------------------

class TestRenameChIP:
    def _argv(self, base):
        return ['cat_rename.py', base + '/', os.path.join(base, 'ChIPseq_PE')]

    def test_renames_with_antibody(self):
        """Antibody column is included in the output filename."""
        with tempfile.TemporaryDirectory() as base:
            fastq_dir = os.path.join(base, 'ChIPseq_PE', 'inputs', 'fastq')
            os.makedirs(fastq_dir)
            create_empty_fastq(os.path.join(fastq_dir, 'sample1_R1.fastq.gz'))
            create_empty_fastq(os.path.join(fastq_dir, 'sample1_R2.fastq.gz'))
            with patch.object(sys, 'argv', self._argv(base)):
                cat_rename.rename_ChIP(make_sample_table_ChIP(fastq_dir), [])
            assert os.path.exists(os.path.join(fastq_dir, 'Sample1_BL_1_H3K27ac_R1.fastq.gz'))
            assert os.path.exists(os.path.join(fastq_dir, 'Sample1_BL_1_H3K27ac_R2.fastq.gz'))


# ---------------------------------------------------------------------------
# main() argument validation
# ---------------------------------------------------------------------------

class TestMainArgValidation:
    def test_exits_with_no_args(self):
        """main() exits with usage message when no arguments are provided."""
        with patch.object(sys, 'argv', ['cat_rename.py']):
            with pytest.raises(SystemExit) as exc_info:
                cat_rename.main()
        assert exc_info.value.code == 1

    def test_exits_with_one_arg(self):
        """main() exits when only one argument is provided."""
        with patch.object(sys, 'argv', ['cat_rename.py', '/some/dir']):
            with pytest.raises(SystemExit) as exc_info:
                cat_rename.main()
        assert exc_info.value.code == 1
