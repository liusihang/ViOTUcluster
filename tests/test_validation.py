#!/usr/bin/env python3
"""
Unit tests for ViOTUcluster validation and filtering modules.

Run with: python -m pytest tests/test_validation.py -v
"""

import pytest
import tempfile
import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ViOTUcluster.validation import (
    validate_file_exists,
    validate_fasta,
    validate_directory,
    validate_positive_integer,
)


class TestValidateFileExists:
    """Tests for validate_file_exists function."""
    
    def test_existing_file(self, tmp_path):
        """Test with an existing file."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")
        assert validate_file_exists(str(test_file)) is True
    
    def test_nonexistent_file(self):
        """Test with a non-existent file."""
        with pytest.raises(FileNotFoundError):
            validate_file_exists("/nonexistent/path/file.txt")
    
    def test_directory_instead_of_file(self, tmp_path):
        """Test with a directory path instead of file."""
        with pytest.raises(ValueError):
            validate_file_exists(str(tmp_path))


class TestValidateFasta:
    """Tests for validate_fasta function."""
    
    def test_valid_fasta(self, tmp_path):
        """Test with a valid FASTA file."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nACGTACGT\n>seq2\nTGCATGCA\n")
        assert validate_fasta(str(fasta_file)) is True
    
    def test_valid_fa_extension(self, tmp_path):
        """Test with .fa extension."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">seq1\nACGT\n")
        assert validate_fasta(str(fasta_file)) is True
    
    def test_empty_fasta(self, tmp_path):
        """Test with an empty FASTA file."""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text("")
        with pytest.raises(ValueError, match="empty"):
            validate_fasta(str(fasta_file))
    
    def test_invalid_extension(self, tmp_path):
        """Test with invalid file extension."""
        txt_file = tmp_path / "test.txt"
        txt_file.write_text(">seq1\nACGT\n")
        with pytest.raises(ValueError, match="invalid extension"):
            validate_fasta(str(txt_file))


class TestValidateDirectory:
    """Tests for validate_directory function."""
    
    def test_existing_directory(self, tmp_path):
        """Test with an existing directory."""
        assert validate_directory(str(tmp_path)) is True
    
    def test_nonexistent_directory_no_create(self):
        """Test with non-existent directory, create=False."""
        with pytest.raises(FileNotFoundError):
            validate_directory("/nonexistent/directory")
    
    def test_create_directory(self, tmp_path):
        """Test creating a new directory."""
        new_dir = tmp_path / "new_subdir"
        assert validate_directory(str(new_dir), create=True) is True
        assert os.path.isdir(new_dir)


class TestValidatePositiveInteger:
    """Tests for validate_positive_integer function."""
    
    def test_valid_integer(self):
        """Test with valid integer string."""
        assert validate_positive_integer("42") == 42
    
    def test_zero_with_default_min(self):
        """Test that 0 fails with default min_val=1."""
        with pytest.raises(ValueError):
            validate_positive_integer("0")
    
    def test_zero_with_min_zero(self):
        """Test that 0 passes with min_val=0."""
        assert validate_positive_integer("0", min_val=0) == 0
    
    def test_negative_integer(self):
        """Test with negative integer."""
        with pytest.raises(ValueError):
            validate_positive_integer("-5")
    
    def test_non_integer_string(self):
        """Test with non-integer string."""
        with pytest.raises(ValueError):
            validate_positive_integer("abc")
    
    def test_max_value(self):
        """Test max value constraint."""
        assert validate_positive_integer("5", max_val=10) == 5
        with pytest.raises(ValueError):
            validate_positive_integer("15", max_val=10)


class TestFilterContigs:
    """Tests for filter_contigs functionality."""
    
    def test_filter_by_length(self, tmp_path):
        """Test filtering FASTA sequences by length."""
        from ViOTUcluster.filter_contigs import _filter_single_fasta
        
        # Create test FASTA with short and long sequences
        input_file = tmp_path / "input.fasta"
        input_file.write_text(
            ">short_seq\n"
            "ACGT\n"
            ">long_seq\n"
            + "ACGT" * 100 + "\n"  # 400 bp
        )
        
        output_file = tmp_path / "output.fasta"
        _filter_single_fasta(100, str(input_file), str(output_file))
        
        # Read output and verify only long sequence remains
        from Bio import SeqIO
        records = list(SeqIO.parse(str(output_file), "fasta"))
        assert len(records) == 1
        assert records[0].id == "long_seq"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
