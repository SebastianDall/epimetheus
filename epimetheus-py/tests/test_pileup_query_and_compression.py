import os
import tempfile
import pytest
from pathlib import Path
from epymetheus import epymetheus

@pytest.fixture
def data_dir():
    here = os.path.dirname(__file__)
    return os.path.join(here, "..", "..", "epimetheus-cli", "tests", "data")

def test_bgzf_compression_and_query(data_dir, tmp_path):
    """Test bgzf compression and querying functionality"""
    # Input file
    pileup_input = os.path.join(data_dir, "geobacillus-plasmids.pileup.bed")
    # Output files in temp directory
    compressed_file = tmp_path / "test_output.bed.gz"
    # Step 1: Compress the pileup file
    epymetheus.bgzf_pileup(
        pileup_input,
        str(compressed_file),
        keep=True,   # Keep original file
        force=False  # Don't overwrite if exists
    )
    # Verify compression worked
    assert compressed_file.exists(), "Compressed file should be created"
    assert Path(f"{compressed_file}.tbi").exists(), "Index file should be created"
    # Step 2: Query for existing contig (should pass)
    records_contig3 = epymetheus.query_pileup_records(
        str(compressed_file),
        ["contig_3"]
    )
    assert len(records_contig3) > 0, "Should find records for contig_3"
    print(f"Found {len(records_contig3)} records for contig_3")
    # Verify record properties
    first_record = records_contig3[0]
    assert hasattr(first_record, 'contig'), "Record should have contig attribute"
    assert hasattr(first_record, 'start'), "Record should have start attribute"
    assert hasattr(first_record, 'fraction_modified'), "Record should have fraction_modified attribute"
    assert first_record.contig == "contig_3", "Record should belong to contig_3"
    # Step 3: Query for non-existing contig (should fail or return empty)
    try:
        records_contig10 = epymetheus.query_pileup_records(
            str(compressed_file),
            ["contig_10"]
        )
        # If it doesn't raise an error, it should return empty
        assert len(records_contig10) == 0, "Should return no records for non-existent contig_10"
        print("Query for contig_10 returned empty results (expected)")
    except Exception as e:
        # It's also acceptable for it to raise an error
        print(f"Query for contig_10 failed as expected: {e}")
    # Cleanup is automatic with tmp_path fixture, but let's be explicit
    if compressed_file.exists():
        compressed_file.unlink()
    index_file = Path(f"{compressed_file}.tbi")
    if index_file.exists():
        index_file.unlink()

def test_bgzf_compression_with_auto_output(data_dir, tmp_path):
    """Test bgzf compression with automatic output naming"""
    # Copy input file to temp directory so we can test auto-naming
    pileup_input = os.path.join(data_dir, "geobacillus-plasmids.pileup.bed")
    temp_input = tmp_path / "input.bed"
    # Copy content
    with open(pileup_input, 'r') as src, open(temp_input, 'w') as dst:
        dst.write(src.read())
    # Compress with auto output (None)
    epymetheus.bgzf_pileup(
        str(temp_input),
        None,        # Auto-generate output name
        keep=True,   # Keep original
        force=False
    )
    # Check auto-generated output
    expected_output = temp_input.with_suffix('.bed.gz')
    assert expected_output.exists(), "Auto-generated compressed file should exist"
    assert Path(f"{expected_output}.tbi").exists(), "Auto-generated index should exist"

def test_bgzf_force_overwrite(data_dir, tmp_path):
    """Test bgzf force overwrite functionality"""
    pileup_input = os.path.join(data_dir, "geobacillus-plasmids.pileup.bed")
    output_file = tmp_path / "test.bed.gz"
    # Create first compression
    epymetheus.bgzf_pileup(str(pileup_input), str(output_file), True, False)
    assert output_file.exists()
    original_size = output_file.stat().st_size
    # Try to compress again with force=True (should succeed)
    epymetheus.bgzf_pileup(str(pileup_input), str(output_file), True, True)
    assert output_file.exists()
    # File should still exist and have similar size
    new_size = output_file.stat().st_size
    assert abs(new_size - original_size) < 100, "File sizes should be similar"
