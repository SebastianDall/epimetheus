import os
import tempfile
import pytest
import polars as pl
from pathlib import Path
from epymetheus.epymetheus import *

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
    bgzf_pileup(
        pileup_input,
        str(compressed_file),
        keep=True,   # Keep original file
        force=False  # Don't overwrite if exists
    )
    # Verify compression worked
    assert compressed_file.exists(), "Compressed file should be created"
    assert Path(f"{compressed_file}.tbi").exists(), "Index file should be created"
    # Step 2: Query for existing contig (should pass)
    records_contig3 = query_pileup_records(
        str(compressed_file),
        ["contig_3"]
    )
    assert len(records_contig3) > 0, "Should find records for contig_3"
    print(f"Found {len(records_contig3)} records for contig_3")
    
    # Cleanup is automatic with tmp_path fixture, but let's be explicit
    if compressed_file.exists():
        compressed_file.unlink()
    index_file = Path(f"{compressed_file}.tbi")
    if index_file.exists():
        index_file.unlink()

def test_query_data(data_dir, tmp_path):
    """Test querying multiple contigs and the data"""
    pileup_input = os.path.join(data_dir, "geobacillus.bed.gz")

    df = query_pileup_records(pileup_input, contigs=["contig_2","contig_3"])

    assert df.columns == [
        "contig",
        "start",
        "end",
        "mod_type",
        "score",
        "strand",
        "start_pos",
        "end_pos",
        "color",
        "n_valid_cov",
        "fraction_modified",
        "n_modified",
        "n_canonical",
        "n_other_mod",
        "n_delete",
        "n_fail",
        "n_diff",
        "n_no_call",
    ], "Colums do not match"

    assert len(df.filter(pl.col("contig") == "contig_3")) > 0, "No records matched contig_3"
    assert len(df.filter(pl.col("contig") == "contig_2")) > 0, "No records matched contig_2"
    assert len(df.filter(pl.col("contig") == "contig_10")) == 0, "Records matched contig_10 contig 10 should not be present"

def test_query_with_col_filter(data_dir, tmp_path):
    """Test querying multiple contigs and the data"""
    pileup_input = os.path.join(data_dir, "geobacillus.bed.gz")

    df = query_pileup_records(pileup_input, contigs=["contig_2","contig_3"], columns = [PileupColumn.Contig, PileupColumn.Start, PileupColumn.Score])

    assert df.columns == [
        "contig",
        "start",
        "score",
    ], "Colums do not match"

    assert len(df.filter(pl.col("contig") == "contig_3")) > 0, "No records matched contig_3"
    assert len(df.filter(pl.col("contig") == "contig_2")) > 0, "No records matched contig_2"
    assert len(df.filter(pl.col("contig") == "contig_10")) == 0, "Records matched contig_10 contig 10 should not be present"
    



def test_bgzf_compression_with_auto_output(data_dir, tmp_path):
    """Test bgzf compression with automatic output naming"""
    # Copy input file to temp directory so we can test auto-naming
    pileup_input = os.path.join(data_dir, "geobacillus-plasmids.pileup.bed")
    temp_input = tmp_path / "input.bed"
    # Copy content
    with open(pileup_input, 'r') as src, open(temp_input, 'w') as dst:
        dst.write(src.read())
    # Compress with auto output (None)
    bgzf_pileup(
        str(temp_input),
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
    bgzf_pileup(str(pileup_input), str(output_file), True, False)
    assert output_file.exists()
    original_size = output_file.stat().st_size
    # Try to compress again with force=True (should succeed)
    bgzf_pileup(str(pileup_input), str(output_file), True, True)
    assert output_file.exists()
    # File should still exist and have similar size
    new_size = output_file.stat().st_size

def test_bgzf_compression_from_lines(data_dir, tmp_path):
    """Test bgzf streaming compression with contig name transformation"""
    pileup_input = os.path.join(data_dir, "geobacillus.bed.gz")
    temp_output = tmp_path / "output_lines.bed.gz"

    # Verify input file exists
    assert os.path.exists(pileup_input), f"Input file not found: {pileup_input}"

    writer = BgzfWriter(str(temp_output), force=True)

    original_record_counts = {}
    processed_record_counts = {}

    try:
        for c in ["contig_2", "contig_3"]:
            records = query_pileup_records(pileup_input, [c], None)
            assert len(records) > 0, f"No records found for contig {c} in input file"

            required_columns = ["contig","start","end","mod_type","score","strand","start_pos","end_pos","color","n_valid_cov","fraction_modified","n_modified","n_canonical","n_other_mod","n_delete","n_fail","n_diff","n_no_call"]
            assert records.columns == required_columns, "Columns does not match"
            
            original_record_counts[c] = len(records)

            mut_records = records.with_columns(
                pl.when(pl.col("contig") == "contig_2")
                    .then(pl.lit("contig_4"))
                    .otherwise(pl.lit("contig_5")).alias("contig")
            )

            csv_string = mut_records.write_csv(separator="\t", include_header=False)
            lines = csv_string.strip().split('\n')

            writer.write_lines(lines)
            processed_record_counts[c] = len(lines)

        writer.finish()

        # Verify output files exist
        assert temp_output.exists(), "Compressed output file should exist"
        assert Path(f"{temp_output}.tbi").exists(), "Tabix index file should exist"

        # Verify transformed data
        expected_contigs = ["contig_4", "contig_5"]
        contigs_in_output = query_pileup_records(str(temp_output), expected_contigs, None)
        
        assert contigs_in_output.shape[0] > 0, "No records found in compressed output"
        
        unique_contigs = contigs_in_output["contig"].unique().to_list()
        assert sorted(unique_contigs) == sorted(expected_contigs), f"Expected contigs {expected_contigs}, got {unique_contigs}"
        
        # Verify record counts are preserved
        output_counts = {}
        for contig in expected_contigs:
            contig_records = contigs_in_output.filter(pl.col("contig") == contig)
            output_counts[contig] = contig_records.shape[0]
        
        # Map back to original contigs for comparison
        original_contig_2_count = original_record_counts["contig_2"]
        original_contig_3_count = original_record_counts["contig_3"]
        
        assert output_counts["contig_4"] == original_contig_2_count, f"Record count mismatch for contig_4: expected {original_contig_2_count}, got {output_counts['contig_4']}"
        assert output_counts["contig_5"] == original_contig_3_count, f"Record count mismatch for contig_5: expected {original_contig_3_count}, got {output_counts['contig_5']}"
        
        print(f"Successfully processed {sum(original_record_counts.values())} records")
        print(f"Original counts: {original_record_counts}")
        print(f"Output counts: {output_counts}")

    except Exception as e:
        # Cleanup on failure
        if temp_output.exists():
            temp_output.unlink()
        tbi_file = Path(f"{temp_output}.tbi")
        if tbi_file.exists():
            tbi_file.unlink()
        raise e
