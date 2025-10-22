import os
import tempfile
import pytest
import polars as pl
import pickle
from epymetheus import epymetheus
from epymetheus.epymetheus import MethylationOutput

def _normalize(s: str) -> str:
    return s.replace("\r\n", "\n").strip()

@pytest.fixture
def data_dir():
    here = os.path.dirname(__file__)
    return os.path.join(here, "..", "..", "epimetheus-cli", "tests", "data")


def test_methylation_pattern_median(data_dir, tmp_path):
    pileup = os.path.join(data_dir, "geobacillus-plasmids.pileup.bed")
    assembly = os.path.join(data_dir, "geobacillus-plasmids.assembly.fasta")
    expected = os.path.join(data_dir, "expected_out_median.tsv")

    outfile = tmp_path / "out.tsv"

    epymetheus.methylation_pattern(
        pileup,
        assembly,
        motifs = ["GATC_a_1", "GATC_m_3", "RGATCY_a_2"],
        output = str(outfile),
        threads = 1,
        output_type=MethylationOutput.Median
    )

    actual = outfile.read_text()
    expected_text = open(expected).read()
    assert _normalize(actual) == _normalize(expected_text)




def test_methylation_pattern_weighted_mean(data_dir, tmp_path):
    pileup = os.path.join(data_dir, "geobacillus-plasmids.pileup.bed")
    assembly = os.path.join(data_dir, "geobacillus-plasmids.assembly.fasta")
    expected = os.path.join(data_dir, "expected_out_weighted_mean.tsv")

    outfile = tmp_path / "out.tsv"

    epymetheus.methylation_pattern(
        pileup,
        assembly,
        output=str(outfile),
        threads = 1,
        motifs = ["GATC_a_1", "GATC_m_3", "RGATCY_a_2"],
        output_type=MethylationOutput.WeightedMean
    )

    actual = outfile.read_text()
    expected_text = open(expected).read()
    assert _normalize(actual) == _normalize(expected_text)


def test_methylation_pattern_weighted_mean_from_df(data_dir, tmp_path):
    pileup = os.path.join(data_dir, "geobacillus.bed.gz")
    assembly = os.path.join(data_dir, "geobacillus-plasmids.assembly.fasta")
    expected = os.path.join(data_dir, "expected_out_weighted_mean.tsv")

    outfile = tmp_path / "out.tsv"

    df = epymetheus.query_pileup_records(pileup, ["contig_2", "contig_3"])
    motifs = ["GATC_a_1", "GATC_m_3", "RGATCY_a_2"]

    result = epymetheus.methylation_pattern_from_dataframe(
        df,
        assembly,
        motifs = motifs,
        threads = 1,
        min_valid_read_coverage=3,
        min_valid_cov_to_diff_fraction=0.8,
        output_type=MethylationOutput.WeightedMean
    )

    result = result.sort(["contig", "motif", "mod_type"])
    result.write_csv(outfile, separator = "\t")

    actual = outfile.read_text()
    expected_text = open(expected).read()
    assert _normalize(actual) == _normalize(expected_text)


def test_methylation_output_pickle():
    """Test that MethylationOutput enum variants can be pickled and unpickled correctly."""
    # Test all enum variants
    variants = [
        MethylationOutput.Raw,
        MethylationOutput.Median,
        MethylationOutput.WeightedMean
    ]

    # Check module information
    print(f"MethylationOutput.__module__: {getattr(MethylationOutput, '__module__', 'NOT SET')}")
    print(f"MethylationOutput.__qualname__: {getattr(MethylationOutput, '__qualname__', 'NOT SET')}")

    for variant in variants:
        print(f"Testing variant: {variant}")
        print(f"Variant type: {type(variant)}")
        print(f"Variant module: {getattr(type(variant), '__module__', 'NOT SET')}")

        # Pickle and unpickle
        pickled = pickle.dumps(variant)
        unpickled = pickle.loads(pickled)

        # Verify the unpickled object is equal to the original
        assert unpickled == variant
        assert type(unpickled) == type(variant)

        # Verify string representation is preserved
        assert str(unpickled) == str(variant)
