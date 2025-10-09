import os
import tempfile
import pytest
import polars as pl
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
