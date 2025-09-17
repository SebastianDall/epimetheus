import os
import tempfile
import pytest
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
        str(outfile),
        1,
        motifs = ["GATC_a_1", "GATC_m_3", "RGATCY_a_2"],
        batch_size=1000,
        min_valid_read_coverage=3,
        min_valid_cov_to_diff_fraction=0.8,
        allow_assembly_pileup_mismatch=False,
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
        str(outfile),
        1,
        motifs = ["GATC_a_1", "GATC_m_3", "RGATCY_a_2"],
        batch_size=1000,
        min_valid_read_coverage=3,
        min_valid_cov_to_diff_fraction=0.8,
        allow_assembly_pileup_mismatch=False,
        output_type=MethylationOutput.WeightedMean
    )

    actual = outfile.read_text()
    expected_text = open(expected).read()
    assert _normalize(actual) == _normalize(expected_text)
