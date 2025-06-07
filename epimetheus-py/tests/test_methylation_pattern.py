import os
import tempfile
import pytest
from epymetheus import epymetheus

def _normalize(s: str) -> str:
    return s.replace("\r\n", "\n").strip()

@pytest.fixture
def data_dir():
    here = os.path.dirname(__file__)
    return os.path.join(here, "..", "..", "epimetheus-cli", "tests", "data")


def test_methylation_pattern(data_dir, tmp_path):
    pileup = os.path.join(data_dir, "geobacillus-plasmids.pileup.bed")
    assembly = os.path.join(data_dir, "geobacillus-plasmids.assembly.fasta")
    expected = os.path.join(data_dir, "expected_out.tsv")

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
        allow_assembly_pilup_mismatch=False,
    )

    actual = outfile.read_text()
    expected_text = open(expected).read()
    assert _normalize(actual) == _normalize(expected_text)   
