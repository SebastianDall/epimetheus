import os
import tempfile
import pytest
from epymetheus import epymetheus

def _normalize(s: str) -> str:
    return s.strip().replace("\r\n", "\n")

def test_methylation_pattern(tmp_path):
    outfile = tmp_path / "out.tsv"

    epymetheus.remove_child_motifs(
        output = str(outfile),
        motifs = ["GATC_a_1", "GATC_m_3", "RGATCY_a_2"],
    )

    actual = outfile.read_text()
    expected = "motif\tmod_type\tmod_position\nGATC\ta\t1\nGATC\tm\t3"    
    assert _normalize(actual) == _normalize(expected)   

