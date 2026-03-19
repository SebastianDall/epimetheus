[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/epimetheus-py/README.html)
![PyPI - Version](https://img.shields.io/pypi/v/epimetheus-py)

# epymetheus

Python bindings for [epimetheus](https://github.com/SebastianDall/epimetheus) — a tool for converting Nanopore pileup files (from [modkit](https://github.com/nanoporetech/modkit)) to per-motif methylation data.

## Installation

```bash
pip install epimetheus-py

# or via conda

conda install -c conda-forge -c bioconda epimetheus-py
```

### Build from source

> Requires the Rust/Cargo toolchain.

```bash
git clone https://github.com/SebastianDall/epimetheus.git
cd epimetheus

python -m venv epimetheus-venv
source epimetheus-venv/bin/activate

pip install maturin
maturin develop -m epimetheus-py/Cargo.toml
```

## API

### `methylation_pattern`

Extract per-motif methylation statistics from a pileup file and assembly.

```python
import epymetheus
from epymetheus.epymetheus import MethylationOutput

df = epymetheus.methylation_pattern(
    pileup="pileup.bed.gz",
    assembly="assembly.fasta",
    motifs=["GATC_a_1", "CCWGG_m_1"],
    output_type=MethylationOutput.Median,
    threads=4,
)
```

The `assembly` argument accepts either a file path or a `dict[str, str | SeqRecord]` (e.g. loaded with `Bio.SeqIO`).

**Output columns (Median / WeightedMean):**

| Column | Description |
|--------|-------------|
| `contig` | Contig ID |
| `motif` | Motif sequence |
| `mod_type` | Modification type (pileup code) |
| `mod_position` | Modified base position in motif |
| `methylation_value` | Median or weighted-mean methylation fraction |
| `mean_read_cov` | Mean read coverage at used positions |
| `n_motif_obs` | Motif positions above `min_valid_read_coverage` |
| `motif_occurences_total` | Total motif occurrences in contig |

**Output columns (Raw):**

| Column | Description |
|--------|-------------|
| `contig` | Contig ID |
| `start` | Position (0-based) |
| `strand` | Strand (`+` / `-`) |
| `motif` | Motif sequence |
| `mod_type` | Modification type |
| `mod_position` | Modified base position in motif |
| `n_modified` | Number of modified reads |
| `n_valid_cov` | Number of valid coverage reads |
| `n_diff` | Number of reads with a base other than the non-canonical |
| `n_fail` | Number of calls where the base was below the threshold |


**Output types:**

| Value | Description |
|-------|-------------|
| `MethylationOutput.Median` | Median of per-position methylation fractions |
| `MethylationOutput.WeightedMean` | Coverage-weighted mean methylation |
| `MethylationOutput.Raw` | All positions with raw counts |

**Key parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `threads` | `1` | Parallel threads |
| `min_valid_read_coverage` | `3` | Minimum reads required at a position |
| `min_valid_cov_to_diff_fraction` | `0.8` | Minimum `n_valid_cov / (n_valid_cov + n_diff)` |
| `contigs` | `None` | Optional list of contig IDs to restrict processing |
| `output` | `None` | Optional path to write output TSV |
| `allow_assembly_pileup_mismatch` | `False` | Continue if a pileup contig is absent from the assembly |

---

### `methylation_pattern_from_dataframe`

Same as `methylation_pattern` but accepts a Polars DataFrame (e.g. from `query_pileup_records`) instead of a file path.

```python
df = epymetheus.methylation_pattern_from_dataframe(
    pileup_df=pileup_df,
    assembly="assembly.fasta",
    motifs=["GATC_a_1"],
    output_type=MethylationOutput.Median,
    threads=4,
)
```

---

### `query_pileup_records`

Query specific contigs from a BGZF-compressed pileup file and return a Polars DataFrame.

```python
df = epymetheus.query_pileup_records(
    pileup_path="pileup.bed.gz",
    contigs=["contig_1", "contig_2"],
)
```

Use the `columns` parameter to select a subset of columns and reduce memory usage:

```python
from epymetheus.epymetheus import PileupColumn

df = epymetheus.query_pileup_records(
    pileup_path="pileup.bed.gz",
    contigs=["contig_1"],
    columns=[PileupColumn.NValidCov, PileupColumn.NModified],
)
```

---

### `bgzf_pileup`

Compress a pileup BED file to BGZF format for fast random access. Strongly recommended when the pileup is queried multiple times (speeds up `methylation_pattern` by ~6x).

```python
epymetheus.bgzf_pileup(
    input="pileup.bed",
    output="pileup.bed.gz",
    keep=True,
    force=False,
)
```

---

### `BgzfWriter`

Write pileup lines incrementally to a BGZF file.

```python
writer = epymetheus.BgzfWriter("pileup.bed.gz", force=True)
writer.write_lines(["line1", "line2"])
writer.finish()  # writes the tabix index (.tbi) and finalises the file
```

---

### `remove_child_motifs`

Remove redundant child motifs from a results file. A child motif is one whose sequence is contained within a parent motif (e.g. `GATC_a_1` is a child of `RGATCY_a_2`).

```python
epymetheus.remove_child_motifs(
    output="methylation_pattern.tsv",
    motifs=["GATC_a_1", "RGATCY_a_2"],
)
```

---

### Motif format

Motifs are specified as `<sequence>_<mod_type>_<mod_position>`:

- `GATC_a_1` — adenine methylation at position 1 of GATC
- `CCWGG_m_1` — cytosine methylation at position 1 of CCWGG (W = A or T)

IUPAC ambiguity codes are supported in motif sequences.

---

### Version

```python
import epymetheus
print(epymetheus.__version__)
```
