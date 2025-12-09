[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/epimetheus-py/README.html)
![PyPI - Version](https://img.shields.io/pypi/v/epimetheus-py)
# epimetheus
`epimetheus` is a CLI tool for converting a pileup file ([modkit](https://github.com/nanoporetech/modkit)), where reads with modification calls are mapped to an assembly, to motif methylation instead.

To use `epimetheus` it is therefore necessary to know which motifs to look for beforehand (epimetheus meaning hindsight or afterthought). [Nanomotif](https://github.com/MicrobialDarkMatter/nanomotif) can find motifs in metagenomic data using Nanopore methylation calls.
## Installation
`epimetheus` is distributed as a python package and CLI tool.

> If using nixos the entire project should be contained in the `flake.nix` + `.envrc`.
### Python package
This project is distributed as a python package using `pyo3`, which can be build using `maturin develop -m epimetheus-py/Cargo.toml`.
Can be installed via `conda` or `pypi` using:

```bash
pip install epimetheus-py

# or

conda install -c conda-forge -c bioconda epimetheus-py
```

#### Develop from source
To create the python wheel from source follow the steps below:

> Note this requires the cargo tool chain


```bash
git clone https://github.com/SebastianDall/epimetheus.git

# Create a python venv or install maturin
python -m venv epimetheus-venv

pip install maturing develop

cd epimetheus-py
maturin develop # This will build and install the package.
```


### CLI tool
The CLI tool comes pre-compiled for musl-linux. Just download the cli tool in release.

To build from source run:

```bash
git clone https://github.com/SebastianDall/epimetheus.git

cargo install --locked --path epimetheus-cli
```



## CLI Usage:
```bash
Usage: epimetheus <COMMAND>

Commands:
  methylation-pattern  
  motif-cluster        
  bgzip                
  help                 Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

### BGZF (recommended)
The bed files can be compressed (often a factor of 8-10) to a tabbed gz file (BGZF) which allows for fast lookup with the hts-lib.
To compress the file use:
```bash
Usage: epimetheus bgzip compress [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>    Path to output pileup file. [.bed].
  -o, --output <OUTPUT>  Path to output pileup file [.bed.gz].
      --keep             Setting flag will keep the original uncompressed file.
      --force            Setting flag will override the file if exists.
  -h, --help             Print help
```


This will allow for fast lookup, which speeds up the `methylation-pattern` by a factor of 6.
This is highly recommended if pileup is accessed multiple times which it will in `nanomotif`.

To look up a specific contig use the `decompress` command:

```bash
Usage: epimetheus bgzip decompress [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>                Path to output pileup file. [.bed.gz].
  -o, --output <OUTPUT>              Path to output pileup file [.bed].
      --ls                           list contig names in pileup.
      --contigs <CONTIGS>...         Optional vector of contig ids to query. Left empty the whole pileup will be read.
      --contigs-file <CONTIGS_FILE>  File with contig names in it.
  -h, --help                         Print help
```

### methylation pattern
The motif methylation can be searched for on read and contig level.

```bash
Usage: epimetheus methylation-pattern <COMMAND>

Commands:
  contig  
  read    
  help    Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help

```
#### Contig level

Efficient processing of a pileup file for finding the read methylation degree of a motif for all contigs. Supply the assembly, the pileup and the motifs of interest. The tool will:
 - Find motif occurences
 - Find the number of reads and mean read methylation at each position
 - calculate the median of mean methylated positions.

The return is a dataframe with:
- contig: The contig id
- motif: The motif sequence
- mod_type: The modification type [6mA, 5mC, 4mC (as pileup codes)]
- mod_position: The modification position in the motif sequence
- methylation_value: see below
- mean_read_cov: The mean read coverage for positions used in the median calculation
- n_motif_obs: The number of motifs with methylation information above `min-valid-read-coverage`
- motif_occcurences_total: The total number of occurences of motif in contig.


Three output types are available:
- median: Firstly the fraction of reads at motif positions is calculated and the median of these are returned.
- weighted-mean: the fraction of reads modified weighted by the n_valid_coverage at those positions.
- raw: Outputs the all motif positions and their n_modified and n_valid_cov
```bash
Usage: epimetheus methylation-pattern contig [OPTIONS] --pileup <PILEUP> --assembly <ASSEMBLY> --output <OUTPUT> --motifs <MOTIFS>...

Options:
  -p, --pileup <PILEUP>
          Path to pileup. Can be .bed.gz (recommended see bgzip command) or .bed
  -a, --assembly <ASSEMBLY>
          Path to assembly.
  -o, --output <OUTPUT>
          Path to output file. Must be .tsv.
  -t, --threads <THREADS>
          Number of parallel tasks. [default: 1]
  -m, --motifs <MOTIFS>...
          Supply chain of motifs as <motif>_<mod_type>_<mod_position>. Example: '-m GATC_a_1 RGATCY_a_2'
      --min-valid-read-coverage <MIN_VALID_READ_COVERAGE>
          Minimum valid read coverage for calculating methylation. [default: 3]
      --batch-size <BATCH_SIZE>
          Number of contigs to process at a time. Higher number will use more RAM. [default: 1000]
      --min-valid-cov-to-diff-fraction <MIN_VALID_COV_TO_DIFF_FRACTION>
          Required fraction of valid coverage relative to different read mapping. N_valid_cov / (N_valid_cov + N_diff) [default: 0.8]
      --allow-mismatch
          Allow epimetheus to continue if a contig in the pileup is not present in the assembly
      --output-type <OUTPUT_TYPE>
          Specify the type of methylation output type. Raw will give all motif methylations for each contig. [default: median] [possible values: raw, median, weighted-mean]
  -h, --help
          Print help
```


#### Read level
This mode first searches for motif occurences in reads and then returns the quality of the methylation call from the basecaller at that position [0-255]

The input required is an indexed bam file.

The output is:
- contig_id: The contig id where the read is mapped.
- read_id
- read_length
- motif
- mod_type
- mod_position
- quality: u8, 255: 100% confidence [0-255]

> If the contig has no mapped reads, currently no warning is produced.

```bash
Usage: epimetheus methylation-pattern read [OPTIONS] --input <INPUT> --output <OUTPUT> --motifs <MOTIFS>...

Options:
  -i, --input <INPUT>            Path to bam file.
      --contig-ids <CONTIG_IDS>  File with specific contig ids to process.
  -o, --output <OUTPUT>          Path to output file. Must be .tsv.
  -t, --threads <THREADS>        Number of parallel tasks. [default: 1]
  -m, --motifs <MOTIFS>...       Supply chain of motifs as <motif>_<mod_type>_<mod_position>. Example: '-m GATC_a_1 RGATCY_a_2'
  -h, --help                     Print help

```


### motif-cluster
Motif-cluster will collapse a list of provided motifs to a set of "parent" motifs.
A `parent` motif is contained within another motif which would be the `child` motif. For instance
provided a list of `GATC_m_3` and `RGATCY_m_4`, the resulting list of motifs will be `GATC_m_3` because
`RGATCY_m_4` is a `child` of  `GATC_m_3`.

For binning this can be beneficial to reduce correlated features.

```bash
Usage: epimetheus motif-cluster --output <OUTPUT> --motifs <MOTIFS>...

Options:
  -o, --output <OUTPUT>     Path to output file. Must be .tsv.
  -m, --motifs <MOTIFS>...  Supply chain of motifs as <motif>_<mod_type>_<mod_position>. Example: '-m GATC_a_1 RGATCY_a_2'
  -h, --help                Print help
```


