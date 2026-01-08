# epimetheus_methylome

A Rust library for DNA motif representation and methylation analysis, with support for IUPAC nucleotide codes and Nanopore methylation data.

## Features

- **DNA Motif Handling**: Create, manipulate, and search for DNA motifs with modification sites
- **IUPAC Support**: Full support for IUPAC nucleotide codes (including ambiguous bases like R, Y, N, etc.)
- **Methylation Types**: Support for common DNA modifications:
  - 6mA (N6-methyladenine)
  - 5mC (5-methylcytosine)
  - 4mC (4-methylcytosine)
- **Motif Operations**: Reverse complement, regex conversion, motif clustering/hierarchy
- **Nanopore Integration**: Parse methylation data from FASTQ records with MM/ML tags

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
epimetheus_methylome = "1.1.0"
```

## Usage

### Creating and Working with Motifs

```rust
use epimetheus_methylome::{Motif, ModType};

// Create a GATC motif with 6mA modification at position 1
let motif = Motif::new("GATC", "a", 1).unwrap();
assert_eq!(motif.mod_type, ModType::SixMA);

// Get reverse complement
let rev_comp = motif.reverse_complement();
assert_eq!(rev_comp.sequence_to_string(), "GATC");

// Convert to regex pattern
let regex = motif.to_regex();
assert_eq!(regex, "GATC");
```

### IUPAC Nucleotide Codes

```rust
use epimetheus_methylome::{IupacBase, Motif};

// Create motif with ambiguous bases
let motif = Motif::new("RGATCY", "m", 4).unwrap();
// R = A or G, Y = C or T

// Convert to regex
let regex = motif.to_regex();
assert_eq!(regex, "[AG]GATC[CT]");

// Get all possible DNA sequences
let sequences = motif.possible_dna_sequences();
// Returns all combinations: AGATCC, AGATCT, GGATCC, GGATCT
```

### Methylation Types

```rust
use epimetheus_methylome::ModType;

// Parse from pileup codes
let mod_type = "a".parse::<ModType>().unwrap();
assert_eq!(mod_type, ModType::SixMA);

// Get pileup code
assert_eq!(ModType::FiveMC.to_pileup_code(), "m");
assert_eq!(ModType::FourMC.to_pileup_code(), "21839");
```

### Finding Motifs in Sequences

```rust
use epimetheus_methylome::{Motif, find_motif_indices_in_sequence};
use epimetheus_methylome::sequence::Sequence;

let sequence = Sequence::from_str("GGATCTCCATGATC").unwrap();
let motif = Motif::new("GATC", "m", 3).unwrap();

let indices = find_motif_indices_in_sequence(&sequence, &motif);
assert_eq!(indices, vec![4, 13]); // Returns modification positions
```

### Motif Hierarchy

```rust
use epimetheus_methylome::Motif;

// Check parent-child relationships
let parent = Motif::new("GATC", "a", 1).unwrap();
let child = Motif::new("RGATCY", "a", 2).unwrap();

assert!(parent.is_child_motif(&child));
```

### Working with Nanopore Data

```rust
use epimetheus_methylome::read::Read;
use noodles_fastq as fastq;

// Parse FASTQ record with MM/ML tags
let record = fastq::Record::new(
    fastq::record::Definition::new("read_id", "MM:Z:A+a.,0,1; ML:B:C,204,255"),
    "GGGCGGATCAGATC",
    "JIGSSGIH=GLLML"
);

let read = Read::from_fastq_record(record).unwrap();
let modifications = read.get_modifications();
// Access methylation calls at specific positions
// first and third A should have the methylation scores 204 and 255, respectively
```

## API Documentation

For detailed API documentation, see [docs.rs/epimetheus_methylome](https://docs.rs/epimetheus_methylome).

## Use Cases

This library is designed for:
- Bioinformatics pipelines processing Nanopore methylation data
- Restriction enzyme motif analysis
- Methylation pattern detection and quantification
- Motif clustering and hierarchical analysis

## License

MIT

## Contributing

This crate is part of the [Epimetheus](https://github.com/SebastianDall/epimetheus) project for analyzing methylation patterns from pileup files.
