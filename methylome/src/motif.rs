use crate::{IupacBase, ModType, sequence::Sequence};
use anyhow::{Result, bail};
use std::str::FromStr;

pub type Position = u8;

/// Represents a biological motif, which includes a nucleotide sequence,
/// its modification type, and the position of the modification.
///
/// # Fields
/// - `sequence`: A vector of IUPAC bases representing the motif sequence.
/// - `mod_type`: The type of modification (e.g., 6mA, 5mC).
/// - `mod_position`: The position of the modification within the sequence (0-indexed).
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Motif {
    pub sequence: Sequence,
    pub mod_type: ModType,
    pub mod_position: Position,
}

impl Motif {
    /// Constructs a new `Motif` from a string sequence, modification type, and modification position.
    ///
    /// # Arguments
    /// - `sequence`: A string representing the nucleotide sequence (using IUPAC codes).
    /// - `mod_type`: A string representing the modification type (e.g., "a" (6mA), "m" (5mC), "21839" (4mC)0).
    /// - `mod_position`: The 0-indexed position of the modification in the sequence.
    ///
    /// # Errors
    /// Returns an error if:
    /// - The `sequence` contains invalid IUPAC codes.
    /// - The `mod_position` is out of bounds for the sequence.
    /// - The `mod_type` does not match the base at `mod_position` (e.g., 6mA must modify an 'A').
    ///
    /// # Examples
    /// ```
    /// use methylome::{Motif, ModType};
    ///
    /// let motif = Motif::new("GATC", "a", 1).unwrap();
    /// assert_eq!(motif.mod_type, ModType::SixMA);
    /// ```
    pub fn new(sequence_str: &str, mod_type: &str, mod_position: u8) -> Result<Self> {
        let mod_type = ModType::from_str(mod_type)?;

        let parsed_sequence = Sequence::from_str(sequence_str)?;

        if mod_position as usize > parsed_sequence.len() - 1 {
            bail!(
                "mod_position {} is out of bounds for sequence of length {}. Note mod_position is 0-indexed.",
                mod_position,
                parsed_sequence.len()
            );
        }

        let base_at_position = &parsed_sequence[mod_position as usize];
        match mod_type {
            ModType::SixMA => {
                if *base_at_position != IupacBase::A {
                    bail!(
                        "mod_position {} points to base '{}' which is invalid for 6mA.",
                        mod_position,
                        base_at_position
                    );
                }
            }
            ModType::FiveMC | ModType::FourMC => {
                if *base_at_position != IupacBase::C {
                    bail!(
                        "mod_position {} points to base '{}' which is invalid for {} modification type.",
                        mod_position,
                        base_at_position,
                        mod_type
                    );
                }
            }
        }

        if parsed_sequence.first() == Some(&IupacBase::N)
            || parsed_sequence.last() == Some(&IupacBase::N)
        {
            bail!(
                "Motif sequence starts or ends with N, which is invalid: {}",
                sequence_str
            );
        }

        Ok(Self {
            sequence: parsed_sequence,
            mod_type,
            mod_position,
        })
    }

    /// Returns the reverse complement of the motif.
    ///
    /// The reverse complement reverses the sequence and replaces each base
    /// with its complement (e.g., A ↔ T, C ↔ G). The modification position
    /// is adjusted to reflect its position in the reverse-complemented sequence.
    ///
    /// # Examples
    /// ```
    /// use methylome::Motif;
    ///
    /// let motif = Motif::new("TCCCG", "m", 1).unwrap();
    /// let rev_comp = motif.reverse_complement();
    /// assert_eq!(rev_comp.sequence_to_string(), "CGGGA");
    /// assert_eq!(rev_comp.mod_position, 3);
    /// ```
    pub fn reverse_complement(&self) -> Self {
        Self {
            // sequence: (&self.sequence.chars().rev().collect::<String>()).to_string(),
            sequence: self
                .sequence
                .iter()
                .rev()
                .map(IupacBase::to_complement_base)
                .collect(),
            mod_type: self.mod_type.clone(),
            mod_position: self.sequence.len() as u8 - self.mod_position - 1,
        }
    }

    /// Converts the motif sequence into a regular expression string.
    ///
    /// Each base in the sequence is mapped to its corresponding regex
    /// pattern based on IUPAC codes. For example, `R` (purine) becomes `[AG]`.
    ///
    /// # Examples
    /// ```
    /// use methylome::Motif;
    ///
    /// let motif = Motif::new("RGATCY", "a", 2).unwrap();
    /// let regex = motif.to_regex();
    /// assert_eq!(regex, "[AG]GATC[CT]");
    /// ```
    pub fn to_regex(&self) -> String {
        self.sequence.iter().map(IupacBase::to_regex).collect()
    }

    /// Converts the motif sequence into a plain string representation.
    ///
    /// This method maps each IUPAC base in the sequence to its corresponding character.
    ///
    /// # Examples
    /// ```
    /// use methylome::Motif;
    ///
    /// let motif = Motif::new("GATC", "m", 3).unwrap();
    /// let sequence = motif.sequence_to_string();
    /// assert_eq!(sequence, "GATC");
    /// ```
    pub fn sequence_to_string(&self) -> String {
        self.sequence.iter().map(IupacBase::to_string).collect()
    }

    /// Converts a motif sequence into its possible raw DNA sequences.
    ///
    /// Degenerate bases are converted to its corresponding actual nucleotides
    /// and a vector stores all the possible DNA sequences.
    ///
    /// # Examples
    /// ```
    /// use methylome::{IupacBase, Motif};
    ///
    /// let motif = Motif::new("GATCY", "m", 3).unwrap();
    /// let sequences = motif.possible_dna_sequences();
    /// let possible_seq = vec![
    ///     vec![IupacBase::G, IupacBase::A, IupacBase::T, IupacBase::C, IupacBase::C],
    ///     vec![IupacBase::G, IupacBase::A, IupacBase::T, IupacBase::C, IupacBase::T]
    /// ];
    /// assert_eq!(sequences, possible_seq);
    /// ```
    pub fn possible_dna_sequences(&self) -> Vec<Vec<IupacBase>> {
        let mut sequences = vec![Vec::new()];

        for base in &self.sequence {
            let nucleotides = base.to_possible_nucleotides();
            let mut new_sequences = Vec::new();
            for seq in &sequences {
                for nuc in &nucleotides {
                    let mut new_seq = seq.clone();
                    new_seq.push(*nuc);
                    new_sequences.push(new_seq);
                }
            }
            sequences = new_sequences;
        }
        sequences
    }

    /// Checks if current motif is the parent motif of another motif
    ///
    /// # Examples
    /// ```
    /// use methylome::{IupacBase, Motif};
    ///
    /// let parent = Motif::new("GATC", "a", 1).unwrap();
    /// let child = Motif::new("RGATCY", "a", 2).unwrap();
    ///
    /// let parent2 = Motif::new("CCNGG", "m", 0).unwrap();
    /// let not_child = Motif::new("CCWGG", "m", 1).unwrap();
    ///
    /// assert!(parent.is_child_motif(&child));
    /// assert!(!parent2.is_child_motif(&not_child));
    /// ```
    pub fn is_child_motif(&self, child: &Motif) -> bool {
        if self.mod_type != child.mod_type {
            return false;
        };
        if self.sequence.len() > child.sequence.len() {
            return false;
        }

        // The offset that aligns the two modification sites
        let mod_offset = child.mod_position as isize - self.mod_position as isize;
        if mod_offset < 0 {
            return false;
        } // parent (self) would stick out of child.

        if mod_offset + self.sequence.len() as isize > child.sequence.len() as isize {
            return false;
        }

        self.sequence
            .iter()
            .zip(child.sequence[(mod_offset as usize)..].iter())
            .all(|(p, c)| p.mask() & c.mask() != 0)
    }

    /// Extend motif with N's
    ///
    /// # Examples
    /// ```
    /// use methylome::{IupacBase, Motif};
    ///
    /// let mut motif = Motif::new("GATC", "a", 1).unwrap();
    ///
    /// motif.extend_motif_with_n(2);
    /// assert_eq!(motif.sequence_to_string(), "GATCNN");
    /// assert_eq!(motif.mod_position, 1);
    /// ```
    pub fn extend_motif_with_n(&mut self, n: usize) -> &mut Self {
        self.sequence
            .extend(std::iter::repeat(IupacBase::N).take(n));
        self
    }
    /// Extend motif with N's
    ///
    /// # Examples
    /// ```
    /// use methylome::{IupacBase, Motif};
    ///
    /// let mut motif = Motif::new("GATC", "a", 1).unwrap();
    ///
    /// motif.prepend_n(2);
    /// assert_eq!(motif.sequence_to_string(), "NNGATC");
    /// assert_eq!(motif.mod_position, 3);
    /// ```
    pub fn prepend_n(&mut self, n: usize) -> &mut Self {
        let ns = vec![IupacBase::N; n];

        self.sequence.splice(0..0, ns.iter().cloned());
        self.mod_position = self.mod_position + n as u8;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_motif_creation() {
        let motif = Motif::new("GATC", "a", 1).unwrap();
        assert_eq!(motif.sequence, Sequence::from_str("GATC").unwrap());
        assert_eq!(motif.mod_type, ModType::SixMA);
        assert_eq!(motif.mod_position, 1);
    }

    #[test]
    fn test_out_of_bounds() {
        let result = Motif::new("GATC", "m", 4);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "mod_position 4 is out of bounds for sequence of length 4. Note mod_position is 0-indexed."
        );
    }

    #[test]
    fn test_unidentified_motif_type() {
        let result = Motif::new("GATC", "d", 1);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().to_string(), "Unsupported mod type: d");
    }

    #[test]
    fn test_invalid_mod_position_base() {
        let result = Motif::new("ATCG", "m", 3); // 'G' is invalid for 5mC
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "mod_position 3 points to base 'G' which is invalid for 5mC (m) modification type."
        );
    }

    #[test]
    fn test_motif_starts_with_n() {
        let result = Motif::new("NATGC", "m", 4);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Motif sequence starts or ends with N, which is invalid: NATGC",
        );
    }
    #[test]
    fn test_motif_ends_with_n() {
        let result = Motif::new("ATGCN", "m", 3);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Motif sequence starts or ends with N, which is invalid: ATGCN",
        );
    }

    #[test]
    fn test_invalid_iupac_base() {
        let result = Motif::new("ATZG", "a", 0); // 'G' is invalid for 5mC
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Base 'Z' in sequence 'ATZG' is not a valid IUPAC code"
        );
    }

    #[test]
    fn test_motif_reverse_complement() {
        let motif1 = Motif::new("GATC", "m", 3).unwrap();
        let motif2 = Motif::new("TCCCG", "m", 1).unwrap();
        let motif3 = Motif::new("RGATCY", "a", 2).unwrap();
        assert_eq!(
            motif1.reverse_complement().sequence,
            Sequence::from_str("GATC").unwrap()
        );
        assert_eq!(
            motif2.reverse_complement().sequence,
            Sequence::from_str("CGGGA").unwrap()
        );
        assert_eq!(
            motif3.reverse_complement().sequence,
            Sequence::from_str("RGATCY").unwrap()
        );
        assert_eq!(
            motif1.reverse_complement().mod_type,
            ModType::from_str("m").unwrap()
        );
        assert_eq!(
            motif2.reverse_complement().mod_type,
            ModType::from_str("m").unwrap()
        );
        assert_eq!(
            motif3.reverse_complement().mod_type,
            ModType::from_str("a").unwrap()
        );
        assert_eq!(motif1.reverse_complement().mod_position, 0);
        assert_eq!(motif2.reverse_complement().mod_position, 3);
        assert_eq!(motif3.reverse_complement().mod_position, 3);
    }

    #[test]
    fn test_to_regex() {
        let motif1 = Motif::new("GATC", "m", 3).unwrap();
        let motif2 = Motif::new("RGATCY", "m", 4).unwrap();

        assert_eq!(motif1.to_regex(), "GATC");
        assert_eq!(motif2.to_regex(), "[AG]GATC[CT]");
    }

    #[test]
    fn test_is_child_motif() {
        let parent = Motif::new("GATC", "m", 3).unwrap();
        let child = Motif::new("RGATCY", "m", 4).unwrap();

        assert!(parent.is_child_motif(&child));
        assert!(!child.is_child_motif(&parent));
    }
}
