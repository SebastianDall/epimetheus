// use regex::Regex;

pub mod iupac;
pub mod modtype;
pub mod motif;
pub mod read;
pub mod sequence;
pub mod strand;

pub use iupac::IupacBase;
pub use modtype::ModType;
pub use motif::Motif;
pub use strand::Strand;

use crate::sequence::Sequence;

pub fn find_motif_indices_in_sequence(sequence: &Sequence, motif: &Motif) -> Vec<usize> {
    // let regex_str = motif.to_regex();
    // let re = Regex::new(&regex_str).expect("Expected regex pattern");

    // let indices = re
    //     .find_iter(sequence)
    //     .map(|m| m.start() as usize + motif.mod_position as usize)
    //     .collect();

    let motif_bases = motif.sequence.clone();
    let motif_len = motif_bases.len();
    let mut indices = Vec::new();

    if sequence.len() < motif_len {
        return indices;
    }

    for i in 0..=(sequence.len() - motif_len) {
        let mut matches = true;

        for (j, &motif_base) in motif_bases.iter().enumerate() {
            let seq_base = sequence[i + j];
            if (seq_base.mask() & motif_base.mask()) == 0 {
                matches = false;
                break;
            }
        }

        if matches {
            indices.push(i + motif.mod_position as usize);
        }
    }

    indices
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_find_motif_indices_in_contig() {
        let contig = Sequence::from_str("GGATCTCCATGATC").unwrap();
        let contig2 = Sequence::from_str("TGGACGATCCCGATC").unwrap();
        let motif1 = Motif::new("GATC", "m", 3).unwrap();
        let motif2 = Motif::new("RGATCY", "m", 4).unwrap();
        let motif3 = Motif::new("GATC", "a", 1).unwrap();
        let motif4 = Motif::new("GGANNNTCC", "a", 2).unwrap();

        println!("{}", &motif4.to_regex());
        assert_eq!(
            find_motif_indices_in_sequence(&contig, &motif1),
            vec![4, 13]
        );
        assert_eq!(find_motif_indices_in_sequence(&contig, &motif2), vec![4]);

        assert_eq!(
            find_motif_indices_in_sequence(&contig2, &motif3),
            vec![6, 12]
        );
        assert_eq!(
            find_motif_indices_in_sequence(&contig2, &motif3.reverse_complement()),
            vec![7, 13]
        );

        assert_eq!(find_motif_indices_in_sequence(&contig2, &motif4), vec![3])
    }
}
