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

    use crate::read::{MethBase, MethQual, Read};
    use noodles_fastq::{self as fastq, record::Definition};

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

    #[test]
    fn test_find_motif_indices_in_read() {
        let description = "MM:Z:A+a.,0,0,0,0,0,2,0,9,0,0,0,0,1,0,0,0,0,2,0,0,0,0,16,0,0,0,4,0,0,0,1,11,1,0,1,0,0,0,0,0,4,0,0,0,2,0,10,6,5,11,0,11,1,6,0,0,0,0,0,2,3,12,0,4,16,0,0,1,0,1,4,0,0,0,0,0;C+21839.,6,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,2,11,0,0,0,0,6,0,5,4,2,9,0,1,3,0,0,0,5,2,1,11,1,0,3,0,0;C+m.,6,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,2,11,0,0,0,0,6,0,5,4,2,9,0,1,3,0,0,0,5,2,1,11,1,0,3,0,0; ML:B:C,204,119,22,36,26,40,16,20,15,25,97,104,150,20,112,20,16,34,81,66,52,12,30,67,20,155,15,21,28,20,85,22,13,14,13,19,13,17,24,12,12,14,30,13,20,20,147,16,17,22,36,41,37,163,29,14,71,28,58,12,12,14,14,12,12,15,64,25,137,42,19,34,29,23,231,46,6,16,17,30,9,41,40,25,27,26,14,179,86,24,8,23,42,15,48,12,16,13,15,14,10,16,162,21,9,3,16,14,8,31,3,2,7,4,6,21,3,15,12,19,20,12,83,45,12,18,10,26,17,33,68,70,49,53,23,13,23,21,48,40,83,5,5,5,5,2,11,13,29,60,7,24,12,16,2,3,14,14,44,12,20,13,13,9,14,6,10,6,7,4,2,34";
        let sequence = "ACTATAAATCATTTATTTTATATTTAATGTAAACATTTCTTCACCTTCTAAGGTGCCACAAAGATAATCATTAGCATCTACCCGTCCTACACCTGCTGGTGTACCTGTAAATATGACATCTCCTTTTTTTAACATAAAGTATTGAGAAACATACGATATAAGCTCGTCTATTTTCCACAACATTAAATTTGTATTTCCCTTTTGTACTATTTCTTCATTTTTTAACAATGAAAAATTTATATTATCTACAGAAGAAAATTTACTTTTTGGCAACCATTTACCTATTACTGCCGCACCATCAAAACCTTTTGCTTTTTCCCAAGGTAATCCTTTTTCTTTTAACTTAGATTGAAGATCACGTGCTGTAAAGTCTATTCCTAGACCAATTTCATCATAATAATTGGCAGCAAATTTTTGTTCTATATGTTTTCCTACTTTTTTAATTTTAACTAAAACCTCTACTTCATAATGTATGTTATTAGAGAATTCTGGTATGTAAAAATCTTGTTCTTTTGGTAGCACAGCAGAATCTGGCTTAATAAAAACAACAGGATCTGTAGGCTTTTCATTGGCCAATTCTTTAATATGATCTGTGTAATTACGACCTATACAAATTATCTTCATAATAAACTTTATTAGGTTATTGCTTTTCTTTAGACAATTTTAAACTTACCGCTGGTACTAAATGCTTGGTATTATCTTTAGACTTTACCTTAAATCTATCAT";
        let quals = "JIGSSGIH=GLLMLJFB;=>@?HSSNKJSSMKIKDCD?>??<>?JMSJJGGFSSLLJSSSOIJLSONMJSSRSSJSSSSLSLIIHFGSQSLSSSLSSRSSSSKSQMIHSSISSSSLIQSLLK88:1111EDMQIIOROPMS::22336C@CICCBCDSSNJJEG?;;:<;<>@;8B44HSSSSSPPISHHILGKJQAMLRPJJ:SLMSMQS99999KHKISSSSMJIKSIMSSSSSSSSSIKGJJMSHM50))))GSMSLSS>;;<PSSSSJMPMSSSSNLSSSKIKNSSNKKKONMJSSQLSSB@PKHGJFBSSSJSSQDGSSSQKKSSSSSSSSMMNJSSSSSSSSSSORISSSSSSSSSSSOSSSSSSMSKSSNSRQSSSSSSSSNISSSSSNSSSSMJKNSSQMMPSSSSSSSSSPSNSSSOJJMSSSSJJMKSSSSSSSSSSSSSSSSSSSSSQSSSSSNSSSSSOSNMSHSGFIDMIIFFIG@@ILSMINSLKSSSSLGSBSNSRPSSSSSPPSFCA@@CLSSSSSKNSSHEEEHHSSSSRNLLOSSSSSSSPSSSNSSSSSSRLOLILSSSSSSSQPSSSHJMKJSSSSRSSSIINSLNSRSSSMLOCCHILKJKJKLOSSSSSSSSSSSPSSSSPSSMIAGSA>EEE;AAABBI@>>==@SLLKIKEIHLLLOJHHKMSHISSHQSKQNSPSJFD>@>>>950///.-./(,*((&&&";
        let fastq = fastq::Record::new(
            Definition::new("@8c32d39d-6b88-480f-ac6d-e8061e4d674b", description),
            sequence,
            quals,
        );

        let read = Read::from_fastq_record(fastq).unwrap();
        let motif = Motif::new("ACTATA", "a", 0).unwrap();

        let indices = find_motif_indices_in_sequence(read.get_sequence(), &motif);
        assert_eq!(indices, vec![0]);
    }
    #[test]
    fn test_read_construction_simple() {
        let description = "MM:Z:A+a.,0,1; ML:B:C,204,255";
        let sequence = "GGGCGGATCAGATC";
        let quals = "JIGSSGIH=";
        let fastq = fastq::Record::new(
            Definition::new("8c32d39d-6b88-480f-ac6d-e8061e4d674b", description),
            sequence,
            quals,
        );

        let read = Read::from_fastq_record(fastq).unwrap();
        let modifcations = read.get_modifications();
        let motif = Motif::new("GATC", "a", 1).unwrap();

        assert_eq!(
            *modifcations.0.get(&6).unwrap(),
            MethBase {
                base: ModType::SixMA,
                quality: MethQual(204)
            }
        );
        assert_eq!(modifcations.0.get(&9), None);
        assert_eq!(
            *modifcations.0.get(&11).unwrap(),
            MethBase {
                base: ModType::SixMA,
                quality: MethQual(255)
            }
        );
        let indices = find_motif_indices_in_sequence(read.get_sequence(), &motif);
        assert_eq!(indices, vec![6, 11]);
    }
}
