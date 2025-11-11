use anyhow::{Result, anyhow};
use std::collections::HashMap;

use crate::{IupacBase, ModType, sequence::Sequence};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MethQual(pub u8);

impl MethQual {
    pub fn new(_0: u8) -> Self {
        Self(_0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MethBase {
    pub base: ModType,
    pub quality: MethQual,
}

impl MethBase {
    pub fn new(base: ModType, quality: MethQual) -> Self {
        Self { base, quality }
    }
}

type Position = usize;

#[derive(Debug)]
pub struct BaseModifications(pub HashMap<Position, MethBase>);

impl BaseModifications {
    pub fn new() -> Self {
        Self(HashMap::new())
    }
}

pub type ReadId = String;

#[derive(Debug)]
pub struct Read {
    name: ReadId,
    sequence: Sequence,
    modifications: BaseModifications,
}

impl Read {
    pub fn new(name: ReadId, sequence: Sequence, modifications: BaseModifications) -> Self {
        Self {
            name,
            sequence,
            modifications,
        }
    }

    pub fn from_fastq_record(record: noodles_fastq::Record) -> Result<Self> {
        let sequence: Result<Vec<IupacBase>, String> = record
            .sequence()
            .iter()
            .enumerate()
            .map(|(i, &byte)| {
                IupacBase::from_ascii(byte)
                    .ok_or_else(|| format!("Invalid base '{}' at position {}", byte as char, i))
            })
            .collect();

        let sequence = Sequence::from_iupac(sequence.map_err(|e| anyhow!("{}", e.to_string()))?);

        let description = record.description().to_string();
        let tags = parse_header_tags(&description);

        let mm_string = tags.get("MM").unwrap_or(&"".to_string()).clone();
        let ml_string = tags.get("ML").unwrap_or(&"".to_string()).clone();
        let quality_scores = parse_ml_records(&ml_string)?;

        let skip_distances = MethSkipDistances::from_meth_tags(mm_string, quality_scores)?;

        let mods = convert_skip_distances_to_positions(&sequence, skip_distances)?;

        Ok(Self::new(record.name().to_string(), sequence, mods))
    }

    pub fn get_name(&self) -> &String {
        &self.name
    }
    pub fn get_sequence(&self) -> &Sequence {
        &self.sequence
    }
    pub fn get_modifications(&self) -> &BaseModifications {
        &self.modifications
    }
}

pub fn convert_skip_distances_to_positions(
    seq: &Sequence,
    skip_distances: MethSkipDistances,
) -> Result<BaseModifications> {
    let mut base_mods = BaseModifications::new();

    for (mod_type, skips_with_qual) in skip_distances.distances {
        let target_base = IupacBase::from_mod_type(&mod_type);

        let target_positions: Vec<usize> = seq
            .0
            .iter()
            .enumerate()
            .filter_map(|(pos, base)| {
                if *base == target_base {
                    Some(pos)
                } else {
                    None
                }
            })
            .collect();

        let mut target_index = 0;
        for (skip_distance, quality) in skips_with_qual {
            target_index += skip_distance.0;

            if let Some(&seq_pos) = target_positions.get(target_index) {
                let new_meth_base = MethBase {
                    base: mod_type,
                    quality,
                };

                match base_mods.0.get(&seq_pos) {
                    Some(existing) => {
                        if new_meth_base.quality.0 > existing.quality.0 {
                            base_mods.0.insert(seq_pos, new_meth_base);
                        } else if new_meth_base.quality.0 == existing.quality.0 {
                            base_mods.0.remove(&seq_pos);
                        }
                    }
                    None => {
                        base_mods.0.insert(seq_pos, new_meth_base);
                    }
                }
                target_index += 1;
            }
        }
    }

    Ok(base_mods)
}

fn parse_header_tags(header: &str) -> HashMap<String, String> {
    let mut tags = HashMap::new();

    for part in header.split_whitespace() {
        if let Some((tag, value)) = part.split_once(':') {
            if let Some((_, actual_value)) = value.split_once(':') {
                // MM:Z:A+a.,0,0... -> ("MM", "A+a.,0,0...")
                tags.insert(tag.to_string(), actual_value.to_string());
            }
        }
    }

    tags
}

#[derive(Debug)]
pub struct SkipDistance(pub usize);
#[derive(Debug)]
pub struct MethSkipDistances {
    pub distances: HashMap<ModType, Vec<(SkipDistance, MethQual)>>,
}

impl MethSkipDistances {
    /// Parses the methylation skip distances from the Nanopore FASTQ header
    ///
    /// # Examples
    /// ```
    /// use methylome::read::{MethSkipDistances, SkipDistance, MethQual};
    /// use methylome::ModType;
    /// use noodles_fastq::record::Definition;
    ///
    /// let definition = Definition::new("read-id", "MM:Z:A+a.,0,1;C+m.,2; ML:B:C,255,204,180");
    /// let mm_string = "A+a.,0,1;C+m.,2;".to_string();
    /// let quality_scores = vec![MethQual::new(255), MethQual::new(204), MethQual::new(180)];
    /// let distances = MethSkipDistances::from_meth_tags(mm_string, quality_scores).unwrap();
    ///
    /// // Check that we parsed SixMA modifications correctly
    /// let sixma = distances.distances.get(&ModType::SixMA).unwrap();
    /// assert_eq!(sixma.len(), 2);
    /// assert_eq!(sixma[0].0.0, 0);   // First skip distance: 0
    /// assert_eq!(sixma[0].1.0, 255); // First quality: 255
    /// assert_eq!(sixma[1].0.0, 1);   // Second skip distance: 1
    /// assert_eq!(sixma[1].1.0, 204); // Second quality: 204
    /// ```
    pub fn from_meth_tags(mm_string: String, quality_scores: Vec<MethQual>) -> Result<Self> {
        let mut distances = HashMap::new();

        if mm_string.chars().collect::<Vec<char>>().len() == 0 {
            return Ok(Self { distances });
        }

        let total_modifications: usize = mm_string
            .split(';')
            .map(|segment| segment.split(',').skip(1).count())
            .sum();

        if quality_scores.len() != total_modifications {
            return Err(anyhow!(
                "MM/ML length mismatch: {} modifications but {} quality scores",
                total_modifications,
                quality_scores.len()
            ));
        }

        let mut quality_iter = quality_scores.into_iter();

        for segment in mm_string.split(";") {
            if let Some((mod_info, distances_str)) = segment.split_once(',') {
                if let Some(mod_type) = parse_mod_type(mod_info) {
                    let skip_distances: Vec<SkipDistance> = distances_str
                        .split(',')
                        .filter_map(|s| s.parse().ok().map(SkipDistance))
                        .collect();

                    let distances_with_qual: Vec<(SkipDistance, MethQual)> = skip_distances
                        .into_iter()
                        .filter_map(|skip| quality_iter.next().map(|qual| (skip, qual)))
                        .collect();
                    distances.insert(mod_type, distances_with_qual);
                }
            }
        }

        Ok(Self { distances })
    }
}

fn parse_mod_type(mod_info: &str) -> Option<ModType> {
    let base = mod_info.chars().nth(0)?;
    let mod_code = &mod_info[2..mod_info.len() - 1];

    ModType::from_sam_code(base, mod_code)
}

fn parse_ml_records(ml_string: &str) -> Result<Vec<MethQual>> {
    let qualities: Result<Vec<MethQual>> = ml_string
        .split(',')
        .skip(1) // Skip the first C
        .map(|s| {
            s.parse::<u8>()
                .map(MethQual)
                .map_err(|e| anyhow!("Parse error: '{}'", e.to_string()))
        })
        .collect();

    Ok(qualities?)
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use noodles_fastq::{self as fastq, record::Definition};

    #[test]
    fn test_parse_header_tags() {
        let header = "@8c32d39d-6b88-480f-ac6d-e8061e4d674b   MM:Z:A+a.,0,0,0,0;C+21839.,6,0,1,1;C+m.,6,0,1,1; ML:B:C,204,119,22,36";

        let tags = parse_header_tags(header);
        println!("{:?}", tags);
        assert_eq!(
            tags.get("MM"),
            Some(&"A+a.,0,0,0,0;C+21839.,6,0,1,1;C+m.,6,0,1,1;".to_string())
        );
        assert_eq!(tags.get("ML"), Some(&"C,204,119,22,36".to_string()));
    }

    #[test]
    fn test_meth_skip_distances() {
        use noodles_fastq::record::Definition;

        let definition = Definition::new(
            "read-id",
            "MM:Z:A+a.,0,1,0;C+m.,3,0,0; ML:B:C,255,255,255,204,255,255",
        );
        let tags = parse_header_tags(&definition.description().to_string());

        let mm_string = tags.get("MM").unwrap_or(&"".to_string()).clone();
        println!("{mm_string}");
        let ml_string = tags.get("ML").unwrap_or(&"".to_string()).clone();
        let quality_scores = parse_ml_records(&ml_string).unwrap();
        let distances = MethSkipDistances::from_meth_tags(mm_string, quality_scores).unwrap();

        // Debug: print what we actually parsed
        println!("Parsed distances: {:?}", distances);
        println!("Keys: {:?}", distances.distances.keys().collect::<Vec<_>>());

        // Check SixMA modifications
        if let Some(sixma_distances) = distances.distances.get(&ModType::SixMA) {
            println!("SixMA distances: {:?}", sixma_distances);
            assert_eq!(sixma_distances.len(), 3);
            assert_eq!(sixma_distances[0].0.0, 0); // First skip distance
            assert_eq!(sixma_distances[0].1.0, 255); // First quality
        } else {
            panic!("No SixMA data found!");
        }

        // Check FiveMC modifications
        if let Some(fivemc_distances) = distances.distances.get(&ModType::FiveMC) {
            println!("FiveMC distances: {:?}", fivemc_distances);
            assert_eq!(fivemc_distances.len(), 3);
            assert_eq!(fivemc_distances[0].0.0, 3); // First skip distance
            assert_eq!(fivemc_distances[0].1.0, 204); // First quality
        } else {
            panic!("No FiveMC data found!");
        }
    }

    #[test]
    fn test_read_construction() {
        let description = "MM:Z:A+a.,0,0,0,0,0,2,0,9,0,0,0,0,1,0,0,0,0,2,0,0,0,0,16,0,0,0,4,0,0,0,1,11,1,0,1,0,0,0,0,0,4,0,0,0,2,0,10,6,5,11,0,11,1,6,0,0,0,0,0,2,3,12,0,4,16,0,0,1,0,1,4,0,0,0,0,0;C+21839.,6,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,2,11,0,0,0,0,6,0,5,4,2,9,0,1,3,0,0,0,5,2,1,11,1,0,3,0,0;C+m.,6,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,2,11,0,0,0,0,6,0,5,4,2,9,0,1,3,0,0,0,5,2,1,11,1,0,3,0,0; ML:B:C,204,119,22,36,26,40,16,20,15,25,97,104,150,20,112,20,16,34,81,66,52,12,30,67,20,155,15,21,28,20,85,22,13,14,13,19,13,17,24,12,12,14,30,13,20,20,147,16,17,22,36,41,37,163,29,14,71,28,58,12,12,14,14,12,12,15,64,25,137,42,19,34,29,23,231,46,6,16,17,30,9,41,40,25,27,26,14,179,86,24,8,23,42,15,48,12,16,13,15,14,10,16,162,21,9,3,16,14,8,31,3,2,7,4,6,21,3,15,12,19,20,12,83,45,12,18,10,26,17,33,68,70,49,53,23,13,23,21,48,40,83,5,5,5,5,2,11,13,29,60,7,24,12,16,2,3,14,14,44,12,20,13,13,9,14,6,10,6,7,4,2,34";
        let sequence = "ACTATAAATCATTTATTTTATATTTAATGTAAACATTTCTTCACCTTCTAAGGTGCCACAAAGATAATCATTAGCATCTACCCGTCCTACACCTGCTGGTGTACCTGTAAATATGACATCTCCTTTTTTTAACATAAAGTATTGAGAAACATACGATATAAGCTCGTCTATTTTCCACAACATTAAATTTGTATTTCCCTTTTGTACTATTTCTTCATTTTTTAACAATGAAAAATTTATATTATCTACAGAAGAAAATTTACTTTTTGGCAACCATTTACCTATTACTGCCGCACCATCAAAACCTTTTGCTTTTTCCCAAGGTAATCCTTTTTCTTTTAACTTAGATTGAAGATCACGTGCTGTAAAGTCTATTCCTAGACCAATTTCATCATAATAATTGGCAGCAAATTTTTGTTCTATATGTTTTCCTACTTTTTTAATTTTAACTAAAACCTCTACTTCATAATGTATGTTATTAGAGAATTCTGGTATGTAAAAATCTTGTTCTTTTGGTAGCACAGCAGAATCTGGCTTAATAAAAACAACAGGATCTGTAGGCTTTTCATTGGCCAATTCTTTAATATGATCTGTGTAATTACGACCTATACAAATTATCTTCATAATAAACTTTATTAGGTTATTGCTTTTCTTTAGACAATTTTAAACTTACCGCTGGTACTAAATGCTTGGTATTATCTTTAGACTTTACCTTAAATCTATCAT";
        let quals = "JIGSSGIH=GLLMLJFB;=>@?HSSNKJSSMKIKDCD?>??<>?JMSJJGGFSSLLJSSSOIJLSONMJSSRSSJSSSSLSLIIHFGSQSLSSSLSSRSSSSKSQMIHSSISSSSLIQSLLK88:1111EDMQIIOROPMS::22336C@CICCBCDSSNJJEG?;;:<;<>@;8B44HSSSSSPPISHHILGKJQAMLRPJJ:SLMSMQS99999KHKISSSSMJIKSIMSSSSSSSSSIKGJJMSHM50))))GSMSLSS>;;<PSSSSJMPMSSSSNLSSSKIKNSSNKKKONMJSSQLSSB@PKHGJFBSSSJSSQDGSSSQKKSSSSSSSSMMNJSSSSSSSSSSORISSSSSSSSSSSOSSSSSSMSKSSNSRQSSSSSSSSNISSSSSNSSSSMJKNSSQMMPSSSSSSSSSPSNSSSOJJMSSSSJJMKSSSSSSSSSSSSSSSSSSSSSQSSSSSNSSSSSOSNMSHSGFIDMIIFFIG@@ILSMINSLKSSSSLGSBSNSRPSSSSSPPSFCA@@CLSSSSSKNSSHEEEHHSSSSRNLLOSSSSSSSPSSSNSSSSSSRLOLILSSSSSSSQPSSSHJMKJSSSSRSSSIINSLNSRSSSMLOCCHILKJKJKLOSSSSSSSSSSSPSSSSPSSMIAGSA>EEE;AAABBI@>>==@SLLKIKEIHLLLOJHHKMSHISSHQSKQNSPSJFD>@>>>950///.-./(,*((&&&";
        let fastq = fastq::Record::new(
            Definition::new("@8c32d39d-6b88-480f-ac6d-e8061e4d674b", description),
            sequence,
            quals,
        );

        let read = Read::from_fastq_record(fastq).unwrap();
        let modifcations = read.get_modifications();
        println!("{:?}", *modifcations.0.get(&47).unwrap());

        assert_eq!(
            *modifcations.0.get(&0).unwrap(),
            MethBase {
                base: ModType::SixMA,
                quality: MethQual(204)
            }
        );
        assert_eq!(
            *modifcations.0.get(&3).unwrap(),
            MethBase {
                base: ModType::SixMA,
                quality: MethQual(119)
            }
        );
        assert_eq!(
            *modifcations.0.get(&47).unwrap(),
            MethBase {
                base: ModType::FiveMC,
                quality: MethQual(18)
            }
        );
    }
}
