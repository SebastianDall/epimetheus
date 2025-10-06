use methylome::{ModType, Strand};
use std::fmt;

use crate::models::methylation::{MethylationCoverage, MethylationRecord};

// pub struct Pileup {
//     records: Vec<PileupRecord>,
// }

// impl Pileup {
//     pub fn new(records: Vec<PileupRecord>) -> Self {
//         Self { records }
//     }
// }

#[derive(Clone)]
pub struct PileupRecordString(pub String);

impl PileupRecordString {
    pub fn new(_0: String) -> Self {
        Self(_0)
    }
}

#[derive(Clone)]
pub struct PileupRecord {
    pub contig: String,
    pub start: u32,
    pub end: u32,
    pub mod_type: ModType,
    pub score: u32,
    pub strand: Strand,
    pub start_pos: u32,
    pub end_pos: u32,
    pub color: String,
    pub n_valid_cov: u32,
    pub fraction_modified: f64,
    pub n_modified: u32,
    pub n_canonical: u32,
    pub n_other_mod: u32,
    pub n_delete: u32,
    pub n_fail: u32,
    pub n_diff: u32,
    pub n_no_call: u32,
}

impl PileupRecord {
    pub fn new(
        contig: String,
        start: u32,
        end: u32,
        mod_type: ModType,
        score: u32,
        strand: Strand,
        start_pos: u32,
        end_pos: u32,
        color: String,
        n_valid_cov: u32,
        fraction_modified: f64,
        n_modified: u32,
        n_canonical: u32,
        n_other_mod: u32,
        n_delete: u32,
        n_fail: u32,
        n_diff: u32,
        n_no_call: u32,
    ) -> Self {
        Self {
            contig,
            start,
            end,
            mod_type,
            score,
            strand,
            start_pos,
            end_pos,
            color,
            n_valid_cov,
            fraction_modified,
            n_modified,
            n_canonical,
            n_other_mod,
            n_delete,
            n_fail,
            n_diff,
            n_no_call,
        }
    }

    pub fn to_methylation_record(&self) -> anyhow::Result<MethylationRecord> {
        let methylation_coverage =
            MethylationCoverage::new(self.n_modified, self.n_valid_cov, self.n_other_mod)?;

        Ok(MethylationRecord {
            contig: self.contig.clone(),
            position: self.start as usize,
            strand: self.strand,
            mod_type: self.mod_type,
            methylation: methylation_coverage,
        })
    }
}

impl TryFrom<PileupRecordString> for PileupRecord {
    type Error = anyhow::Error;

    fn try_from(value: PileupRecordString) -> std::result::Result<Self, Self::Error> {
        let fields: Vec<&str> = value.0.trim().split('\t').collect();

        Ok(Self {
            contig: fields[0].to_string(),
            start: fields[1].parse()?,
            end: fields[2].parse()?,
            mod_type: fields[3].parse()?,
            score: fields[4].parse()?,
            strand: fields[5].parse()?,
            start_pos: fields[6].parse()?,
            end_pos: fields[7].parse()?,
            color: fields[8].to_string(),
            n_valid_cov: fields[9].parse()?,
            fraction_modified: fields[10].parse()?,
            n_modified: fields[11].parse()?,
            n_canonical: fields[12].parse()?,
            n_other_mod: fields[13].parse()?,
            n_delete: fields[14].parse()?,
            n_fail: fields[15].parse()?,
            n_diff: fields[16].parse()?,
            n_no_call: fields[17].parse()?,
        })
    }
}

impl fmt::Display for PileupRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.contig,
            self.start,
            self.end,
            self.mod_type.to_pileup_code(),
            self.score,
            self.strand,
            self.start_pos,
            self.end_pos,
            self.color,
            self.n_valid_cov,
            self.fraction_modified,
            self.n_modified,
            self.n_canonical,
            self.n_other_mod,
            self.n_delete,
            self.n_fail,
            self.n_diff,
            self.n_no_call,
        )
    }
}
