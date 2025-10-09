use anyhow::anyhow;
use methylome::{ModType, Strand};
use std::{fmt, str::FromStr};

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

#[derive(Debug, Clone)]
#[cfg_attr(feature = "python", pyo3::pyclass)]
pub enum PileupColumn {
    Contig,
    Start,
    End,
    ModType,
    Score,
    Strand,
    StartPos,
    EndPos,
    Color,
    NValidCov,
    FractionModified,
    NModified,
    NCanonical,
    NOtherMod,
    NDelete,
    NFail,
    NDiff,
    NNoCall,
}

impl FromStr for PileupColumn {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "contig" => Ok(PileupColumn::Contig),
            "start" => Ok(PileupColumn::Start),
            "end" => Ok(PileupColumn::End),
            "mod_type" => Ok(PileupColumn::ModType),
            "score" => Ok(PileupColumn::Score),
            "strand" => Ok(PileupColumn::Strand),
            "start_pos" => Ok(PileupColumn::StartPos),
            "end_pos" => Ok(PileupColumn::EndPos),
            "color" => Ok(PileupColumn::Color),
            "n_valid_cov" => Ok(PileupColumn::NValidCov),
            "fraction_modified" => Ok(PileupColumn::FractionModified),
            "n_modified" => Ok(PileupColumn::NModified),
            "n_canonical" => Ok(PileupColumn::NCanonical),
            "n_other_mod" => Ok(PileupColumn::NOtherMod),
            "n_delete" => Ok(PileupColumn::NDelete),
            "n_fail" => Ok(PileupColumn::NFail),
            "n_diff" => Ok(PileupColumn::NDiff),
            "n_no_call" => Ok(PileupColumn::NNoCall),
            _ => Err(anyhow!(
                "Could not convert '{}' to pileup column",
                s.to_string()
            )),
        }
    }
}

impl ToString for PileupColumn {
    fn to_string(&self) -> String {
        match self {
            PileupColumn::Contig => "contig".to_string(),
            PileupColumn::Start => "start".to_string(),
            PileupColumn::End => "end".to_string(),
            PileupColumn::ModType => "mod_type".to_string(),
            PileupColumn::Score => "score".to_string(),
            PileupColumn::Strand => "strand".to_string(),
            PileupColumn::StartPos => "start_pos".to_string(),
            PileupColumn::EndPos => "end_pos".to_string(),
            PileupColumn::Color => "color".to_string(),
            PileupColumn::NValidCov => "n_valid_cov".to_string(),
            PileupColumn::FractionModified => "fraction_modified".to_string(),
            PileupColumn::NModified => "n_modified".to_string(),
            PileupColumn::NCanonical => "n_canonical".to_string(),
            PileupColumn::NOtherMod => "n_other_mod".to_string(),
            PileupColumn::NDelete => "n_delete".to_string(),
            PileupColumn::NFail => "n_fail".to_string(),
            PileupColumn::NDiff => "n_diff".to_string(),
            PileupColumn::NNoCall => "n_no_call".to_string(),
        }
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
