use methylome::{ModType, Strand};
use std::fmt;

#[cfg(feature = "pyo3")]
use pyo3::prelude::*;

// pub struct Pileup {
//     records: Vec<PileupRecord>,
// }

// impl Pileup {
//     pub fn new(records: Vec<PileupRecord>) -> Self {
//         Self { records }
//     }
// }

pub struct PileupRecordString(pub String);

impl PileupRecordString {
    pub fn new(_0: String) -> Self {
        Self(_0)
    }
}

#[cfg_attr(feature = "pyo3", pyclass)]
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

#[cfg(feature = "pyo3")]
#[pymethods]
impl PileupRecord {
    #[getter]
    fn contig(&self) -> &str {
        &self.contig
    }
    #[getter]
    fn start(&self) -> u32 {
        self.start
    }
    #[getter]
    fn end(&self) -> u32 {
        self.end
    }
    #[getter]
    fn mod_type(&self) -> String {
        self.mod_type.to_pileup_code().to_string()
    }
    #[getter]
    fn score(&self) -> u32 {
        self.score
    }
    #[getter]
    fn strand(&self) -> String {
        self.strand.to_string()
    }
    #[getter]
    fn start_pos(&self) -> u32 {
        self.start_pos
    }
    #[getter]
    fn end_pos(&self) -> u32 {
        self.end_pos
    }
    #[getter]
    fn color(&self) -> &str {
        &self.color
    }
    #[getter]
    fn n_valid_cov(&self) -> u32 {
        self.n_valid_cov
    }
    #[getter]
    fn fraction_modified(&self) -> f64 {
        self.fraction_modified
    }
    #[getter]
    fn n_modified(&self) -> u32 {
        self.n_modified
    }
    #[getter]
    fn n_canonical(&self) -> u32 {
        self.n_canonical
    }
    #[getter]
    fn n_other_mod(&self) -> u32 {
        self.n_other_mod
    }
    #[getter]
    fn n_delete(&self) -> u32 {
        self.n_delete
    }
    #[getter]
    fn n_fail(&self) -> u32 {
        self.n_fail
    }
    #[getter]
    fn n_diff(&self) -> u32 {
        self.n_diff
    }
    #[getter]
    fn n_no_call(&self) -> u32 {
        self.n_no_call
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
