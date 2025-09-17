use std::path::Path;

use ahash::AHashMap;
use anyhow::Result;

use crate::models::{contig::Contig, pileup::PileupRecordString};

pub trait BatchLoader<T> {
    fn new(
        reader: std::io::BufReader<std::fs::File>,
        assembly: AHashMap<String, Contig>,
        batch_size: usize,
        min_valid_read_coverage: u32,
        min_valid_cov_to_diff_fraction: f32,
        allow_mismatch: bool,
    ) -> Self;
    fn next_batch(&mut self) -> Option<Result<T>>;
}

pub trait PileupReader {
    fn from_path(path: &Path) -> Result<Self>
    where
        Self: Sized;
    fn query_contig(&mut self, contig: &str) -> Result<Vec<PileupRecordString>>;
    fn available_contigs(&self) -> Vec<String>;
}

impl PileupReader for Box<dyn PileupReader> {
    fn from_path(_path: &Path) -> Result<Self>
    where
        Self: Sized,
    {
        unimplemented!("Cannot create Box<dyn PileupReader> from path. Use concrete type.")
    }

    fn query_contig(&mut self, contig: &str) -> Result<Vec<PileupRecordString>> {
        (**self).query_contig(contig)
    }

    fn available_contigs(&self) -> Vec<String> {
        (**self).available_contigs()
    }
}

pub trait FastaReader {
    fn read_fasta(path: &Path) -> Result<AHashMap<String, Contig>>;
}
