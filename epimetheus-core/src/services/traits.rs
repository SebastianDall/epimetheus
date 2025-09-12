use std::path::Path;

use anyhow::Result;

use crate::models::pileup::PileupRecordString;

pub trait BatchLoader<T> {
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
