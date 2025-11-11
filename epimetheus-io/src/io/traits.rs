use std::path::Path;

use ahash::AHashMap;
use anyhow::Result;
use epimetheus_core::models::{contig::Contig, pileup::PileupRecordString};
use methylome::read::Read;

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
    fn read_fasta(
        path: &Path,
        contig_filter: Option<Vec<String>>,
    ) -> Result<AHashMap<String, Contig>>;
}

pub trait FastqReader {
    fn read_fastq(path: &Path, read_filter: Option<Vec<String>>) -> Result<Vec<Read>>;
}
