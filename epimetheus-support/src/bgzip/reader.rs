use anyhow::{anyhow, Result};
use epimetheus_core::data::contig::ContigId;

use crate::bgzip::args::BgzipReaderArgs;

#[derive(Debug)]
pub struct ReaderConfig {
    pub input: String,
    pub output: Option<String>,
    pub contigs: Option<Vec<ContigId>>,
}

pub fn pileup_reader(args: &ReaderConfig) -> Result<()> {
    Ok(())
}
