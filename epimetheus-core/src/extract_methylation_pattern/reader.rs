use std::path::Path;

use crate::{
    data::{GenomeWorkspace, contig::Contig},
    extract_methylation_pattern::{
        batch_loader::SequentialBatchLoader, parallel_batch_loader::ParallelBatchLoader,
    },
};
use ahash::AHashMap;
use anyhow::Result;

pub trait BatchReader {
    fn next_batch(&mut self) -> Option<Result<GenomeWorkspace>>;
}

pub fn set_reader_strategy(
    file: &Path,
    assembly: AHashMap<String, Contig>,
    batch_size: usize,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,
    threads: Option<usize>,
) -> Result<Box<dyn BatchReader>> {
    let reader: Box<dyn BatchReader> = if file.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(ParallelBatchLoader::new(
            file,
            assembly,
            batch_size,
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
            allow_mismatch,
            threads.unwrap_or(1),
        )?)
    } else {
        let file = std::fs::File::open(file)?;
        let buf_reader = std::io::BufReader::new(file);
        Box::new(SequentialBatchLoader::new(
            buf_reader,
            assembly,
            batch_size,
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
            allow_mismatch,
        ))
    };
    Ok(reader)
}
