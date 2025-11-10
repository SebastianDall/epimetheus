use ahash::AHashMap;
use anyhow::Result;

use crate::models::contig::Contig;

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
