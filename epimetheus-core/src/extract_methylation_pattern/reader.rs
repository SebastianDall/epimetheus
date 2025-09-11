use humantime::format_duration;
use log::{debug, error, info};
use methylome::Motif;
use rayon::prelude::*;
use std::{path::Path, time::Instant};

use crate::{
    data::{GenomeWorkspace, contig::Contig},
    extract_methylation_pattern::{
        batch_loader::SequentialBatchLoader, parallel_batch_loader::process_contig,
    },
    processing::{
        MotifMethylationDegree, calculate_contig_read_methylation_pattern,
        calculate_contig_read_methylation_single,
    },
};
use ahash::AHashMap;
use anyhow::{Result, bail};

pub trait BatchReader {
    fn next_batch(&mut self) -> Option<Result<GenomeWorkspace>>;
}

pub fn parallel_processer(
    file: &Path,
    contigs: &AHashMap<String, Contig>,
    motifs: Vec<Motif>,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
) -> Result<Vec<MotifMethylationDegree>> {
    let methylation = contigs
        .par_iter()
        .map(
            |(_contig_id, contig)| -> Result<Vec<MotifMethylationDegree>> {
                let mut reader = epimetheus_support::zipper::reader::PileupReader::from_path(file)?;

                let contig_w_meth = process_contig(
                    &mut reader,
                    contig,
                    min_valid_read_coverage,
                    min_valid_cov_to_diff_fraction,
                )?;

                let meth_pattern =
                    calculate_contig_read_methylation_single(&contig_w_meth, motifs.clone())?;

                Ok(meth_pattern)
            },
        )
        .collect::<Result<Vec<Vec<_>>>>()?
        .into_iter()
        .flatten()
        .collect();

    Ok(methylation)
}

pub fn sequential_processer(
    file: &Path,
    contigs: AHashMap<String, Contig>,
    motifs: Vec<Motif>,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,
    threads: usize,
) -> Result<Vec<MotifMethylationDegree>> {
    let file = std::fs::File::open(file)?;
    let buf_reader = std::io::BufReader::new(file);
    let mut batch_loader = SequentialBatchLoader::new(
        buf_reader,
        contigs,
        100,
        min_valid_read_coverage,
        min_valid_cov_to_diff_fraction,
        allow_mismatch,
    );
    let mut methylation_pattern_results: Vec<MotifMethylationDegree> = Vec::new();

    let mut batch_processing_time = Instant::now();
    let mut contigs_processed = 0;
    loop {
        match batch_loader.next_batch() {
            Some(ws_result) => match ws_result {
                Ok(workspace) => {
                    debug!("Workspace initialized");
                    let contigs_in_batch = workspace.get_workspace().len() as u32;
                    let mut methylation_pattern = calculate_contig_read_methylation_pattern(
                        workspace,
                        motifs.clone(),
                        threads,
                    )?;
                    methylation_pattern_results.append(&mut methylation_pattern);

                    contigs_processed += contigs_in_batch;
                    let elapsed_batch_processing_time = batch_processing_time.elapsed();
                    if contigs_processed % 100 == 0 {
                        info!(
                            "Finished processing {} contigs. Processing time: {}",
                            contigs_processed,
                            format_duration(elapsed_batch_processing_time).to_string()
                        );
                    }
                    batch_processing_time = Instant::now();
                }
                Err(e) => {
                    error!("Error reading batch: {e}");
                    bail!("Processing terminated due to error: {e}")
                }
            },
            None => break,
        }
    }

    Ok(methylation_pattern_results)
}
