use anyhow::{Result, bail};
use humantime::format_duration;
use log::{debug, error, info};
use methylome::Motif;
use std::time::Instant;

use crate::{
    algorithms::methylation_pattern::calculate_contig_read_methylation_pattern,
    models::{genome_workspace::GenomeWorkspace, methylation::MotifMethylationDegree},
    services::traits::BatchLoader,
};

pub fn sequential_processer<L: BatchLoader<GenomeWorkspace>>(
    loader: &mut L,
    motifs: Vec<Motif>,
    threads: usize,
) -> Result<Vec<MotifMethylationDegree>> {
    let mut methylation_pattern_results: Vec<MotifMethylationDegree> = Vec::new();

    let mut batch_processing_time = Instant::now();
    let mut contigs_processed = 0;
    loop {
        match loader.next_batch() {
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
