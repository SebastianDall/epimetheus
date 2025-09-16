use ahash::AHashMap;
use anyhow::{Result, bail};
use humantime::format_duration;
use log::{debug, error, info};
use methylome::Motif;
use rayon::prelude::*;
use std::time::Instant;

use crate::{
    algorithms::methylation_pattern::calculate_contig_read_methylation_pattern,
    models::{
        genome_workspace::GenomeWorkspace,
        methylation::{MethylationOutput, MethylationPatternVariant, MotifMethylationPositions},
    },
    services::traits::BatchLoader,
};

pub fn sequential_processer<L: BatchLoader<GenomeWorkspace>>(
    loader: &mut L,
    motifs: Vec<Motif>,
    threads: usize,
    output: &MethylationOutput,
) -> Result<MethylationPatternVariant> {
    let mut methylation_pattern_results: Vec<MethylationPatternVariant> = Vec::new();

    let mut batch_processing_time = Instant::now();
    let mut contigs_processed = 0;
    loop {
        match loader.next_batch() {
            Some(ws_result) => match ws_result {
                Ok(workspace) => {
                    debug!("Workspace initialized");
                    let contigs_in_batch = workspace.get_workspace().len() as u32;
                    let methylation_pattern = calculate_contig_read_methylation_pattern(
                        workspace,
                        motifs.clone(),
                        threads,
                    )?;

                    let merged_results = match output {
                        MethylationOutput::Raw => {
                            MethylationPatternVariant::Raw(methylation_pattern)
                        }
                        MethylationOutput::Median => MethylationPatternVariant::Median(
                            methylation_pattern.to_median_degrees(),
                        ),

                        MethylationOutput::WeightedMean => MethylationPatternVariant::WeightedMean(
                            methylation_pattern.to_weighted_mean_degress(),
                        ),
                    };
                    methylation_pattern_results.push(merged_results);

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

    let merged_results = match output {
        MethylationOutput::Raw => {
            let mut all_results = AHashMap::new();
            for res in methylation_pattern_results {
                if let MethylationPatternVariant::Raw(positions) = res {
                    all_results.extend(positions.methylation);
                }
            }
            MethylationPatternVariant::Raw(MotifMethylationPositions::new(all_results))
        }
        MethylationOutput::Median => {
            let collected = methylation_pattern_results
                .into_par_iter()
                .flat_map(|meth| {
                    if let MethylationPatternVariant::Median(median) = meth {
                        median
                    } else {
                        Vec::new()
                    }
                })
                .collect();

            MethylationPatternVariant::Median(collected)
        }

        MethylationOutput::WeightedMean => {
            let collected = methylation_pattern_results
                .into_par_iter()
                .flat_map(|meth| {
                    if let MethylationPatternVariant::WeightedMean(weighted_mean) = meth {
                        weighted_mean
                    } else {
                        Vec::new()
                    }
                })
                .collect();

            MethylationPatternVariant::WeightedMean(collected)
        }
    };

    Ok(merged_results)
}
