use ahash::AHashMap;
use anyhow::{Result, bail};
use epimetheus_core::{
    algorithms::methylation_pattern::calculate_contig_read_methylation_single,
    models::{
        contig::Contig,
        genome_workspace::GenomeWorkspace,
        methylation::{
            MethylationOutput, MethylationPatternVariant, MethylationRecord,
            MotifMethylationPositions,
        },
    },
    services::{
        domain::contig_service::populate_contig_with_methylation,
        traits::{BatchLoader, PileupReader},
    },
};
use epimetheus_io::services::data_loading_service::load_pileup_records_for_contig;
use humantime::format_duration;
use indicatif::ProgressBar;
use log::{debug, info};
use methylome::Motif;
use rayon::prelude::*;
use std::path::Path;
use std::{collections::HashSet, time::Instant};

pub fn extract_methylation_patten_from_gz<R: PileupReader + Clone>(
    contigs: AHashMap<String, Contig>,
    pileup_path: &Path,
    motifs: Vec<Motif>,
    threads: usize,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,
    output_type: &MethylationOutput,
) -> Result<MethylationPatternVariant> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Could not initialize threadpool");

    let contigs_in_index: HashSet<String> = R::from_path(pileup_path)?
        .available_contigs()
        .into_iter()
        .collect();

    let filtered_contigs: Vec<(&String, &Contig)> = if allow_mismatch {
        contigs
            .iter()
            .filter(|(contig_id, _)| contigs_in_index.contains(*contig_id))
            .collect()
    } else {
        let contig_vec = contigs.iter().collect();
        let missing_in_pileup: Vec<&String> = contigs
            .keys()
            .filter(|contig_id| !contigs_in_index.contains(*contig_id))
            .collect();

        if !missing_in_pileup.is_empty() {
            bail!(
                "Contig mismatch detected between pileup and assembly. Use --allow-mismatch to ignore this error. The following contigs are in the assembly but not the pileup: {:?}",
                missing_in_pileup
            );
        }
        contig_vec
    };

    let progress_bar = ProgressBar::new(filtered_contigs.len() as u64);

    let per_contig_results = filtered_contigs
        .par_iter()
        .map(|(contig_id, contig)| -> Result<MethylationPatternVariant> {
            let pileup_records = load_pileup_records_for_contig::<R>(pileup_path, contig_id)?;
            debug!(
                "{}\nPileup records before filtering: {}",
                contig_id,
                pileup_records.len()
            );

            let mut meth_records = Vec::new();
            for rec in pileup_records {
                let meth = MethylationRecord::try_from_with_filters(
                    rec,
                    min_valid_read_coverage,
                    min_valid_cov_to_diff_fraction,
                )?;

                match meth {
                    Some(m) => meth_records.push(m),
                    None => continue,
                }
            }

            debug!(
                "{}\nMethylation records after filtering: {}",
                contig_id,
                meth_records.len()
            );

            let contig_w_meth = populate_contig_with_methylation(contig, meth_records)?;

            let positions =
                calculate_contig_read_methylation_single(&contig_w_meth, motifs.clone())?;

            progress_bar.inc(1);
            match output_type {
                MethylationOutput::Raw => Ok(MethylationPatternVariant::Raw(positions)),
                MethylationOutput::Median => Ok(MethylationPatternVariant::Median(
                    positions.to_median_degrees(),
                )),
                MethylationOutput::WeightedMean => Ok(MethylationPatternVariant::WeightedMean(
                    positions.to_weighted_mean_degress(),
                )),
            }
        })
        .collect::<Result<Vec<MethylationPatternVariant>>>()?;

    let merged_results = match output_type {
        MethylationOutput::Raw => {
            let mut all_results = AHashMap::new();
            for res in per_contig_results {
                if let MethylationPatternVariant::Raw(positions) = res {
                    all_results.extend(positions.methylation);
                }
            }
            MethylationPatternVariant::Raw(MotifMethylationPositions::new(all_results))
        }
        MethylationOutput::Median => {
            let collected = per_contig_results
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
            let collected = per_contig_results
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

pub fn extract_methylation_pattern_bed<L: BatchLoader<GenomeWorkspace>>(
    loader: &mut L,
    motifs: Vec<Motif>,
    threads: usize,
    output_type: &MethylationOutput,
) -> Result<MethylationPatternVariant> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Could not initialize threadpool");

    let mut all_batch_results = Vec::new();
    let mut contigs_processed = 0;
    let mut batch_processing_time = Instant::now();

    for batch_result in
        epimetheus_io::services::data_loading_service::process_batches_from_loader(loader)
    {
        let populated_contigs = batch_result?;
        debug!("Workspace initialized");

        let batch_methylation_patterns: Result<Vec<MethylationPatternVariant>> = populated_contigs
            .par_iter()
            .map(|(_, contig)| {
                let positions = calculate_contig_read_methylation_single(contig, motifs.clone())?;

                match output_type {
                    MethylationOutput::Raw => Ok(MethylationPatternVariant::Raw(positions)),
                    MethylationOutput::Median => Ok(MethylationPatternVariant::Median(
                        positions.to_median_degrees(),
                    )),
                    MethylationOutput::WeightedMean => Ok(MethylationPatternVariant::WeightedMean(
                        positions.to_weighted_mean_degress(),
                    )),
                }
            })
            .collect();

        let batch_patterns = batch_methylation_patterns?;
        all_batch_results.extend(batch_patterns);

        contigs_processed += populated_contigs.len();
        let elapsed = batch_processing_time.elapsed();
        if contigs_processed % 100 == 0 {
            info!(
                "Finished processing {} contigs. Processing time: {}",
                contigs_processed,
                format_duration(elapsed)
            );
        }
        batch_processing_time = Instant::now();
    }

    let merged_results = match output_type {
        MethylationOutput::Raw => {
            let mut all_results = AHashMap::new();
            for res in all_batch_results {
                if let MethylationPatternVariant::Raw(positions) = res {
                    all_results.extend(positions.methylation);
                }
            }
            MethylationPatternVariant::Raw(MotifMethylationPositions::new(all_results))
        }
        MethylationOutput::Median => {
            let collected = all_batch_results
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
            let collected = all_batch_results
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
