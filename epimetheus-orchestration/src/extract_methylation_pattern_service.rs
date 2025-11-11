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
        pileup::PileupRecord,
    },
    services::{domain::contig_service::populate_contig_with_methylation, traits::BatchLoader},
};
use epimetheus_io::{
    io::traits::PileupReader, loaders::sequential_batch_loader::SequentialBatchLoader,
    services::data_loading_service::load_pileup_records_for_contig,
};
use humantime::format_duration;
use indicatif::ProgressBar;
use log::{debug, info};
use methylome::Motif;
use polars::prelude::*;
use rayon::prelude::*;
use std::{collections::HashSet, io::BufReader, time::Instant};
use std::{
    fs::File,
    path::{Path, PathBuf},
};

#[derive(Debug)]
pub enum MethylationInput {
    GzFile(PathBuf),
    BedFile(PathBuf, usize),
    DataFrame(DataFrame),
}

fn merge_methylation_results(
    results: Vec<MethylationPatternVariant>,
    output_type: &MethylationOutput,
) -> MethylationPatternVariant {
    match output_type {
        MethylationOutput::Raw => {
            let mut all_results = AHashMap::new();
            for res in results {
                if let MethylationPatternVariant::Raw(positions) = res {
                    all_results.extend(positions.methylation);
                }
            }
            MethylationPatternVariant::Raw(MotifMethylationPositions::new(all_results))
        }
        MethylationOutput::Median => {
            let collected = results
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
            let collected = results
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
    }
}

pub fn extract_methylation_pattern(
    input: MethylationInput,
    contigs: AHashMap<String, Contig>,
    motifs: Vec<Motif>,
    threads: usize,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,
    output_type: &MethylationOutput,
) -> Result<MethylationPatternVariant> {
    match input {
        MethylationInput::GzFile(path) => {
            extract_methylation_patten_from_gz::<epimetheus_io::io::readers::bgzf_bed::Reader>(
                contigs,
                &path,
                motifs,
                threads,
                min_valid_read_coverage,
                min_valid_cov_to_diff_fraction,
                allow_mismatch,
                output_type,
            )
        }
        MethylationInput::BedFile(path, batch_size) => {
            let file = File::open(&path)?;
            let buf_reader = BufReader::new(file);
            let mut loader = SequentialBatchLoader::new(
                buf_reader,
                contigs,
                batch_size,
                min_valid_read_coverage,
                min_valid_cov_to_diff_fraction,
                allow_mismatch,
            );
            extract_methylation_pattern_bed(&mut loader, motifs, threads, output_type)
        }
        MethylationInput::DataFrame(df) => extract_methylation_pattern_polars(
            contigs,
            df,
            motifs,
            threads,
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
            output_type,
        ),
    }
}

fn extract_methylation_patten_from_gz<R: PileupReader + Clone>(
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

    let merged_results = merge_methylation_results(per_contig_results, output_type);

    Ok(merged_results)
}

fn extract_methylation_pattern_bed<L: BatchLoader<GenomeWorkspace>>(
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

    let merged_results = merge_methylation_results(all_batch_results, output_type);

    Ok(merged_results)
}

fn extract_methylation_pattern_polars(
    contigs: AHashMap<String, Contig>,
    pileup_df: DataFrame,
    motifs: Vec<Motif>,
    threads: usize,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
    output_type: &MethylationOutput,
) -> Result<MethylationPatternVariant> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Could not initialize threadpool");

    let pileup_records: Result<Vec<PileupRecord>, _> = (0..pileup_df.height())
        .map(|i| -> Result<PileupRecord, anyhow::Error> {
            let row = pileup_df.get_row(i)?;

            Ok(PileupRecord::new(
                row.0[0].get_str().unwrap().to_string(),
                row.0[1].try_extract::<u32>()?,
                row.0[2].try_extract::<u32>()?,
                row.0[3].get_str().unwrap().parse()?,
                row.0[4].try_extract::<u32>()?,
                row.0[5].get_str().unwrap().parse()?,
                row.0[6].try_extract::<u32>()?,
                row.0[7].try_extract::<u32>()?,
                row.0[8].get_str().unwrap().to_string(),
                row.0[9].try_extract::<u32>()?,
                row.0[10].try_extract::<f64>()?,
                row.0[11].try_extract::<u32>()?,
                row.0[12].try_extract::<u32>()?,
                row.0[13].try_extract::<u32>()?,
                row.0[14].try_extract::<u32>()?,
                row.0[15].try_extract::<u32>()?,
                row.0[16].try_extract::<u32>()?,
                row.0[17].try_extract::<u32>()?,
            ))
        })
        .collect();
    let pileup_records = pileup_records?;

    let mut meth_records = Vec::new();
    for rec in &pileup_records {
        match MethylationRecord::try_from_with_filters(
            rec.clone(),
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
        )? {
            Some(m) => meth_records.push(m),
            None => continue,
        }
    }

    let records_by_contig: AHashMap<String, Vec<MethylationRecord>> = meth_records
        .into_iter()
        .fold(AHashMap::new(), |mut acc, record| {
            acc.entry(record.get_contig_id()).or_default().push(record);
            acc
        });

    let per_contig_results = records_by_contig
        .par_iter()
        .filter_map(|(contig_id, meth_records)| {
            contigs
                .get(contig_id)
                .map(|contig| -> Result<MethylationPatternVariant> {
                    let contig_w_meth =
                        populate_contig_with_methylation(contig, meth_records.clone())?;
                    let positions =
                        calculate_contig_read_methylation_single(&contig_w_meth, motifs.clone())?;

                    match output_type {
                        MethylationOutput::Raw => Ok(MethylationPatternVariant::Raw(positions)),
                        MethylationOutput::Median => Ok(MethylationPatternVariant::Median(
                            positions.to_median_degrees(),
                        )),
                        MethylationOutput::WeightedMean => {
                            Ok(MethylationPatternVariant::WeightedMean(
                                positions.to_weighted_mean_degress(),
                            ))
                        }
                    }
                })
        })
        .collect::<Result<Vec<MethylationPatternVariant>>>()?;

    let merged_results = merge_methylation_results(per_contig_results, output_type);

    Ok(merged_results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_pileup() {
        let contig_vec = ["contig_2"];
        let start_vec = [0];
        let end_vec = [1];
        let mod_type_vec = ["a"];
        let score_vec = [20];
        let strand_vec = ["+"];
        let start_pos_vec = [0];
        let end_pos_vec = [1];
        let color_vec = ["255,0,0"];
        let n_valid_cov_vec = [20];
        let fraction_modified_vec = [1.0];
        let n_modified_vec = [20];
        let n_canonical_vec = [0];
        let n_other_mod_vec = [0];
        let n_delete_vec = [0];
        let n_fail_vec = [0];
        let n_diff_vec = [0];
        let n_no_call_vec = [0];

        let pileup_df = df!(
            "contig" => contig_vec,
            "start" => start_vec,
            "end" => end_vec,
            "mod_type" => mod_type_vec,
            "score" => score_vec,
            "strand" => strand_vec,
            "start_pos" => start_pos_vec,
            "end_pos" => end_pos_vec,
            "color" => color_vec,
            "n_valid_cov" => n_valid_cov_vec,
            "fraction_modified" => fraction_modified_vec,
            "n_modified" => n_modified_vec,
            "n_canonical" => n_canonical_vec,
            "n_other_mod" => n_other_mod_vec,
            "n_delete" => n_delete_vec,
            "n_fail" => n_fail_vec,
            "n_diff" => n_diff_vec,
            "n_no_call" => n_no_call_vec,
        )
        .unwrap();

        let pileup_records: Result<Vec<PileupRecord>, _> = (0..pileup_df.height())
            .map(|i| -> Result<PileupRecord, anyhow::Error> {
                let row = pileup_df.get_row(i)?;
                println!("{:?}", row);

                Ok(PileupRecord::new(
                    row.0[0].get_str().unwrap().to_string(),
                    row.0[1].try_extract::<u32>()?,
                    row.0[2].try_extract::<u32>()?,
                    row.0[3].get_str().unwrap().parse()?,
                    row.0[4].try_extract::<u32>()?,
                    row.0[5].get_str().unwrap().parse()?,
                    row.0[6].try_extract::<u32>()?,
                    row.0[7].try_extract::<u32>()?,
                    row.0[8].get_str().unwrap().to_string(),
                    row.0[9].try_extract::<u32>()?,
                    row.0[10].try_extract::<f64>()?,
                    row.0[11].try_extract::<u32>()?,
                    row.0[12].try_extract::<u32>()?,
                    row.0[13].try_extract::<u32>()?,
                    row.0[14].try_extract::<u32>()?,
                    row.0[15].try_extract::<u32>()?,
                    row.0[16].try_extract::<u32>()?,
                    row.0[17].try_extract::<u32>()?,
                ))
            })
            .collect();
        let pileup_records = pileup_records.unwrap();

        assert_eq!(pileup_records[0].contig, "contig_2");
    }
}
