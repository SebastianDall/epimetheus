use indicatif::ProgressBar;
// use log::{debug, error, info};
use methylome::Motif;
use rayon::prelude::*;
use std::{collections::HashSet, path::Path};

use ahash::AHashMap;
use anyhow::Result;

use crate::{
    algorithms::methylation_pattern::calculate_contig_read_methylation_single,
    models::{
        contig::Contig,
        methylation::{MethylationOutput, MethylationPatternVariant, MotifMethylationPositions},
        pileup::PileupRecord,
    },
    services::{domain::methylation_processor::process_contig, traits::PileupReader},
};

pub fn parallel_processer<R: PileupReader + Clone>(
    file: &Path,
    contigs: &AHashMap<String, Contig>,
    motifs: Vec<Motif>,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,
    output: &MethylationOutput,
) -> Result<MethylationPatternVariant> {
    let reader = R::from_path(&file)?;
    let contigs_in_index: HashSet<String> = reader.available_contigs().into_iter().collect();

    let filtered_contigs: Vec<(&String, &Contig)> = if allow_mismatch {
        contigs
            .iter()
            .filter(|(contig_id, _)| contigs_in_index.contains(*contig_id))
            .collect()
    } else {
        contigs.iter().collect()
    };

    let progress_bar = ProgressBar::new(filtered_contigs.len() as u64);

    let per_contig_results = filtered_contigs
        .par_iter()
        .map(
            |(_contig_id, contig)| -> Result<MethylationPatternVariant> {
                let mut reader = R::from_path(file)?;

                let contig_w_meth = process_contig(
                    &mut reader,
                    contig,
                    min_valid_read_coverage,
                    min_valid_cov_to_diff_fraction,
                )?;

                let positions =
                    calculate_contig_read_methylation_single(&contig_w_meth, motifs.clone())?;

                progress_bar.inc(1);
                match output {
                    MethylationOutput::Raw => Ok(MethylationPatternVariant::Raw(positions)),
                    MethylationOutput::Median => Ok(MethylationPatternVariant::Median(
                        positions.to_median_degrees(),
                    )),
                    MethylationOutput::WeightedMean => Ok(MethylationPatternVariant::WeightedMean(
                        positions.to_weighted_mean_degress(),
                    )),
                }
            },
        )
        .collect::<Result<Vec<MethylationPatternVariant>>>()?;

    let merged_results = match output {
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

pub fn query_pileup<R: PileupReader>(
    reader: &mut R,
    contigs: &[String],
) -> Result<Vec<PileupRecord>> {
    let mut all_records = Vec::new();

    for c in contigs {
        let records = reader.query_contig(&c)?;

        for rec in records {
            let pileup_rec = PileupRecord::try_from(rec)?;
            all_records.push(pileup_rec);
        }
    }
    Ok(all_records)
}
