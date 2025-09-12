use indicatif::ProgressBar;
// use log::{debug, error, info};
use methylome::Motif;
use rayon::prelude::*;
use std::{collections::HashSet, path::Path};

use ahash::AHashMap;
use anyhow::Result;

use crate::{
    algorithms::methylation_pattern::calculate_contig_read_methylation_single,
    models::{contig::Contig, methylation::MotifMethylationDegree},
    services::{domain::methylation_processor::process_contig, traits::PileupReader},
};

pub fn parallel_processer<R: PileupReader + Clone>(
    file: &Path,
    contigs: &AHashMap<String, Contig>,
    motifs: Vec<Motif>,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,
) -> Result<Vec<MotifMethylationDegree>> {
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

    let methylation = filtered_contigs
        .par_iter()
        .map(
            |(_contig_id, contig)| -> Result<Vec<MotifMethylationDegree>> {
                let mut reader = R::from_path(file)?;

                let contig_w_meth = process_contig(
                    &mut reader,
                    contig,
                    min_valid_read_coverage,
                    min_valid_cov_to_diff_fraction,
                )?;
                progress_bar.inc(1);
                Ok(calculate_contig_read_methylation_single(
                    &contig_w_meth,
                    motifs.clone(),
                )?)
            },
        )
        .collect::<Result<Vec<Vec<_>>>>()?
        .into_iter()
        .flatten()
        .collect();

    Ok(methylation)
}
