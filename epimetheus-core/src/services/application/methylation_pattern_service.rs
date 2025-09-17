use crate::{
    models::{
        genome_workspace::GenomeWorkspace,
        methylation::{MethylationOutput, MethylationPatternVariant},
    },
    services::{
        domain::{
            motif_processor::create_motifs, parallel_processer::parallel_processer,
            sequential_processer::sequential_processer,
        },
        traits::{BatchLoader, FastaReader, PileupReader},
    },
};
use anyhow::{Context, Result};
use log::{info, warn};
use std::path::Path;

pub fn extract_methylation_pattern<R, A, B>(
    pileup: &Path,
    assembly: &Path,
    threads: usize,
    motifs: &Vec<String>,
    min_valid_read_coverage: u32,
    batch_size: usize,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,
    output_type: &MethylationOutput,
) -> Result<MethylationPatternVariant>
where
    R: PileupReader + Clone,
    A: FastaReader,
    B: BatchLoader<GenomeWorkspace>,
{
    info!(
        "Running epimetheus 'methylation-pattern' with {} threads",
        threads
    );
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Could not initialize threadpool");

    let motifs = create_motifs(&motifs).context("Failed to parse motifs")?;
    info!("Successfully parsed motifs.");

    info!("Loading assembly");
    let contigs = A::read_fasta(assembly)
        .with_context(|| format!("Error loading assembly from path: '{:#?}'", assembly))?;

    if contigs.len() == 0 {
        anyhow::bail!("No contigs are loaded!");
    }
    info!("Total contigs in assembly: {}", contigs.len());

    info!("Processing Pileup");
    if allow_mismatch {
        warn!("Mismatch between contigs in pileup and assembly is allowed.");
    }

    let methylation_pattern_results = if pileup.extension().and_then(|s| s.to_str()) == Some("gz") {
        parallel_processer::<R>(
            pileup,
            &contigs,
            motifs,
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
            allow_mismatch,
            output_type,
        )?
    } else {
        let file = std::fs::File::open(pileup)?;
        let buf_reader = std::io::BufReader::new(file);
        let mut batch_loader = B::new(
            buf_reader,
            contigs,
            batch_size,
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
            allow_mismatch,
        );
        sequential_processer(&mut batch_loader, motifs, threads, output_type)?
    };

    // methylation_pattern_results.sort_meth();

    Ok(methylation_pattern_results)
}
