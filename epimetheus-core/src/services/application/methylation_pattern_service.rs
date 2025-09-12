use crate::{
    models::genome_workspace::GenomeWorkspace,
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
use std::{
    io::{BufWriter, Write},
    path::Path,
};

pub fn extract_methylation_pattern<R, A, B>(
    pileup: &Path,
    assembly: &Path,
    output: &Path,
    threads: usize,
    motifs: &Vec<String>,
    min_valid_read_coverage: u32,
    batch_size: usize,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,
) -> Result<()>
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

    let mut methylation_pattern_results =
        if pileup.extension().and_then(|s| s.to_str()) == Some("gz") {
            parallel_processer::<R>(
                pileup,
                &contigs,
                motifs,
                min_valid_read_coverage,
                min_valid_cov_to_diff_fraction,
                allow_mismatch,
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
            sequential_processer(&mut batch_loader, motifs, threads)?
        };

    methylation_pattern_results.sort_by(|a, b| a.contig.cmp(&b.contig));

    let outfile = std::fs::File::create(output)
        .with_context(|| format!("Failed to create file at: {:?}", output))?;
    let mut writer = BufWriter::new(outfile);

    writeln!(
        writer,
        "contig\tmotif\tmod_type\tmod_position\tmedian\tmean_read_cov\tN_motif_obs\tmotif_occurences_total"
    )?;

    for entry in &methylation_pattern_results {
        let motif_sequence = entry.motif.sequence_to_string();
        let mod_type_str = entry.motif.mod_type.to_pileup_code();
        let mod_position = entry.motif.mod_position;

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            entry.contig,
            motif_sequence,
            mod_type_str,
            mod_position,
            entry.median,
            entry.mean_read_cov,
            entry.n_motif_obs,
            entry.motif_occurences_total
        )?;

        writer.flush()?;
    }

    Ok(())
}
