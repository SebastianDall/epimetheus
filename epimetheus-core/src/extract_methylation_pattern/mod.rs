use anyhow::{Context, Result};
use log::{info, warn};
use std::{
    io::{BufWriter, Write},
    path::Path,
};

use crate::{
    data_load::load_contigs,
    extract_methylation_pattern::reader::{parallel_processer, sequential_processer},
    processing::create_motifs,
    utils::create_output_file,
};

pub mod args;
pub mod batch_loader;
pub mod parallel_batch_loader;
mod reader;
pub mod utils;

pub use args::MethylationPatternArgs;
pub use utils::parse_to_methylation_record;

pub fn extract_methylation_pattern(args: &MethylationPatternArgs) -> Result<()> {
    info!(
        "Running epimetheus 'methylation-pattern' with {} threads",
        &args.threads
    );
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .expect("Could not initialize threadpool");

    let outpath = Path::new(&args.output);
    create_output_file(&outpath)?;

    let motifs = match &args.motifs {
        Some(motifs) => {
            info!("Motifs loaded");
            motifs
        }
        _ => {
            anyhow::bail!("No motifs found");
        }
    };

    let motifs = create_motifs(motifs.clone()).context("Failed to parse motifs")?;
    info!("Successfully parsed motifs.");

    info!("Loading assembly");
    let contigs = load_contigs(&args.assembly)
        .with_context(|| format!("Error loading assembly from path: '{}'", args.assembly))?;

    if contigs.len() == 0 {
        anyhow::bail!("No contigs are loaded!");
    }
    info!("Total contigs in assembly: {}", contigs.len());

    info!("Processing Pileup");
    if args.allow_assembly_pileup_mismatch {
        warn!("Mismatch between contigs in pileup and assembly is allowed.");
    }
    let file = Path::new(&args.pileup);

    let mut methylation_pattern_results = if file.extension().and_then(|s| s.to_str()) == Some("gz")
    {
        parallel_processer(
            file,
            &contigs,
            motifs,
            args.min_valid_read_coverage,
            args.min_valid_cov_to_diff_fraction,
            args.allow_assembly_pileup_mismatch,
        )?
    } else {
        sequential_processer(
            file,
            contigs,
            motifs,
            args.min_valid_read_coverage,
            args.min_valid_cov_to_diff_fraction,
            args.allow_assembly_pileup_mismatch,
            args.threads,
        )?
    };

    methylation_pattern_results.sort_by(|a, b| a.contig.cmp(&b.contig));

    let outfile = std::fs::File::create(outpath)
        .with_context(|| format!("Failed to create file at: {:?}", outpath))?;
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
