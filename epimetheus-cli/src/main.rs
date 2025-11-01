use anyhow::{Result, bail};
use clap::Parser;
use epimetheus_core::services::traits::{FastaReader, FastqReader};
use epimetheus_core::services::{
    application::motif_clustering_service::motif_clustering, domain::motif_processor::create_motifs,
};

use epimetheus_io::io::readers::fastq;
use epimetheus_io::services::compression_service::CompressorService;
use epimetheus_io::services::decompression_service::extract_from_pileup;

use epimetheus_orchestration::extract_methylation_pattern_service::{
    MethylationInput, extract_methylation_pattern,
};
use epimetheus_orchestration::extract_read_methylation_service::extract_read_methylation_pattern;
use humantime::format_duration;
use indicatif::HumanDuration;
use log::{info, warn};
use polars::io::csv::write::CsvWriter;
use polars::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

mod argparser;
mod commands;
mod utils;
use argparser::Args;

pub use crate::commands::compression::args::BgZipCommands;
use crate::commands::extract_methylation_pattern::SequenceCommand;
use crate::utils::create_output_file;

fn main() -> Result<()> {
    // let guard = pprof::ProfilerGuard::new(1000).unwrap();
    let total_duration = Instant::now();
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = Args::parse();

    match args.command {
        argparser::Commands::MethylationPattern(generic_methyl_args) => {
            match &generic_methyl_args.commands {
                SequenceCommand::Contig(methyl_args) => {
                    create_output_file(&methyl_args.output)?;

                    let motifs = create_motifs(&methyl_args.motifs)?;

                    if methyl_args.contigs.is_some() {
                        methyl_args.validate_filter()?;
                    }
                    let contigs = if let Some(contigs_filter) = &methyl_args.contigs {
                        info!("Loading assembly - specified contigs provided");
                        epimetheus_io::io::readers::fasta::Reader::read_fasta(
                            &methyl_args.assembly,
                            Some(contigs_filter.clone()),
                        )?
                    } else {
                        info!("Loading assembly");
                        epimetheus_io::io::readers::fasta::Reader::read_fasta(
                            &methyl_args.assembly,
                            None,
                        )?
                    };

                    if contigs.len() == 0 {
                        bail!("No contigs found in assembly");
                    }

                    let ext = methyl_args.pileup.extension().and_then(|s| s.to_str());
                    let input = if ext == Some("gz") {
                        MethylationInput::GzFile(methyl_args.pileup.clone())
                    } else if ext == Some("bed") {
                        MethylationInput::BedFile(
                            methyl_args.pileup.clone(),
                            methyl_args.batch_size,
                        )
                    } else {
                        bail!("Unsupported file type")
                    };

                    info!("Finding methylation");
                    let meth_pattern = extract_methylation_pattern(
                        input,
                        contigs,
                        motifs,
                        methyl_args.threads,
                        methyl_args.min_valid_read_coverage,
                        methyl_args.min_valid_cov_to_diff_fraction,
                        methyl_args.allow_mismatch,
                        &methyl_args.output_type,
                    )?;

                    info!("Writing output to: {}", &methyl_args.output.display());
                    meth_pattern.write_output(&methyl_args.output)?;
                }
                SequenceCommand::Read(methyl_args) => {
                    create_output_file(&methyl_args.output)?;

                    let motifs = create_motifs(&methyl_args.motifs)?;

                    let read_ids_filter = if let Some(file) = &methyl_args.read_ids_file {
                        let mut ids = Vec::new();
                        let reader = BufReader::new(File::open(file)?);

                        for result in reader.lines() {
                            let id = result?;
                            ids.push(id);
                        }
                        info!("Found {} read ids in file.", ids.len());
                        Some(ids)
                    } else {
                        None
                    };

                    info!("Reading reads");
                    let reads = fastq::Reader::read_fastq(&methyl_args.input, read_ids_filter)?;

                    info!(
                        "{} reads loaded. Searching methylation pattern.",
                        reads.len()
                    );
                    let mut meth_pattern = extract_read_methylation_pattern(
                        reads,
                        motifs,
                        methyl_args.min_meth_quality.clone(),
                        methyl_args.threads.clone(),
                    )?;

                    info!("Writing methylation pattern");
                    let mut file = File::create(&methyl_args.output)?;
                    CsvWriter::new(&mut file)
                        .with_separator(b'\t')
                        .finish(&mut meth_pattern)?;
                    info!(
                        "Written read methylation pattern to: {}",
                        methyl_args.output.display()
                    );
                }
            }
        }
        argparser::Commands::MotifCluster(motif_cluster_args) => {
            create_output_file(&motif_cluster_args.output)?;

            motif_clustering(&motif_cluster_args.output, &motif_cluster_args.motifs)?;
        }
        argparser::Commands::Bgzip(bgzip_args) => match &bgzip_args.commands {
            BgZipCommands::Compress(compress_args) => {
                let input_reader = compress_args.validate_input()?;

                if compress_args.should_remove_input_file() {
                    warn!("'--keep' not set. This will remove the input file after compression.");
                }

                if let Some(out) = &compress_args.output {
                    if compress_args.force & out.exists() {
                        warn!("'--force' is set. This will overwrite the original file");
                    } else if !compress_args.force & out.exists() {
                        bail!(
                            "Output file '{}' already exist. Set '--force' to override.",
                            out.display()
                        );
                    }
                }

                let output = compress_args.set_output()?;
                if let Some(ref out_path) = output {
                    info!("Writing to: {}", &out_path.display());
                } else {
                    info!("Writing to stdout");
                }

                CompressorService::compress_pileup(input_reader, output.as_deref())?;

                if compress_args.should_remove_input_file() {
                    info!(
                        "Removing file: {}",
                        &compress_args.input.as_ref().unwrap().display()
                    );
                    std::fs::remove_file(&compress_args.input.as_ref().unwrap())?;
                }
            }
            BgZipCommands::Decompress(decompress_args) => {
                let contigs = decompress_args.resolve_contigs()?;
                extract_from_pileup(
                    &decompress_args.input,
                    decompress_args.output.as_deref(),
                    decompress_args.ls,
                    contigs,
                )?;
            }
        },
    }

    let elapsed_total_duration = total_duration.elapsed();
    info!(
        "Total time: {} - ({})",
        HumanDuration(elapsed_total_duration).to_string(),
        format_duration(elapsed_total_duration).to_string()
    );
    // if let Ok(report) = guard.report().build() {
    //     use std::fs::File;
    //     use std::io::Write;

    //     let mut file = File::create("flamegraph.svg").unwrap();
    //     report.flamegraph(&mut file).unwrap();
    // }
    Ok(())
}
