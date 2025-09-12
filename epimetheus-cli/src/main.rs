use anyhow::Result;
use clap::Parser;
use epimetheus_core::services::application::{
    methylation_pattern_service::extract_methylation_pattern,
    motif_clustering_service::motif_clustering,
};

use epimetheus_io::compression::bgzip::compressor::zip_pileup;
use epimetheus_io::compression::bgzip::decompressor::extract_from_pileup;
use epimetheus_io::loaders::sequential_batch_loader::SequentialBatchLoader;
use epimetheus_io::readers::bedgz::Reader as GzPileupReader;
use epimetheus_io::readers::fasta::Reader as FastaReader;

use humantime::format_duration;
use indicatif::HumanDuration;
use log::info;
use std::time::Instant;

mod argparser;
mod commands;
mod utils;
use argparser::Args;

pub use crate::commands::compression::args::BgZipCommands;
use crate::utils::create_output_file;

fn main() -> Result<()> {
    // let guard = pprof::ProfilerGuard::new(1000).unwrap();
    let total_duration = Instant::now();
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = Args::parse();

    match &args.command {
        argparser::Commands::MethylationPattern(methyl_args) => {
            create_output_file(&methyl_args.output)?;

            extract_methylation_pattern::<
                GzPileupReader,
                FastaReader,
                SequentialBatchLoader<std::io::BufReader<std::fs::File>>,
            >(
                &methyl_args.pileup,
                &methyl_args.assembly,
                &methyl_args.output,
                methyl_args.threads,
                &methyl_args.motifs,
                methyl_args.min_valid_read_coverage,
                methyl_args.batch_size,
                methyl_args.min_valid_cov_to_diff_fraction,
                methyl_args.allow_mismatch,
            )?;
        }
        argparser::Commands::MotifCluster(motif_cluster_args) => {
            create_output_file(&motif_cluster_args.output)?;

            motif_clustering(&motif_cluster_args.output, &motif_cluster_args.motifs)?;
        }
        argparser::Commands::Bgzip(bgzip_args) => match &bgzip_args.commands {
            BgZipCommands::Compress(compress_args) => {
                zip_pileup(
                    &compress_args.input,
                    compress_args.output.as_deref(),
                    compress_args.keep,
                    compress_args.force,
                )?;
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
