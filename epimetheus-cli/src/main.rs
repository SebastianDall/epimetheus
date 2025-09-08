use anyhow::Result;
use clap::Parser;
use epimetheus_support::bgzip::args::BgZipCommands;
use epimetheus_support::bgzip::reader::ReaderConfig;
use epimetheus_support::bgzip::reader::read_pileup;
use epimetheus_support::bgzip::writer::zip_pileup;
use humantime::format_duration;
use indicatif::HumanDuration;
use log::info;
use std::time::Instant;

use epimetheus_core::extract_methylation_pattern;
use epimetheus_core::motif_clustering;

mod argparser;
use argparser::Args;

fn main() -> Result<()> {
    // let guard = pprof::ProfilerGuard::new(1000).unwrap();
    let total_duration = Instant::now();
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = Args::parse();

    match &args.command {
        argparser::Commands::MethylationPattern(methyl_args) => {
            extract_methylation_pattern(methyl_args)?;
        }
        argparser::Commands::MotifCluster(motif_cluster_args) => {
            motif_clustering(motif_cluster_args)?;
        }
        argparser::Commands::Bgzip(bgzip_args) => match &bgzip_args.commands {
            BgZipCommands::Compress(compress_args) => {
                zip_pileup(compress_args)?;
            }
            BgZipCommands::Decompress(decompress_args) => {
                let reader_conf = ReaderConfig {
                    input: decompress_args.input.clone(),
                    output: decompress_args.output.clone(),
                    contigs: decompress_args.contigs.clone(),
                };
                read_pileup(&reader_conf)?;
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
