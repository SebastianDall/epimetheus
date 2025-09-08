use anyhow::Result;
use clap::Parser;
use epimetheus_support::bgzip::zip_pileup;
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
            let methyl_args = methyl_args.clone();
            extract_methylation_pattern(methyl_args)?;
        }
        argparser::Commands::MotifCluster(motif_cluster_args) => {
            let motif_cluster_args = motif_cluster_args.clone();
            motif_clustering(motif_cluster_args)?;
        }
        argparser::Commands::Bgzip(bgzip_args) => {
            let bgzip_args = bgzip_args.clone();
            zip_pileup(bgzip_args)?;
        }
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
