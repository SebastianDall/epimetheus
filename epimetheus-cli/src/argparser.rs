use clap::{Parser, Subcommand};

use crate::commands::{
    compression::args::BgZipArgs, extract_methylation_pattern::MethylationInput,
    motif_clustering::MotifClusteringArgs,
};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    MethylationPattern(MethylationInput),
    MotifCluster(MotifClusteringArgs),
    Bgzip(BgZipArgs),
}
