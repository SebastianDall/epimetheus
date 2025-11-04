pub mod args;

use clap::{Args, Subcommand};

use crate::commands::extract_methylation_pattern::args::{
    ContigMethylationPatternArgs, ReadMethylationPatternArgs,
};

#[derive(Args, Debug)]
pub struct MethylationInput {
    #[command(subcommand)]
    pub commands: SequenceCommand,
}

#[derive(Subcommand, Debug)]
pub enum SequenceCommand {
    Contig(ContigMethylationPatternArgs),
    Read(ReadMethylationPatternArgs),
}
