use clap::{Args, Parser, Subcommand};
use epimetheus_core::data::contig::ContigId;

#[derive(Args, Debug)]
pub struct BgZipArgs {
    #[command(subcommand)]
    pub commands: BgZipCommands,
}

#[derive(Subcommand, Debug)]
pub enum BgZipCommands {
    Compress(BgzipWriterArgs),
    Decompress(BgzipReaderArgs),
}

#[derive(Parser, Debug, Clone)]
pub struct BgzipWriterArgs {
    #[arg(short, long, required = true, help = "Path to output pileup file.")]
    pub input: String,

    #[arg(
        short,
        long,
        required = false,
        help = "Path to output pileup file [.bed.gz]."
    )]
    pub output: Option<String>,

    #[arg(
        long,
        default_value_t = false,
        help = "Setting flag will keep the original uncompressed file."
    )]
    pub keep: bool,

    #[arg(
        short,
        long,
        required = false,
        default_value_t = 1,
        help = "Set threads for parallel compression"
    )]
    pub threads: usize,
}

#[derive(Parser, Debug, Clone)]
pub struct BgzipReaderArgs {
    #[arg(short, long, required = true, help = "Path to output pileup file.")]
    pub input: String,

    #[arg(
        short,
        long,
        required = false,
        help = "Path to output pileup file [.bed.gz]."
    )]
    pub output: Option<String>,

    #[arg(
        long,
        help = "Optional vector of contig ids to query. Left empty the whole pileup will be read."
    )]
    pub contigs: Option<Vec<ContigId>>,
}
