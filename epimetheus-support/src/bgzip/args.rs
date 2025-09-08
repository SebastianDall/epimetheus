use clap::{Args, Parser, Subcommand};

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
        default_value_t = false,
        help = "Setting flag will keep the original uncompressed file."
    )]
    pub keep: bool,
}
