use clap::Parser;

#[derive(Parser, Debug, Clone)]
pub struct BgzipArgs {
    #[arg(
        short,
        long,
        required = true,
        help = "Path to output file. Must be .tsv."
    )]
    pub input: String,

    #[arg(
        short,
        long,
        required = false,
        help = "Path to output file. Must be .tsv."
    )]
    pub output: String,
}
