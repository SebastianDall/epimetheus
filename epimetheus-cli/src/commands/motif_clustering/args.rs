use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug, Clone)]
pub struct MotifClusteringArgs {
    #[arg(
        short,
        long,
        required = true,
        help = "Path to output file. Must be .tsv."
    )]
    pub output: PathBuf,

    #[arg(short, long, required = true, num_args(1..), help = "Supply chain of motifs as <motif>_<mod_type>_<mod_position>. Example: '-m GATC_a_1 RGATCY_a_2'")]
    pub motifs: Vec<String>,
}
