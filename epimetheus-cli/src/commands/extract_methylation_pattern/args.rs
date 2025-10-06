use std::path::PathBuf;

use anyhow::anyhow;
use clap::Parser;
use epimetheus_core::models::methylation::MethylationOutput;

#[derive(Parser, Debug, Clone)]
pub struct MethylationPatternArgs {
    #[arg(
        short,
        long,
        required = true,
        help = "Path to pileup. Can be .bed.gz (recommended see bgzip command) or .bed"
    )]
    pub pileup: PathBuf,

    #[arg(short, long, required = true, help = "Path to assembly.")]
    pub assembly: PathBuf,

    #[arg(long, num_args(1..), help = "Specific contigs to process. Requires that a pileup is a .bed.gz file")]
    pub contigs: Option<Vec<String>>,

    #[arg(
        short,
        long,
        required = true,
        help = "Path to output file. Must be .tsv."
    )]
    pub output: PathBuf,

    #[arg(short, long, default_value_t = 1, help = "Number of parallel tasks.")]
    pub threads: usize,

    #[arg(short, long, required = true, num_args(1..), help = "Supply chain of motifs as <motif>_<mod_type>_<mod_position>. Example: '-m GATC_a_1 RGATCY_a_2'")]
    pub motifs: Vec<String>,

    #[arg(
        long,
        default_value_t = 3,
        help = "Minimum valid read coverage for calculating methylation."
    )]
    pub min_valid_read_coverage: u32,

    #[arg(
        long,
        default_value_t = 1000,
        help = "Number of contigs to process at a time. Higher number will use more RAM."
    )]
    pub batch_size: usize,

    #[arg(
        long,
        default_value_t = 0.8,
        help = "Required fraction of valid coverage relative to different read mapping. N_valid_cov / (N_valid_cov + N_diff)"
    )]
    pub min_valid_cov_to_diff_fraction: f32,
    // #[arg(long, default_value_t = 0.9, help = "Maximum failed fraction relative to valid coverage. N_valid_cov / (N_valid_cov + N_diff)")]
    // pub : f32,
    #[arg(
        long,
        default_value_t = false,
        help = "Allow epimetheus to continue if a contig in the pileup is not present in the assembly"
    )]
    pub allow_mismatch: bool,

    #[arg(
        long,
        default_value_t = MethylationOutput::Median,
        help = "Specify the type of methylation output type. Raw will give all motif methylations for each contig."
    )]
    pub output_type: MethylationOutput,
}

impl MethylationPatternArgs {
    pub fn validate_filter(&self) -> anyhow::Result<()> {
        if let Some(_contigs) = &self.contigs {
            if self.pileup.extension().and_then(|s| s.to_str()) != Some("gz") {
                return Err(anyhow!(
                    "Pileup must be tabix compressed to use the contig filter."
                ));
            }
        }

        Ok(())
    }
}
