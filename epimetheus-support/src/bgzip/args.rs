use std::{fs::File, io::{BufRead, BufReader}, path::PathBuf};

use anyhow::bail;
use clap::{Args, Parser, Subcommand};

#[derive(Args, Debug)]
pub struct BgZipArgs {
    #[command(subcommand)]
    pub commands: BgZipCommands,
}

#[derive(Subcommand, Debug)]
pub enum BgZipCommands {
    Compress(BgzipWriterArgs),
    Decompress(BgzipExtractArgs),
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
        long,
        default_value_t = false,
        help = "Setting flag will override the file if exists."
    )]
    pub force: bool,
}

#[derive(Parser, Debug, Clone)]
pub struct BgzipExtractArgs {
    #[arg(short, long, required = true, help = "Path to output pileup file.")]
    pub input: String,

    #[arg(
        short,
        long,
        required = false,
        help = "Path to output pileup file [.bed]."
    )]
    pub output: Option<String>,

    #[arg(
        long,
        default_value_t = false,
        help = "list contig names in pileup."
    )]
    pub ls: bool,

    #[arg(
        long,
        num_args(1..), 
        required = false,
        help = "Optional vector of contig ids to query. Left empty the whole pileup will be read."
    )]
    pub contigs: Option<Vec<String>>,

    #[arg(
        long,
        required = false,
        help = "File with contig names in it."
    )]
    pub contigs_file: Option<PathBuf>
}

impl BgzipExtractArgs {
    pub fn resolve_contigs(&self) -> anyhow::Result<Vec<String>> {
        match (&self.contigs, &self.contigs_file) {
            (Some(contigs), None) => Ok(contigs.clone()),
            (None, Some(contig_file)) => {
                let file = File::open(contig_file.as_path())?;
                let reader = BufReader::new(file);

                let mut contigs = Vec::new();
                for line_result in reader.lines() {
                    let line = line_result?;
                    let trimmed = line.trim();
                    if !trimmed.is_empty() && !trimmed.starts_with('#') {
                        contigs.push(trimmed.to_string());
                    }
                }

                if contigs.is_empty() {
                    bail!("No contigs found in file");
                }
                Ok(contigs)
            }
            (Some(_), Some(_)) => bail!("Cannot specify both --contigs and --contigs-file"),
            (None, None) => {
                if self.ls {
                    Ok(Vec::new())
                } else {
                    bail!("Must specify either --contigs or --contigs-file")
                }
            }
        }
    }
}
