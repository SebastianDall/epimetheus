use std::{fs::File, io::{BufRead, BufReader}, path::PathBuf};

use anyhow::bail;
use clap::{Args, Parser, Subcommand};
use epimetheus_io::io::readers::bed::{InputReader, LineReader};

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
    #[arg(short, long, help = "Path to output pileup file. [.bed].")]
    pub input: Option<PathBuf>,

    #[arg(long, required = false, default_value_t=false, help = "Read from stdin.")]
    pub stdin: bool,

    #[arg(
        short,
        long,
        required = false,
        help = "Path to output pileup file [.bed.gz]. If not provided the compression will outputted to stdout and not tabix will be created."
    )]
    pub output: Option<PathBuf>,

    #[arg(long, required = false, default_value_t=false, help = "Output to stdout")]
    pub stdout: bool,

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

impl BgzipWriterArgs {
    pub fn validate_input(&self) -> anyhow::Result<InputReader> {
        if self.stdin & self.keep {
            bail!("Cannot set '--keep' with '--stdin'. No file will be removed.")
        }

        let reader = match (self.input.is_some(), self.stdin) {
            (true, false) => {
                let file = File::open(&self.input.as_ref().unwrap())?;
                let rdr = LineReader::new(BufReader::new(file));
                InputReader::File(rdr)
            },
            (false, true) => InputReader::StdIn(LineReader::new(BufReader::new(std::io::stdin()))),
            (false, false) => bail!("Must specify either '--stdin' or '--input'"),
            (true, true) => bail!("Cannot specify both file '{}' and '--stdin'", self.input.as_ref().unwrap().display()),
        };

        Ok(reader)
    }

    pub fn set_output(&self) -> anyhow::Result<Option<PathBuf>> {
        self.validate_input()?;

        let output_path = match (self.stdout, &self.output) {
            (true, None) => Ok(None),
            (false, Some(output)) => Ok(Some(output.clone())),
            (false, None) => {
                match &self.input {
                    Some(input) => Ok(Some(PathBuf::from(format!("{}.gz", input.display())))),
                    None => bail!("Cannot auto-generate output filename from input, when using stdin."),
                }
            },
            (true, Some(_)) => bail!("Cannot speficy both output and stdout."),
        };

        output_path

    }

    pub fn should_remove_input_file(&self) -> bool {
        !self.keep & !self.stdin
    }
}


#[derive(Parser, Debug, Clone)]
pub struct BgzipExtractArgs {
    #[arg(short, long, required = true, help = "Path to output pileup file. [.bed.gz].")]
    pub input: PathBuf,

    #[arg(
        short,
        long,
        required = false,
        help = "Path to output pileup file [.bed]."
    )]
    pub output: Option<PathBuf>,

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
