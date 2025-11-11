use anyhow::Result;
use log::info;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use crate::{
    io::{readers::bgzf_bed::Reader, traits::PileupReader},
    services::file_processing_service::query_pileup,
};

pub fn extract_from_pileup(
    input: &Path,
    output: Option<&Path>,
    ls: bool,
    contigs: Vec<String>,
) -> Result<()> {
    let mut reader = Reader::from_path(input)?;

    if ls {
        let contigs_available = reader.available_contigs();
        for c in contigs_available {
            println!("{}", c);
        }
        return Ok(());
    }

    let mut writer: Box<dyn Write> = match output {
        Some(out) => {
            let file = File::create(out)?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    info!("Writing {} contigs.", &contigs.len());
    for contig in contigs {
        let records = query_pileup(&mut reader, &[contig])?;
        for r in records {
            writeln!(writer, "{}", r)?;
        }
    }

    Ok(())
}
