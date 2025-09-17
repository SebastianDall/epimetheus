use anyhow::Result;
use epimetheus_core::services::{domain::parallel_processer::query_pileup, traits::PileupReader};
use log::info;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use crate::readers::bgzf_bed::Reader;

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

    info!("Writing {} contigs.", &contigs.len());
    let records = query_pileup(&mut reader, &contigs)?;

    let mut writer: Box<dyn Write> = match output {
        Some(out) => {
            let file = File::create(out)?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    for r in records {
        writeln!(writer, "{}", r)?;
    }

    Ok(())
}
