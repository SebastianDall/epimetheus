use anyhow::Result;
use epimetheus_core::{models::pileup::PileupRecord, services::traits::PileupReader};
use log::info;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use crate::readers::bedgz::Reader;

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

    let mut writer: Box<dyn Write> = match output {
        Some(out) => {
            let file = File::create(out)?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    for c in &contigs {
        let records = reader.query_contig(&c)?;

        for rec in records {
            let pileup_rec = PileupRecord::try_from(rec)?;
            writeln!(writer, "{}", pileup_rec)?;
        }
    }

    Ok(())
}
