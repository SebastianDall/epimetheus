use anyhow::Result;
use epimetheus_core::services::{domain::parallel_processer::query_pileup, traits::PileupReader};
use log::info;
use rayon::prelude::*;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    sync::{Arc, Mutex},
};

use crate::io::readers::bgzf_bed::Reader;

pub fn extract_from_pileup(
    input: &Path,
    output: Option<&Path>,
    ls: bool,
    contigs: Vec<String>,
) -> Result<()> {
    let reader = Reader::from_path(input)?;

    if ls {
        let contigs_available = reader.available_contigs();
        for c in contigs_available {
            println!("{}", c);
        }
        return Ok(());
    }

    // let writer: Box<dyn Write> = match output {
    //     Some(out) => {
    //         let file = File::create(out)?;
    //         Box::new(BufWriter::new(file))
    //     }
    //     None => Box::new(BufWriter::new(std::io::stdout())),
    // };

    // let writer = Arc::new(Mutex::new(writer));

    info!("Writing {} contigs.", &contigs.len());
    match output {
        Some(out) => {
            let file = File::create(out)?;
            let writer = Arc::new(Mutex::new(BufWriter::new(file)));

            contigs.par_iter().try_for_each(|contig| -> Result<()> {
                let mut thread_reader = reader.clone();
                let records = query_pileup(&mut thread_reader, &[contig.to_owned()])?;

                let mut writer_guard = writer.lock().unwrap();
                for r in records {
                    writeln!(writer_guard, "{}", r)?;
                }
                Ok(())
            })?;
        }
        None => {
            let writer = Arc::new(Mutex::new(BufWriter::new(std::io::stdout())));

            contigs.par_iter().try_for_each(|contig| -> Result<()> {
                let mut thread_reader = reader.clone();
                let records = query_pileup(&mut thread_reader, &[contig.to_owned()])?;

                let mut writer_guard = writer.lock().unwrap();
                for r in records {
                    writeln!(writer_guard, "{}", r)?;
                }
                Ok(())
            })?;
        }
    }

    Ok(())
}
