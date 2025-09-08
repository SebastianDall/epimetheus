use std::path::Path;

use anyhow::{anyhow, Result};
use epimetheus_core::data::contig::ContigId;
use rust_htslib::tbx::{Read, Reader};

use crate::bgzip::args::BgzipReaderArgs;

#[derive(Debug)]
pub struct ReaderConfig {
    pub input: String,
    pub output: Option<String>,
    pub contigs: Option<Vec<ContigId>>,
}

impl From<BgzipReaderArgs> for ReaderConfig {
    fn from(value: BgzipReaderArgs) -> Self {
        Self {
            input: value.input,
            output: value.output,
            contigs: value.contigs,
        }
    }
}

pub struct PileupReader {
    reader: Reader,
}

impl PileupReader {
    pub fn from_path(p: &Path) -> Result<Self> {
        let reader = Reader::from_path(p)
            .map_err(|e| anyhow!("Could not open file: {:?}. Error: {}", p, e.to_string()))?;

        Ok(Self { reader })
    }

    pub fn query_contig(&mut self, contig: &ContigId) -> Result<Vec<String>> {
        println!("{}", contig);
        let tid = self.reader.tid(contig).map_err(|e| {
            anyhow!(
                "Failed to fetch contig '{}' in index: {}",
                contig,
                e.to_string()
            )
        })?;

        const UNREASONABLE_MAX: u64 = i32::MAX as u64;

        self.reader
            .fetch(tid, 0, UNREASONABLE_MAX)
            .map_err(|e| anyhow!("Failed to fetch contig '{}': {}", contig, e.to_string()))?;

        let mut results = Vec::new();
        for record in self.reader.records() {
            let record = record?;
            results.push(String::from_utf8_lossy(&record).to_string());
        }

        Ok(results)
    }

    pub fn list_available_contigs(&self) {
        let contigs = self.reader.seqnames();

        println!("Contigs in file");
        for c in contigs.iter() {
            println!(" {}", c);
        }
    }
}

// pub fn show_header(args: &ReaderConfig) -> Result<()> {
//     let pileup_reader = PileupReader::from_path(Path::new(&args.input))?;

//     pileup_reader.list_available_contigs()?;
//     Ok(())
// }

pub fn read_pileup(args: &ReaderConfig) -> Result<()> {
    let mut pileup_reader = PileupReader::from_path(Path::new(&args.input))?;

    match &args.contigs {
        Some(contigs) => {
            let records = pileup_reader.query_contig(&contigs[0])?;
            println!("{:#?}", records)
        }
        None => println!("No records found"),
    };

    Ok(())
}
