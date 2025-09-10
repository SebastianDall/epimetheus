use anyhow::{anyhow, Result};
// use log::{error, info, warn};
use rust_htslib::tbx::{Read, Reader};
use std::path::Path;

use crate::bgzip::args::BgzipReaderArgs;

#[derive(Debug)]
pub struct ReaderConfig {
    pub input: String,
    pub output: Option<String>,
    pub contigs: Vec<String>,
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

    pub fn query_contig(&mut self, contig: &String) -> Result<Vec<String>> {
        let tid = self.reader.tid(contig).map_err(|e| {
            anyhow!(
                "Failed to fetch contig '{}' in index: {}",
                contig,
                e.to_string()
            )
        })?;

        const UNREASONABLE_MAX: u64 = i64::MAX as u64;

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

    pub fn available_contigs(&self) -> Vec<String> {
        let contigs = self.reader.seqnames();

        contigs
    }
}

// pub fn show_header(args: &ReaderConfig) -> Result<()> {
//     let pileup_reader = PileupReader::from_path(Path::new(&args.input))?;

//     pileup_reader.list_available_contigs()?;
//     Ok(())
// }

pub fn extract_from_pileup(args: &ReaderConfig) -> Result<()> {
    todo!()
    // let mut pileup_reader = PileupReader::from_path(Path::new(&args.input))?;

    // let all_contigs = pileup_reader.available_contigs();

    // if !&args.contigs.iter().all(|c| all_contigs.contains(c)) {
    //     error!("Not all contigs are in input: {}", args.input);
    //     for contig in &args.contigs {
    //         if !all_contigs.contains(contig) {
    //             error!("- {}", contig);
    //         }
    //     }
    // }

    // Ok(())
}
