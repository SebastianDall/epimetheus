use anyhow::Context;
use flate2::read::GzDecoder;
use methylome::read::Read;
use noodles_fastq::{self as fastq};

use std::{fs::File, io::BufReader, path::Path};

use crate::io::traits::FastqReader;

pub struct Reader;

impl FastqReader for Reader {
    fn read_fastq(path: &Path, read_filter: Option<Vec<String>>) -> anyhow::Result<Vec<Read>> {
        let file = File::open(path)?;

        let file: Box<dyn std::io::Read> =
            if path.extension().and_then(|s| s.to_str()) == Some("gz") {
                Box::new(GzDecoder::new(file))
            } else {
                Box::new(file)
            };
        let mut reader = fastq::io::Reader::new(BufReader::new(file));

        let mut reads = Vec::new();
        let num_reads_in_filter = if let Some(f) = &read_filter {
            f.len()
        } else {
            0
        };

        let mut filtered_reads = 0;

        for result in reader.records() {
            let record: noodles_fastq::Record =
                result.with_context(|| "Error reading record from fastq file.")?;

            let id = record.name().to_string();

            if let Some(ref read_filter) = read_filter {
                if !read_filter.contains(&id) {
                    continue;
                } else {
                    filtered_reads += 1;
                }
            }

            let read = Read::from_fastq_record(record)?;

            reads.push(read);

            if num_reads_in_filter != 0 && (num_reads_in_filter == filtered_reads) {
                return Ok(reads);
            }
        }
        Ok(reads)
    }
}
