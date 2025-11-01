use anyhow::Context;
use epimetheus_core::services::traits::FastqReader;
use methylome::read::Read;
use noodles_fastq::{self as fastq};

use std::{fs::File, io::BufReader, path::Path};

pub struct Reader;

impl FastqReader for Reader {
    fn read_fastq(path: &Path, read_filter: Option<Vec<String>>) -> anyhow::Result<Vec<Read>> {
        let mut reader = File::open(path)
            .map(BufReader::new)
            .map(fastq::io::Reader::new)?;
        let mut reads = Vec::new();

        let num_reads_in_filter = if let Some(f) = &read_filter {
            f.len()
        } else {
            0
        };

        let mut filtered_reads = 0;

        for result in reader.records() {
            let record = result.with_context(|| "Error reading record from fastq file.")?;

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
