use anyhow::{anyhow, Result};
use epimetheus_core::{models::pileup::PileupRecordString, services::traits::PileupReader};
use rust_htslib::tbx::{Read, Reader as TbxReader};
use std::path::{Path, PathBuf};

pub struct Reader {
    reader: TbxReader,
    records: Vec<PileupRecordString>,
    file_path: PathBuf,
}

impl Clone for Reader {
    fn clone(&self) -> Self {
        Self::from_path(&self.file_path).unwrap()
    }
}

impl PileupReader for Reader {
    fn query_contig(
        &mut self,
        contig: &str,
    ) -> Result<Vec<epimetheus_core::models::pileup::PileupRecordString>> {
        self.records.clear();
        // let io_start = Instant::now();
        let tid = self.reader.tid(contig).map_err(|e| {
            anyhow!(
                "Failed to fetch contig '{}' in index: {}",
                contig,
                e.to_string()
            )
        })?;
        self.reader
            .fetch(tid, 0, i64::MAX as u64)
            .map_err(|e| anyhow!("Failed to fetch contig '{}': {}", contig, e.to_string()))?;
        // let io_duration = io_start.elapsed();

        // let mem_start = Instant::now();
        // let mut record_count = 0;
        for record in self.reader.records() {
            let record = record?;
            let pileup_str = PileupRecordString::new(String::from_utf8_lossy(&record).to_string());
            self.records.push(pileup_str);
            // record_count += 1;
        }
        // let mem_duration = mem_start.elapsed();
        // debug!(
        //     "Contig {}: I/O took {:?}, Processing {} records took {:?}",
        //     &contig, io_duration, record_count, mem_duration
        // );

        Ok(std::mem::take(&mut self.records))
    }

    fn available_contigs(&self) -> Vec<String> {
        self.reader.seqnames()
    }

    fn from_path(path: &Path) -> Result<Self>
    where
        Self: Sized,
    {
        let reader = TbxReader::from_path(path)
            .map_err(|e| anyhow!("Could not open file: {:?}. Error: {}", path, e.to_string()))?;

        Ok(Self {
            reader,
            records: Vec::with_capacity(500_000),
            file_path: path.to_path_buf(),
        })
    }
}
