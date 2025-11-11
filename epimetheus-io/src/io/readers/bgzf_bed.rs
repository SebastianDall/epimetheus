use anyhow::{Result, anyhow};
use epimetheus_core::models::pileup::PileupRecordString;
use noodles_bgzf::VirtualPosition;
use noodles_bgzf::io::Reader as BgzfReader;
use noodles_core::Region;
use noodles_csi::io::IndexedReader;
use noodles_csi::{BinningIndex, binning_index::Index};
use std::{
    fs::File,
    path::{Path, PathBuf},
};

use crate::io::traits::PileupReader;

pub struct Reader {
    reader: IndexedReader<BgzfReader<File>, Index<Vec<VirtualPosition>>>,
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
        let region = Region::new(contig, ..);
        let query = self
            .reader
            .query(&region)
            .map_err(|e| anyhow!("Failed to fetch contig '{}': {}", contig, e.to_string()))?;

        // .(contig).map_err(|e| {
        //     anyhow!(
        //         "Failed to fetch contig '{}' in index: {}",
        //         contig,
        //         e.to_string()
        //     )
        // })?;
        // self.reader
        //     .fetch(tid, 0, i64::MAX as u64)
        //     .map_err(|e| anyhow!("Failed to fetch contig '{}': {}", contig, e.to_string()))?;
        // let io_duration = io_start.elapsed();

        // let mem_start = Instant::now();
        // let mut record_count = 0;
        for record in query {
            let record = record?;
            let pileup_str = PileupRecordString::new(record.as_ref().to_string());
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
        let index = self
            .reader
            .index()
            .header()
            .expect("There should be tabix file [.bed.gz.tbi]");
        index
            .reference_sequence_names()
            .iter()
            .map(|seq| seq.to_string())
            .collect::<Vec<String>>()
    }

    fn from_path(path: &Path) -> Result<Self>
    where
        Self: Sized,
    {
        let reader = noodles_tabix::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| anyhow!("Could not open file: {:?}. Error: {}", path, e.to_string()))?;
        // let reader = File::open(path)
        //     .map(TbxReader::new)
        //     .map_err(|e| anyhow!("Could not open file: {:?}. Error: {}", path, e.to_string()))?;

        Ok(Self {
            reader,
            records: Vec::with_capacity(500_000),
            file_path: path.to_path_buf(),
        })
    }
}
