use std::{fs::File, path::Path};

use anyhow::{Context, Result};
use epimetheus_core::models::contig::ContigId;
use methylome::{
    read::{
        BaseModifications, MethQual, MethSkipDistances, Read, convert_skip_distances_to_positions,
    },
    sequence::Sequence,
};
use noodles_bam::{self as bam};
use noodles_bgzf::{self as bgzf};
use noodles_sam as sam;

pub struct BamReader {
    reader: bam::io::IndexedReader<bgzf::io::Reader<File>>,
}

impl BamReader {
    pub fn new(bam_path: &Path) -> Result<Self> {
        let reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("Could not build bam reader.")?;

        Ok(Self { reader })
    }

    pub fn query_contigs(&mut self) -> Result<Vec<String>> {
        let header = self.reader.read_header()?;
        let reference_sequences = header.reference_sequences();

        let contigs = reference_sequences
            .iter()
            .map(|(name, _)| name.to_string())
            .collect();
        Ok(contigs)
    }

    pub fn query_contig_reads(&mut self, id: &ContigId) -> Result<Vec<Read>> {
        let header = self.reader.read_header()?;
        let region = id.parse()?;
        let query = self.reader.query(&header, &region)?;

        let mut reads = Vec::new();
        for result in query {
            let record = result?;
            let read_id = record.name().unwrap().to_string();
            let bases: Vec<u8> = record.sequence().iter().collect();
            let sequence = Sequence::from_u8(&bases).with_context(|| {
                format!(
                    "Could not parse sequence: {}",
                    String::from_utf8_lossy(&bases)
                )
            })?;

            let data = record.data();
            let mm_tags = data.get(b"MM").and_then(|value| {
                if let Ok(sam::alignment::record::data::field::Value::String(s)) = value {
                    Some(s.to_string())
                } else {
                    None
                }
            });

            let ml_tag = data.get(b"ML").and_then(|value| match value {
                Ok(sam::alignment::record::data::field::Value::Array(
                    sam::alignment::record::data::field::value::Array::UInt8(arr),
                )) => Some(arr.iter().flatten().collect::<Vec<_>>()),
                _ => None,
            });

            let meth_qualities = if let Some(ml) = ml_tag {
                ml.iter()
                    .map(|&s| MethQual::new(s))
                    .collect::<Vec<MethQual>>()
            } else {
                Vec::new()
            };

            let modifications = if let Some(mm) = mm_tags {
                let skip_distances = MethSkipDistances::from_meth_tags(mm, meth_qualities)?;
                let modifications = convert_skip_distances_to_positions(&sequence, skip_distances)?;
                modifications
            } else {
                BaseModifications::new()
            };
            let read = Read::new(read_id, sequence, modifications);

            reads.push(read);
        }

        Ok(reads)
    }
}
