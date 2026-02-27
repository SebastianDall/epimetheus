use anyhow::{Context, Result, anyhow};
use bstr::{BStr, ByteSlice};
use epimetheus_core::models::contig::ContigId;
use epimetheus_methylome::{
    Strand,
    read::{
        BaseModifications, MethQual, MethSkipDistances, Read, ReadMapping,
        convert_skip_distances_to_positions,
    },
    sequence::Sequence,
};
use noodles_bam::{self as bam, record::Data};
use noodles_bgzf::{self as bgzf};
use noodles_sam::Header;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::{self as sam, alignment::record::cigar::Op};
use std::{fs::File, path::Path};

use crate::io::modified_basecalls::record::TagRecord;

pub struct BamReaderIndexed {
    reader: bam::io::IndexedReader<bgzf::io::Reader<File>>,
}

impl BamReaderIndexed {
    pub fn new(bam_path: &Path) -> Result<Self> {
        let reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .context("Could not build bam reader. Did you remember to create the index file?")?;

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
        for result in query.records() {
            let record = result?;
            let flags = record.flags();

            if flags.is_secondary() {
                continue;
            }

            let read_id = record.name().unwrap().to_string();

            let strand = if flags.is_reverse_complemented() {
                Strand::Negative
            } else {
                Strand::Positive
            };
            let bases: Vec<u8> = record.sequence().iter().collect();
            let mut sequence = Sequence::from_u8(&bases).with_context(|| {
                format!(
                    "Could not parse sequence: {}",
                    String::from_utf8_lossy(&bases)
                )
            })?;

            sequence = match strand {
                Strand::Positive => sequence,
                Strand::Negative => sequence.reverse_complement(),
            };

            let alignment_start = if let Some(pos) = record.alignment_start() {
                if let Ok(pos_ok) = pos {
                    pos_ok.get() - 1
                } else {
                    0
                }
            } else {
                return Err(anyhow!("{} not mapped to contig: {}", read_id, id));
            };

            let cigar = record.cigar();
            let cigar_ops: Vec<Op> = cigar.iter().filter_map(|o| o.ok()).collect();
            let mapping_quality = record.mapping_quality().map(|mq| mq.get()).unwrap_or(0);

            let mapping = Some(ReadMapping::new(
                id.clone(),
                alignment_start,
                strand,
                cigar_ops,
                mapping_quality,
            ));

            let data = record.data();
            let mm_tags = extract_mm_tags(&data);
            let ml_tag = extract_ml_tag(&data);

            let meth_qualities = if let Some(ml) = ml_tag {
                ml.iter()
                    .map(|&s| MethQual::new(s))
                    .collect::<Vec<MethQual>>()
            } else {
                Vec::new()
            };

            let modifications = if let Some(mm) = mm_tags {
                let skip_distances =
                    MethSkipDistances::from_meth_tags(mm.to_str()?, meth_qualities)?;
                let modifications = convert_skip_distances_to_positions(&sequence, skip_distances)?;
                modifications
            } else {
                BaseModifications::new()
            };
            let read = Read::new_with_mapping(read_id, sequence, modifications, mapping);

            reads.push(read);
        }

        Ok(reads)
    }
}

pub fn extract_mm_tags<'a>(data: &'a Data) -> Option<&'a BStr> {
    let mm_tags = data.get(&Tag::BASE_MODIFICATIONS).and_then(|value| {
        if let Ok(sam::alignment::record::data::field::Value::String(s)) = value {
            // Some(s.to_string())
            Some(s)
        } else {
            None
        }
    });

    mm_tags
}

pub fn extract_ml_tag(data: &Data) -> Option<Vec<u8>> {
    let ml_tag = data
        .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
        .and_then(|value| match value {
            Ok(sam::alignment::record::data::field::Value::Array(
                sam::alignment::record::data::field::value::Array::UInt8(arr),
            )) => Some(arr.iter().flatten().collect::<Vec<_>>()),
            _ => None,
        });

    ml_tag
}

pub struct BamReader {
    pub reader: bam::io::Reader<bgzf::io::Reader<File>>,
    pub header: Header,
}

impl BamReader {
    pub fn new(path: &Path) -> Result<Self> {
        let mut reader = File::open(path).map(bam::io::Reader::new)?;
        let header = reader.read_header()?;

        Ok(Self { reader, header })
    }

    pub fn iter_tags(&mut self) -> Result<impl Iterator<Item = Result<TagRecord>> + '_> {
        Ok(self.reader.records().filter_map(move |result| {
            let record = result.ok()?;
            let tag_rec = TagRecord::try_from(&record).ok()?;

            Some(Ok(tag_rec))
        }))
    }
}
