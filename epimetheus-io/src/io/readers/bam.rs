use std::{fs::File, path::Path, str::FromStr};

use ahash::HashMap;
use ahash::HashMapExt;
use anyhow::{Context, Result, anyhow};
use epimetheus_core::models::contig::ContigId;
use epimetheus_methylome::{
    IupacBase, Strand,
    read::{
        BaseModifications, MethQual, MethSkipDistances, Read, ReadMapping, SkipDistance,
        convert_skip_distances_to_positions,
    },
    sequence::Sequence,
};
use noodles_bam::{self as bam, record::Data};
use noodles_bgzf::{self as bgzf};
use noodles_sam::Header;
use noodles_sam::alignment::record_buf::data::field::Value as BufValue;
use noodles_sam::alignment::record_buf::data::field::value::Array as BufArray;
use noodles_sam::{self as sam, alignment::record::cigar::Op};

pub struct BamReaderIndexed {
    reader: bam::io::IndexedReader<bgzf::io::Reader<File>>,
}

impl BamReaderIndexed {
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
                let skip_distances = MethSkipDistances::from_meth_tags(mm, meth_qualities)?;
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

fn extract_mm_tags(data: &Data) -> Option<String> {
    let mm_tags = data.get(b"MM").and_then(|value| {
        if let Ok(sam::alignment::record::data::field::Value::String(s)) = value {
            Some(s.to_string())
        } else {
            None
        }
    });

    mm_tags
}

fn extract_ml_tag(data: &Data) -> Option<Vec<u8>> {
    let ml_tag = data.get(b"ML").and_then(|value| match value {
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
            let tag_rec = TagRecord::try_from(record).ok()?;

            Some(Ok(tag_rec))
        }))
    }
}

pub struct TagRecord {
    pub read_id: String,
    pub mm_tags: Option<String>,
    pub ml_tag: Option<Vec<u8>>,
}

impl TryFrom<noodles_bam::Record> for TagRecord {
    type Error = anyhow::Error;

    fn try_from(record: noodles_bam::Record) -> std::result::Result<Self, Self::Error> {
        let read_id = record.name().ok_or(anyhow!("Missing read id"))?.to_string();

        let data = record.data();
        let mm_tags = extract_mm_tags(&data);
        let ml_tag = extract_ml_tag(&data);

        Ok(Self {
            read_id,
            mm_tags,
            ml_tag,
        })
    }
}

#[derive(Default)]
pub struct MethTags {
    pub tags: HashMap<TagKey, MethTag>,
}

impl MethTags {
    pub fn to_sam_value(&self) -> (BufValue, BufValue) {
        let mut mm_vec = Vec::new();
        let mut ml_vec = Vec::new();
        for (_key, tag) in &self.tags {
            let (mm_part, ml_part) = tag.create_tags();
            mm_vec.push(mm_part);
            ml_vec.extend(ml_part);
        }
        let mut mm_string = mm_vec.join(";");
        mm_string.push(';'); // Add a trailing ;

        let mm_value = BufValue::String(mm_string.into());
        let ml_value = BufValue::Array(BufArray::UInt8(ml_vec.into()));

        (mm_value, ml_value)
    }
    pub fn extend_tags(&mut self, other: MethTags) -> Result<()> {
        let overlapping_keys: Vec<String> = other
            .tags
            .keys()
            .filter(|k| self.tags.contains_key(k))
            .map(|k| k.to_string())
            .collect();
        if !overlapping_keys.is_empty() {
            return Err(anyhow!(
                "Following mod code overlap. Consider using --rename-tags\nKeys: {}",
                overlapping_keys.join(" ")
            ));
        }
        self.tags.extend(other.tags.into_iter());
        Ok(())
    }
    pub fn extend(&mut self, mm: String, ml: &[u8]) -> Result<()> {
        let new_meth_tags = MethTags::from_tags(mm, ml)?;
        self.extend_tags(new_meth_tags)?;
        Ok(())
    }
    pub fn from_tags(mm: String, ml: &[u8]) -> Result<Self> {
        let (tags, _) = mm.split(";").filter(|s| !s.is_empty()).try_fold(
            (HashMap::new(), 0_usize),
            |(mut map, meth_counter), m| -> Result<_> {
                let parts: Vec<&str> = m.split(",").collect();
                let key_part = parts.first().ok_or(anyhow!("Missing key part"))?;

                let tag_key = TagKey::from_str(key_part)?;

                let skip_distances: Vec<SkipDistance> = parts
                    .iter()
                    .skip(1)
                    .map(|s| s.parse::<usize>().map(SkipDistance))
                    .collect::<Result<Vec<_>, _>>()
                    .context("Failed to parse skip distances")?;

                let end_idx = meth_counter + skip_distances.len();

                if end_idx > ml.len() {
                    return Err(anyhow!("Not enough ML values"));
                }

                let meth_qualities = ml[meth_counter..end_idx]
                    .iter()
                    .map(|b| MethQual::new(*b))
                    .collect();

                let meth_tag = MethTag::new(tag_key.clone(), skip_distances, meth_qualities)?;

                map.insert(tag_key, meth_tag);

                Ok((map, end_idx))
            },
        )?;

        Ok(Self { tags })
    }

    pub fn rename_tag(&mut self, key: &TagKey, other_key: TagKey) -> Result<()> {
        if self.tags.contains_key(&other_key) {
            return Err(anyhow!(
                "Other Key already present. Could not rename key: {}",
                other_key.to_string()
            ));
        }
        let mut key_values = self
            .tags
            .remove(key)
            .ok_or_else(|| anyhow!("Key '{}' does not exist", key.to_string()))?;
        key_values.tag_key = other_key.clone();
        self.tags.insert(other_key, key_values);
        Ok(())
    }
}

pub struct MethTag {
    pub tag_key: TagKey,
    pub skip_distances: Vec<SkipDistance>,
    pub meth_qualities: Vec<MethQual>,
}

impl MethTag {
    pub fn new(
        tag_key: TagKey,
        skip_distances: Vec<SkipDistance>,
        meth_qualities: Vec<MethQual>,
    ) -> Result<Self> {
        if skip_distances.len() == meth_qualities.len() {
            Ok(Self {
                tag_key,
                skip_distances,
                meth_qualities,
            })
        } else {
            Err(anyhow!(
                "Skip distances and methylation qualities are not the same length"
            ))
        }
    }
    pub fn create_tags(&self) -> (String, Vec<u8>) {
        let mut mm_parts = vec![self.tag_key.to_string()];
        mm_parts.extend(self.skip_distances.iter().map(|s| s.0.to_string()));
        let mm = mm_parts.join(",");

        let ml = self.meth_qualities.iter().map(|q| q.0).collect();

        (mm, ml)
    }
}

/// A tag key is like so:
/// [ATGCU][+-]([a-z]+[0-9]+)[.?]?
/// example: A+a.
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct TagKey {
    fundamental_base: IupacBase,
    strand: Strand,
    base_modification_code: ModCode,
    skip_interpreter: Option<SkipInterpreter>,
}

impl TagKey {
    pub fn new(
        base: char,
        strand: char,
        base_modification_code: String,
        skip_interpreter: Option<char>,
    ) -> Result<Self> {
        let fundamental_base = IupacBase::parse_char(base)?;
        match fundamental_base {
            IupacBase::A | IupacBase::G | IupacBase::C | IupacBase::T => {
                // Valid continue
            }
            _ => {
                return Err(anyhow!(
                    "{} is not a valid fundamental base for MM tags.",
                    fundamental_base.to_string()
                ));
            }
        }

        let strand = Strand::try_from(strand)?;
        let base_modification_code = ModCode::new(base_modification_code)?;
        let skip_interpreter = skip_interpreter.and_then(|i| SkipInterpreter::try_from(i).ok());

        Ok(Self {
            fundamental_base,
            strand,
            base_modification_code,
            skip_interpreter,
        })
    }
    pub fn mutate_mod_code(&self, code: ModCode) -> TagKey {
        let mut new_key = self.clone();
        new_key.base_modification_code = code;
        new_key
    }
}

impl ToString for TagKey {
    fn to_string(&self) -> String {
        let skip_interpreter = match self.skip_interpreter {
            Some(i) => i.to_string(),
            None => "".to_string(),
        };
        format!(
            "{}{}{}{}",
            self.fundamental_base.to_string(),
            self.strand.to_string(),
            self.base_modification_code.0,
            skip_interpreter,
        )
    }
}

impl FromStr for TagKey {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s.trim().is_empty() {
            return Err(anyhow!("Cannot parse empty tag key"));
        }

        let mut chars = s.chars();

        let base = chars
            .next()
            .ok_or_else(|| anyhow!("Missing fundamental base in tag key: '{}'", s))?;

        let strand = chars
            .next()
            .ok_or_else(|| anyhow!("Missing strand in tag key: '{}'", s))?;

        let remaining: String = chars.collect();
        let (code, interpreter) = if let Some(last_char) = remaining.chars().last() {
            if last_char == '?' || last_char == '.' {
                let mod_code = &remaining[..remaining.len() - last_char.len_utf8()];
                (mod_code, Some(last_char))
            } else {
                let mod_code = &remaining[..remaining.len()];
                (mod_code, None)
            }
        } else {
            return Err(anyhow!("Missing base modification code: {}", remaining));
        };

        let tag_key = TagKey::new(base, strand, code.to_string(), interpreter)?;
        Ok(tag_key)
    }
}

/// Has to be either lowercase string ([a-z]+|[0-9]+)
#[derive(Hash, Debug, PartialEq, Eq, Clone)]
pub struct ModCode(pub String);
impl ModCode {
    pub fn new(code: String) -> Result<Self> {
        let all_lowercase = code.chars().all(|c| c.is_ascii_lowercase());
        let all_digits = code.chars().all(|c| c.is_numeric());
        if all_lowercase || all_digits {
            Ok(Self(code))
        } else {
            return Err(anyhow!(
                "Invalid mod base code. Must be lowercase string or all digits."
            ));
        }
    }
}

#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy)]
pub enum SkipInterpreter {
    Questionmark,
    Period,
}

impl FromStr for SkipInterpreter {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "?" => Ok(SkipInterpreter::Questionmark),
            "." => Ok(SkipInterpreter::Period),
            _ => return Err(anyhow!("Invalid skip interpreter: {}", s)),
        }
    }
}

impl TryFrom<char> for SkipInterpreter {
    type Error = anyhow::Error;

    fn try_from(value: char) -> std::result::Result<Self, Self::Error> {
        match value {
            '?' => Ok(SkipInterpreter::Questionmark),
            '.' => Ok(SkipInterpreter::Period),
            _ => return Err(anyhow!("Invalid skip interpreter: {}", value)),
        }
    }
}

impl ToString for SkipInterpreter {
    fn to_string(&self) -> String {
        match self {
            SkipInterpreter::Period => ".".to_string(),
            SkipInterpreter::Questionmark => "?".to_string(),
        }
    }
}
