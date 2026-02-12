use std::str::FromStr;

use crate::io::modified_basecalls::descriptor::ModifiedBaseDescriptor;
use ahash::{HashMap, HashMapExt};
use anyhow::{Context, Result, anyhow};
use epimetheus_methylome::read::{MethQual, SkipDistance};
use noodles_sam::alignment::record_buf::data::field::Value as BufValue;
use noodles_sam::alignment::record_buf::data::field::value::Array as BufArray;

#[derive(Default)]
pub struct ModifiedBasesMap {
    pub tags: HashMap<ModifiedBaseDescriptor, ModifiedBaseSkipValues>,
}

impl ModifiedBasesMap {
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
    pub fn extend_tags(&mut self, other: ModifiedBasesMap) -> Result<()> {
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
        let new_meth_tags = ModifiedBasesMap::from_tags(mm, ml)?;
        self.extend_tags(new_meth_tags)?;
        Ok(())
    }
    pub fn from_tags(mm: String, ml: &[u8]) -> Result<Self> {
        let (tags, _) = mm.split(";").filter(|s| !s.is_empty()).try_fold(
            (HashMap::new(), 0_usize),
            |(mut map, meth_counter), m| -> Result<_> {
                let parts: Vec<&str> = m.split(",").collect();
                let key_part = parts.first().ok_or(anyhow!("Missing key part"))?;

                let tag_key = ModifiedBaseDescriptor::from_str(key_part)?;

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

                let meth_tag =
                    ModifiedBaseSkipValues::new(tag_key.clone(), skip_distances, meth_qualities)?;

                map.insert(tag_key, meth_tag);

                Ok((map, end_idx))
            },
        )?;

        Ok(Self { tags })
    }

    pub fn rename_tag(
        &mut self,
        key: &ModifiedBaseDescriptor,
        other_key: ModifiedBaseDescriptor,
    ) -> Result<()> {
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

pub struct ModifiedBaseSkipValues {
    pub tag_key: ModifiedBaseDescriptor,
    pub skip_distances: Vec<SkipDistance>,
    pub meth_qualities: Vec<MethQual>,
}

impl ModifiedBaseSkipValues {
    pub fn new(
        tag_key: ModifiedBaseDescriptor,
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
