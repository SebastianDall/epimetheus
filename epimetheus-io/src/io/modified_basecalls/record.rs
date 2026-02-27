use std::str::FromStr;

use crate::io::modified_basecalls::descriptor::ModifiedBaseDescriptor;
use crate::io::readers::bam::{extract_ml_tag, extract_mm_tags};
use ahash::{HashSet, HashSetExt};
use anyhow::{Result, anyhow};
use noodles_sam::alignment::record_buf::data::field::Value as BufValue;
use noodles_sam::alignment::record_buf::data::field::value::Array as BufArray;

#[derive(Debug)]
pub struct TagRecord {
    pub read_id: String,
    pub mm_tags: Option<String>,
    pub ml_tag: Option<Vec<u8>>,
}

impl TagRecord {
    pub fn rename_modified_base_descriptor(
        &mut self,
        from: &ModifiedBaseDescriptor,
        to: &ModifiedBaseDescriptor,
    ) {
        let from_str = from.as_str();
        let to_str = to.as_str();
        self.mm_tags = self.mm_tags.as_ref().map(|mm| {
            let new_mm_tag_vec: Vec<String> = mm
                .split(";")
                .filter(|s| !s.is_empty())
                .map(|s| {
                    let modified_bases = s.split_once(",");
                    let new_mm_tag = if let Some((k, v)) = modified_bases {
                        if k == from_str {
                            format!("{},{}", to_str, v)
                        } else {
                            s.to_string()
                        }
                    } else {
                        s.to_string()
                    };
                    new_mm_tag
                })
                .collect();
            let mut new_mm_tags = new_mm_tag_vec.join(";");
            new_mm_tags.push(';');
            new_mm_tags
        });
    }
    pub fn extend_tags_naive(&mut self, other: TagRecord) {
        if let Some(other_mm) = other.mm_tags {
            self.mm_tags = match &mut self.mm_tags {
                Some(existing) => {
                    existing.push_str(&other_mm);
                    Some(existing.clone())
                }
                None => Some(other_mm),
            }
        }
        if let Some(mut other_ml) = other.ml_tag {
            match &mut self.ml_tag {
                Some(ml) => ml.append(&mut other_ml), // No allocation!
                None => self.ml_tag = Some(other_ml),
            }
        }
    }
    /// Create the sam tags. This will consume the record
    pub fn to_sam_value(self) -> Result<Option<(BufValue, BufValue)>> {
        match (self.mm_tags, self.ml_tag) {
            (Some(mm), Some(ml)) => {
                let mm_value = BufValue::String(mm.into());
                let ml_value = BufValue::Array(BufArray::UInt8(ml.into()));
                Ok(Some((mm_value, ml_value)))
            }
            (None, Some(_)) => Err(anyhow!("Missing mm values, but present ml values")),
            (Some(_), None) => Err(anyhow!("Missing ml values, but present mm values")),
            (None, None) => Ok(None),
        }
    }

    pub fn remove_modified_base_descriptor(
        &mut self,
        remove_tags: &Vec<ModifiedBaseDescriptor>,
    ) -> Result<()> {
        if let Some(mm) = &self.mm_tags {
            let (new_mm_vec, ignored_indices, _) =
                mm.split(';').filter(|s| !s.is_empty()).try_fold(
                    (Vec::new(), HashSet::new(), 0_usize),
                    |(mut new_mm_vec, mut ignored_indices, running_counter),
                     mm_part|
                     -> Result<_> {
                        if let Some((mod_descriptor_str, vals)) = mm_part.split_once(',') {
                            let mod_descriptor =
                                ModifiedBaseDescriptor::from_str(mod_descriptor_str)?;
                            let skip_distances_len = vals.split(',').count();

                            if remove_tags.contains(&mod_descriptor) {
                                ignored_indices.extend(
                                    running_counter..(running_counter + skip_distances_len),
                                );
                            } else {
                                new_mm_vec.push(mm_part.to_string());
                            }

                            Ok((
                                new_mm_vec,
                                ignored_indices,
                                running_counter + skip_distances_len,
                            ))
                        } else {
                            Ok((new_mm_vec, ignored_indices, running_counter))
                        }
                    },
                )?;
            self.mm_tags = if new_mm_vec.is_empty() {
                None
            } else {
                Some(new_mm_vec.join(";") + ";") // Need to end with a ;
            };

            if let Some(ml) = &self.ml_tag {
                let new_ml: Vec<u8> = ml
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &val)| {
                        if ignored_indices.contains(&i) {
                            None
                        } else {
                            Some(val)
                        }
                    })
                    .collect();
                self.ml_tag = if new_ml.is_empty() {
                    None
                } else {
                    Some(new_ml)
                }
            }
        }

        Ok(())
    }
}

impl TryFrom<&noodles_bam::Record> for TagRecord {
    type Error = anyhow::Error;

    fn try_from(record: &noodles_bam::Record) -> std::result::Result<Self, Self::Error> {
        let read_id = record.name().ok_or(anyhow!("Missing read id"))?.to_string();

        let data = record.data();
        let mm_tags = extract_mm_tags(&data).map(|mm| mm.to_string());
        let ml_tag = extract_ml_tag(&data);

        Ok(Self {
            read_id,
            mm_tags,
            ml_tag,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_tag_record() -> TagRecord {
        TagRecord {
            read_id: "test_read".to_string(),
            mm_tags: Some("C+m,0,1,2;A+a,3,4;".to_string()),
            ml_tag: Some(vec![100, 150, 200, 50, 75]),
        }
    }

    #[test]
    fn test_rename_modified_base_descriptor() {
        let mut record = create_test_tag_record();
        let from = ModifiedBaseDescriptor::from_str("C+m").unwrap();
        let to = ModifiedBaseDescriptor::from_str("C+h").unwrap();

        record.rename_modified_base_descriptor(&from, &to);

        assert_eq!(record.mm_tags, Some("C+h,0,1,2;A+a,3,4;".to_string()));
        assert_eq!(record.ml_tag, Some(vec![100, 150, 200, 50, 75]));
    }

    #[test]
    fn test_rename_nonexistent_tag() {
        let mut record = create_test_tag_record();
        let from = ModifiedBaseDescriptor::from_str("T+g").unwrap();
        let to = ModifiedBaseDescriptor::from_str("T+h").unwrap();

        record.rename_modified_base_descriptor(&from, &to);

        // Should remain unchanged
        assert_eq!(record.mm_tags, Some("C+m,0,1,2;A+a,3,4;".to_string()));
    }

    #[test]
    fn test_extend_tags_naive_both_present() {
        let mut record1 = TagRecord {
            read_id: "read1".to_string(),
            mm_tags: Some("C+m,0,1;".to_string()),
            ml_tag: Some(vec![100, 150]),
        };

        let record2 = TagRecord {
            read_id: "read2".to_string(),
            mm_tags: Some("A+a,2,3;".to_string()),
            ml_tag: Some(vec![200, 250]),
        };

        record1.extend_tags_naive(record2);

        assert_eq!(record1.mm_tags, Some("C+m,0,1;A+a,2,3;".to_string()));
        assert_eq!(record1.ml_tag, Some(vec![100, 150, 200, 250]));
    }

    #[test]
    fn test_extend_tags_naive_first_empty() {
        let mut record1 = TagRecord {
            read_id: "read1".to_string(),
            mm_tags: None,
            ml_tag: None,
        };

        let record2 = TagRecord {
            read_id: "read2".to_string(),
            mm_tags: Some("A+a,2,3;".to_string()),
            ml_tag: Some(vec![200, 250]),
        };

        record1.extend_tags_naive(record2);

        assert_eq!(record1.mm_tags, Some("A+a,2,3;".to_string()));
        assert_eq!(record1.ml_tag, Some(vec![200, 250]));
    }

    #[test]
    fn test_extend_tags_naive_second_empty() {
        let mut record1 = TagRecord {
            read_id: "read1".to_string(),
            mm_tags: Some("C+m,0,1;".to_string()),
            ml_tag: Some(vec![100, 150]),
        };

        let record2 = TagRecord {
            read_id: "read2".to_string(),
            mm_tags: None,
            ml_tag: None,
        };

        record1.extend_tags_naive(record2);

        assert_eq!(record1.mm_tags, Some("C+m,0,1;".to_string()));
        assert_eq!(record1.ml_tag, Some(vec![100, 150]));
    }

    #[test]
    fn test_to_sam_value_both_present() {
        let record = create_test_tag_record();
        let result = record.to_sam_value().unwrap();

        assert!(result.is_some());
        let (mm_value, ml_value) = result.unwrap();

        if let BufValue::String(mm_str) = mm_value {
            assert_eq!(mm_str, "C+m,0,1,2;A+a,3,4;");
        } else {
            panic!("Expected String value for MM tag");
        }

        if let BufValue::Array(BufArray::UInt8(ml_arr)) = ml_value {
            let ml_vec: Vec<u8> = ml_arr.iter().copied().collect();
            assert_eq!(ml_vec, vec![100, 150, 200, 50, 75]);
        } else {
            panic!("Expected UInt8 Array for ML tag");
        }
    }

    #[test]
    fn test_to_sam_value_both_none() {
        let record = TagRecord {
            read_id: "test".to_string(),
            mm_tags: None,
            ml_tag: None,
        };
        let result = record.to_sam_value().unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_to_sam_value_mm_missing() {
        let record = TagRecord {
            read_id: "test".to_string(),
            mm_tags: None,
            ml_tag: Some(vec![100]),
        };
        let result = record.to_sam_value();
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("Missing mm values")
        );
    }

    #[test]
    fn test_to_sam_value_ml_missing() {
        let record = TagRecord {
            read_id: "test".to_string(),
            mm_tags: Some("C+m,0;".to_string()),
            ml_tag: None,
        };
        let result = record.to_sam_value();
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("Missing ml values")
        );
    }

    #[test]
    fn test_remove_modified_base_descriptor_single_tag() {
        let mut record = TagRecord {
            read_id: "test".to_string(),
            mm_tags: Some("C+m,0,1,2;".to_string()),
            ml_tag: Some(vec![100, 150, 200]),
        };

        let remove_tags = vec![ModifiedBaseDescriptor::from_str("C+m").unwrap()];
        record
            .remove_modified_base_descriptor(&remove_tags)
            .unwrap();

        assert_eq!(record.mm_tags, None);
        assert_eq!(record.ml_tag, None);
    }

    #[test]
    fn test_remove_modified_base_descriptor_one_of_two() {
        let mut record = TagRecord {
            read_id: "test".to_string(),
            mm_tags: Some("C+m,0,1,2;A+a,3,4;".to_string()),
            ml_tag: Some(vec![100, 150, 200, 50, 75]),
        };

        let remove_tags = vec![ModifiedBaseDescriptor::from_str("C+m").unwrap()];
        record
            .remove_modified_base_descriptor(&remove_tags)
            .unwrap();

        assert_eq!(record.mm_tags, Some("A+a,3,4;".to_string()));
        assert_eq!(record.ml_tag, Some(vec![50, 75]));
    }

    #[test]
    fn test_remove_modified_base_descriptor_second_of_two() {
        let mut record = TagRecord {
            read_id: "test".to_string(),
            mm_tags: Some("C+m,0,1,2;A+a,3,4;".to_string()),
            ml_tag: Some(vec![100, 150, 200, 50, 75]),
        };

        let remove_tags = vec![ModifiedBaseDescriptor::from_str("A+a").unwrap()];
        record
            .remove_modified_base_descriptor(&remove_tags)
            .unwrap();

        assert_eq!(record.mm_tags, Some("C+m,0,1,2;".to_string()));
        assert_eq!(record.ml_tag, Some(vec![100, 150, 200]));
    }

    #[test]
    fn test_remove_modified_base_descriptor_multiple_tags() {
        let mut record = TagRecord {
            read_id: "test".to_string(),
            mm_tags: Some("C+m,0,1;A+a,2,3;T+g,4,5;".to_string()),
            ml_tag: Some(vec![10, 20, 30, 40, 50, 60]),
        };

        let remove_tags = vec![
            ModifiedBaseDescriptor::from_str("C+m").unwrap(),
            ModifiedBaseDescriptor::from_str("T+g").unwrap(),
        ];
        record
            .remove_modified_base_descriptor(&remove_tags)
            .unwrap();

        assert_eq!(record.mm_tags, Some("A+a,2,3;".to_string()));
        assert_eq!(record.ml_tag, Some(vec![30, 40]));
    }

    #[test]
    fn test_remove_modified_base_descriptor_nonexistent() {
        let mut record = create_test_tag_record();
        let remove_tags = vec![ModifiedBaseDescriptor::from_str("T+g").unwrap()];

        record
            .remove_modified_base_descriptor(&remove_tags)
            .unwrap();

        // Should remain unchanged
        assert_eq!(record.mm_tags, Some("C+m,0,1,2;A+a,3,4;".to_string()));
        assert_eq!(record.ml_tag, Some(vec![100, 150, 200, 50, 75]));
    }

    #[test]
    fn test_remove_modified_base_descriptor_empty_record() {
        let mut record = TagRecord {
            read_id: "test".to_string(),
            mm_tags: None,
            ml_tag: None,
        };

        let remove_tags = vec![ModifiedBaseDescriptor::from_str("C+m").unwrap()];
        record
            .remove_modified_base_descriptor(&remove_tags)
            .unwrap();

        assert_eq!(record.mm_tags, None);
        assert_eq!(record.ml_tag, None);
    }

    #[test]
    fn test_remove_modified_base_descriptor_complex() {
        let mut record = TagRecord {
            read_id: "test".to_string(),
            mm_tags: Some("C+m,0,1,2,3,4;A+a,5;T+g,6,7,8;G+o,9,10;".to_string()),
            ml_tag: Some(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
        };

        let remove_tags = vec![
            ModifiedBaseDescriptor::from_str("A+a").unwrap(),
            ModifiedBaseDescriptor::from_str("G+o").unwrap(),
        ];
        record
            .remove_modified_base_descriptor(&remove_tags)
            .unwrap();

        // Should keep C+m (indices 0-4) and T+g (indices 6-8)
        assert_eq!(record.mm_tags, Some("C+m,0,1,2,3,4;T+g,6,7,8;".to_string()));
        assert_eq!(record.ml_tag, Some(vec![1, 2, 3, 4, 5, 7, 8, 9]));
    }
}
