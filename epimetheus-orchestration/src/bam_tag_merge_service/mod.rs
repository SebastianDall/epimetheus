use std::{
    io::Write,
    path::{Path, PathBuf},
    str::FromStr,
};

use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use anyhow::{Result, anyhow};
use bstr::ByteSlice;
use epimetheus_io::io::{
    modified_basecalls::{descriptor::ModifiedBaseDescriptor, record::TagRecord},
    readers::bam::BamReader,
    writers::sam::SamStdOutWriter,
};
use indicatif::{HumanCount, ProgressBar, ProgressStyle};
use log::{error, info};
use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::{self as sam, alignment::RecordBuf};
use rayon::prelude::*;

pub struct BamMergeArgs {
    pub from_bam: PathBuf,
    pub to_bam: PathBuf,
    pub rename_tags_from_bam: Option<HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
    pub rename_tags_to_bam: Option<HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
    pub ignore_tags_from_bam: Option<Vec<ModifiedBaseDescriptor>>,
}

fn extract_modified_descriptors_from_first_record(
    bam_path: &Path,
    rename_map: Option<&HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
) -> Result<HashSet<ModifiedBaseDescriptor>> {
    let mut reader = BamReader::new(bam_path)?;

    for result in reader.iter_tags()? {
        let mut tag_record = result?;

        if let Some(_mm) = &tag_record.mm_tags {
            if let Some(mv) = rename_map {
                for (from, to) in mv {
                    tag_record.rename_modified_base_descriptor(from, to);
                }
            }
            let keys = tag_record
                .mm_tags
                .clone()
                .unwrap()
                .split(';')
                .filter(|s| !s.is_empty())
                .filter_map(|s| {
                    s.split_once(',')
                        .map(|(k, _)| ModifiedBaseDescriptor::from_str(k))
                })
                .collect::<Result<HashSet<ModifiedBaseDescriptor>, _>>()?;
            return Ok(keys);
        }
    }
    Ok(HashSet::new())
}

type TagMap = HashMap<String, (String, Vec<u8>)>;
fn extract_tags_to_map(
    from_bam: &Path,
    // rename_tags: Option<&HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
) -> Result<TagMap> {
    let mut map = HashMap::with_capacity(10_000_000);
    let mut rdr_from_bam = BamReader::new(from_bam)?;

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    pb.set_message("Extracting tags from source BAM...");
    let mut count = 0u64;
    for result in rdr_from_bam.iter_tags()? {
        let tag_record = result?;

        // if let Some(mv) = rename_tags {
        //     for (from, to) in mv {
        //         tag_record.rename_modified_base_descriptor(from, to);
        //     }
        // }

        let id = tag_record.read_id;
        let mm = tag_record.mm_tags.unwrap_or_default();
        let ml = tag_record.ml_tag.unwrap_or_default();
        map.insert(id, (mm, ml));
        count += 1;

        if count % 1000 == 0 {
            pb.set_message(format!("Extracted {} reads", HumanCount(count)));
        }
    }
    pb.finish_with_message(format!(
        "Extracted {} reads with Base Modification tags",
        HumanCount(count)
    ));

    Ok(map)
}

// fn extract_tags_to_db(
//     from_bam: &Path,
//     db: &Database,
//     table_def: TableDefinition<&str, (&str, &[u8])>,
// ) -> Result<()> {
//     const COMMIT_BATCH_SIZE: u64 = 100_000;
//
//     let mut rdr_from_bam = BamReader::new(from_bam)?;
//
//     let pb = ProgressBar::new_spinner();
//     pb.set_style(
//         ProgressStyle::default_spinner()
//             .template("{spinner:.green} [{elapsed_precise}] {msg}")
//             .unwrap(),
//     );
//     pb.set_message("Extracting tags from source BAM...");
//
//     let mut write_txn = db.begin_write()?;
//     let mut count = 0u64;
//     let mut batch_count = 0u64;
//
//     for result in rdr_from_bam.reader.records() {
//         let record = result?;
//         let read_id = record.name().context("Missing read ID")?;
//         let data = record.data();
//         let mm_opt = extract_mm_tags(&data);
//         let ml_opt = extract_ml_tag(&data);
//         if let (Some(mm), Some(ml)) = (mm_opt, ml_opt) {
//             {
//                 let mut table = write_txn.open_table(table_def)?;
//                 table.insert(
//                     read_id
//                         .to_str()
//                         .context("Invalid UTF-8 incoding in read ID")?,
//                     (
//                         mm.to_str().context("Invalid UTF-8 in MM tag")?,
//                         ml.as_slice(),
//                     ),
//                 )?;
//             }
//             count += 1;
//             batch_count += 1;
//
//             if count % 1000 == 0 {
//                 pb.set_message(format!("Extracted {} reads", HumanCount(count)));
//             }
//
//             // Commit and start new transaction every COMMIT_BATCH_SIZE records
//             if batch_count >= COMMIT_BATCH_SIZE {
//                 write_txn.commit()?;
//                 write_txn = db.begin_write()?;
//                 batch_count = 0;
//             }
//         }
//     }
//
//     // Commit any remaining records
//     if batch_count > 0 {
//         write_txn.commit()?;
//     }
//
//     pb.finish_with_message(format!(
//         "Extracted {} reads with Base Modification tags",
//         HumanCount(count)
//     ));
//
//     Ok(())
// }

fn merge_tags_for_record(
    record: &bam::Record,
    header: &sam::Header,
    map: &TagMap,
    // table: &ReadOnlyTable<&str, (&str, &[u8])>,
    from_rename_tags: Option<&HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
    to_rename_tags: Option<&HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
    remove_tags_from_bam: &Vec<ModifiedBaseDescriptor>,
) -> Result<RecordBuf> {
    let mut record_buf = RecordBuf::try_from_alignment_record(header, record)?;

    // Extract existing tags from the record
    let mut to_tag_record = TagRecord::try_from(record)?;
    // Apply renaming if specified
    if let Some(mv) = to_rename_tags {
        for (from, to) in mv {
            to_tag_record.rename_modified_base_descriptor(from, to);
        }
    }

    // Get tags from database (from the other BAM file)
    if let Some(read_id_bstr) = record.name() {
        let read_id_string = read_id_bstr.to_str()?.to_string();
        if let Some(value) = map.get(&read_id_string) {
            let (mm_str, ml_bytes) = value;
            let mut from_tag_record = TagRecord {
                read_id: read_id_string,
                mm_tags: Some(mm_str.clone()),
                ml_tag: Some(ml_bytes.clone()),
            };

            if !remove_tags_from_bam.is_empty() {
                from_tag_record.remove_modified_base_descriptor(remove_tags_from_bam)?;
            }
            if let Some(mv) = from_rename_tags {
                for (from, to) in mv {
                    from_tag_record.rename_modified_base_descriptor(from, to);
                }
            }

            to_tag_record.extend_tags_naive(from_tag_record);
        }
    }

    // Update record with merged tags
    // let (mm_value, ml_value) = all_meth_tags.to_sam_value();
    if let Some((mm_value, ml_value)) = to_tag_record.to_sam_value()? {
        let data = record_buf.data_mut();

        let mm_tag = Tag::BASE_MODIFICATIONS;
        let ml_tag = Tag::BASE_MODIFICATION_PROBABILITIES;

        data.remove(&mm_tag);
        data.remove(&ml_tag);
        data.insert(mm_tag, mm_value);
        data.insert(ml_tag, ml_value);
    }

    Ok(record_buf)
}

fn write_merged_bam(
    to_bam: &Path,
    map: TagMap,
    rename_tags_from_bam: Option<&HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
    rename_tags_to_bam: Option<&HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
    remove_tags_from_bam: &Vec<ModifiedBaseDescriptor>,
    mut writer: SamStdOutWriter,
) -> Result<()> {
    let mut rdr = BamReader::new(to_bam)?;

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    pb.set_message("Merging tags...");

    const BATCH_SIZE: usize = 1000;
    let mut batch = Vec::with_capacity(BATCH_SIZE);
    let mut count = 0u64;
    for result in rdr.reader.records() {
        let record = result?;

        // Skip records without names
        if record.name().is_none() {
            continue;
        }

        batch.push(record);

        if batch.len() >= BATCH_SIZE {
            let processed = batch
                .par_iter()
                .map(|record| {
                    merge_tags_for_record(
                        record,
                        &rdr.header,
                        &map,
                        rename_tags_from_bam,
                        rename_tags_to_bam,
                        remove_tags_from_bam,
                    )
                })
                .collect::<Result<Vec<_>>>()?;

            for record_buf in processed {
                writer.write_record(record_buf)?;
            }

            writer.wtr.get_mut().flush()?;

            count += batch.len() as u64;
            pb.set_message(format!("Processed {} reads", HumanCount(count)));
            batch.clear();
        }
    }
    if !batch.is_empty() {
        let processed = batch
            .par_iter()
            .map(|record| {
                merge_tags_for_record(
                    record,
                    &rdr.header,
                    &map,
                    rename_tags_from_bam,
                    rename_tags_to_bam,
                    remove_tags_from_bam,
                )
            })
            .collect::<Result<Vec<_>>>()?;

        for record_buf in processed {
            writer.write_record(record_buf)?;
        }
        count += batch.len() as u64;
    }

    pb.finish_with_message(format!("Merged {} reads", HumanCount(count)));
    writer.try_finish()?;
    // writer.wtr.try_finish()?;
    Ok(())
}

pub fn bam_tag_merge_service(args: &BamMergeArgs) -> Result<()> {
    let rdr_to_bam = BamReader::new(&args.to_bam)?;
    let header = rdr_to_bam.header.clone();
    drop(rdr_to_bam);

    // Step 2: Create writer immediately (this writes header to stdout)
    let writer = SamStdOutWriter::new(header.clone())?;
    use std::io::Write;
    std::io::stdout().flush()?;

    // Check for overlapping keys:
    let from_bam_ignored_keys = args.ignore_tags_from_bam.clone().unwrap_or(Vec::new());
    let mut from_bam_keys = extract_modified_descriptors_from_first_record(
        &args.from_bam,
        args.rename_tags_from_bam.as_ref(),
    )?;
    info!("Keys found in 'from_bam': {}", args.from_bam.display());
    for key in &from_bam_keys {
        if !from_bam_ignored_keys.contains(&key) {
            info!(" - {}", key.to_string());
        } else {
            info!(" - {} (Ignored)", key.to_string());
        }
    }
    from_bam_keys = from_bam_keys
        .into_iter()
        .filter(|m| !from_bam_ignored_keys.contains(m))
        .collect();
    let to_bam_keys = extract_modified_descriptors_from_first_record(
        &args.to_bam,
        args.rename_tags_to_bam.as_ref(),
    )?;
    info!("Keys found in 'to_bam': {}", args.to_bam.display());
    for key in &to_bam_keys {
        info!(" - {}", key.to_string());
    }

    let overlapping_keys: HashSet<ModifiedBaseDescriptor> = from_bam_keys
        .clone()
        .into_iter()
        .filter(|k| to_bam_keys.contains(k))
        .collect();
    if !overlapping_keys.is_empty() {
        error!("The following keys will overlap between in bam and out bam");
        for key in overlapping_keys {
            error!(" - {}", key.to_string());
        }
        return Err(anyhow!(
            "Overlap in modified descriptors found between files. Consider renaming."
        ));
    }

    let map = extract_tags_to_map(&args.from_bam)?;

    // Step 2: Read header from target BAM
    let rdr_to_bam = BamReader::new(&args.to_bam)?;
    drop(rdr_to_bam);

    // Step 3: Write merged BAM to temp file
    // let temp_output = args.to_bam.with_extension("tmp.bam");
    write_merged_bam(
        &args.to_bam,
        map,
        args.rename_tags_from_bam.as_ref(),
        args.rename_tags_to_bam.as_ref(),
        &from_bam_ignored_keys,
        writer,
    )?;

    Ok(())
}
