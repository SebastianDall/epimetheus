use std::{
    path::{Path, PathBuf},
    str::FromStr,
};

use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use anyhow::{Result, anyhow};
use epimetheus_io::io::{
    readers::bam::{BamReader, ModifiedBaseDescriptor, TagRecord},
    writers::bam::BamWriter,
};
use indicatif::{HumanCount, ProgressBar, ProgressStyle};
use log::{error, info};
use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::{self as sam, alignment::RecordBuf};
// use redb::{Database, ReadOnlyTable, ReadableDatabase, TableDefinition};
pub struct BamMergeArgs {
    pub from_bam: PathBuf,
    pub to_bam: PathBuf,
    pub db_path: PathBuf,
    pub rename_tags_from_bam: Option<HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
    pub rename_tags_to_bam: Option<HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
    pub ignore_tags_from_bam: Vec<String>,
    pub keep_db: bool,
    pub keep_outfile: bool,
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

type TagMap = HashMap<String, TagRecord>;
fn extract_tags_to_map(
    from_bam: &Path,
    rename_tags: Option<&HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
) -> Result<TagMap> {
    let mut map = HashMap::new();
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
        let mut tag_record = result?;

        if let Some(mv) = rename_tags {
            for (from, to) in mv {
                tag_record.rename_modified_base_descriptor(from, to);
            }
        }

        map.insert(tag_record.read_id.clone(), tag_record);
        count += 1;

        if count % 1000 == 0 {
            pb.set_message(format!("Extracted {} reads", HumanCount(count)));
        }
        // if let (Some(mm), Some(ml)) = (tag_record.mm_tags, tag_record.ml_tag) {
        //     map.insert(tag_record.read_id, );
        // }
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
//     let write_txn = db.begin_write()?;
//     {
//         let mut table = write_txn.open_table(table_def)?;
//         let mut count = 0u64;
//
//         for result in rdr_from_bam.iter_tags()? {
//             let tag_record = result?;
//             if let (Some(mm), Some(ml)) = (&tag_record.mm_tags, &tag_record.ml_tag) {
//                 table.insert(&tag_record.read_id.as_str(), (mm.as_str(), ml.as_slice()))?;
//                 count += 1;
//
//                 if count % 1000 == 0 {
//                     pb.set_message(format!("Extracted {} reads", HumanCount(count)));
//                 }
//             }
//         }
//         pb.finish_with_message(format!(
//             "Extracted {} reads with Base Modification tags",
//             HumanCount(count)
//         ));
//     }
//     write_txn.commit()?;
//     Ok(())
// }

fn merge_tags_for_record(
    record: &bam::Record,
    header: &sam::Header,
    map: &mut TagMap,
    // table: &ReadOnlyTable<&str, (&str, &[u8])>,
    rename_tags: Option<&HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
) -> Result<RecordBuf> {
    let mut record_buf = RecordBuf::try_from_alignment_record(header, record)?;

    // Extract existing tags from the record
    let mut tag_record = TagRecord::try_from(record)?;

    // Get tags from database (from the other BAM file)
    let read_id_opt = record.name();
    if let Some(read_id_bstr) = read_id_opt {
        let read_id_string = read_id_bstr.to_string();
        if let Some(from_bam_tag_record) = map.remove(&read_id_string) {
            // Apply renaming if specified
            if let Some(mv) = rename_tags {
                for (from, to) in mv {
                    tag_record.rename_modified_base_descriptor(from, to);
                }
            }

            tag_record.extend_tags_naive(from_bam_tag_record);
        }
    }

    // Update record with merged tags
    // let (mm_value, ml_value) = all_meth_tags.to_sam_value();
    let maybe_samtags = tag_record.to_sam_value()?;
    match maybe_samtags {
        Some((mm_value, ml_value)) => {
            let data = record_buf.data_mut();

            let mm_tag = Tag::BASE_MODIFICATION_PROBABILITIES;
            let ml_tag = Tag::BASE_MODIFICATIONS;

            data.remove(&mm_tag);
            data.remove(&ml_tag);
            data.insert(mm_tag, mm_value);
            data.insert(ml_tag, ml_value);
        }
        None => {}
    }

    Ok(record_buf)
}

fn write_merged_bam(
    to_bam: &Path,
    output_path: &Path,
    header: &sam::Header,
    map: &mut TagMap,
    // db: &Database,
    // table_def: TableDefinition<&str, (&str, &[u8])>,
    rename_tags: Option<&HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>>,
) -> Result<()> {
    let mut rdr = BamReader::new(to_bam)?;
    let mut writer = BamWriter::new(output_path, header.clone())?;

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    pb.set_message("Merging tags...");

    // let read_txn = db.begin_read()?;
    // let table = read_txn.open_table(table_def)?;

    let mut count = 0u64;
    for result in rdr.reader.records() {
        let record = result?;

        // Skip records without names
        if record.name().is_none() {
            continue;
        }

        // let record_buf = merge_tags_for_record(&record, &rdr.header, &table, rename_tags)?;
        let record_buf = merge_tags_for_record(&record, &rdr.header, map, rename_tags)?;
        writer.write_record(record_buf)?;

        count += 1;
        if count % 1000 == 0 {
            pb.set_message(format!("Processed {} reads", HumanCount(count)));
        }
    }

    pb.finish_with_message(format!("Merged {} reads", HumanCount(count)));
    Ok(())
}

pub fn bam_tag_merge_service(args: &BamMergeArgs) -> Result<()> {
    // Setup database
    // let mut db_name = args.db_path.clone();
    // db_name.push("tags.redb");
    // let db = Database::create(&db_name).context("Could not create DB")?;
    // let table_definition: TableDefinition<&str, (&str, &[u8])> = TableDefinition::new("tags");
    //
    // // Step 1: Extract tags from source BAM to database
    // extract_tags_to_db(&args.from_bam, &db, table_definition)?;
    let from_bam_keys = extract_modified_descriptors_from_first_record(
        &args.from_bam,
        args.rename_tags_from_bam.as_ref(),
    )?;
    info!("Keys found in 'from_bam': {}", args.from_bam.display());
    for key in &from_bam_keys {
        info!(" - {}", key.to_string());
    }
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

    let mut map = extract_tags_to_map(&args.from_bam, args.rename_tags_from_bam.as_ref())?;

    // Step 2: Read header from target BAM
    let rdr_to_bam = BamReader::new(&args.to_bam)?;
    let header = rdr_to_bam.header.clone();
    drop(rdr_to_bam);

    // Step 3: Write merged BAM to temp file
    let temp_output = args.to_bam.with_extension("tmp.bam");
    write_merged_bam(
        &args.to_bam,
        &temp_output,
        &header,
        // &db,
        // table_definition,
        &mut map,
        args.rename_tags_to_bam.as_ref(),
    )?;

    // TODO: Handle file cleanup and renaming based on args.keep_db and args.keep_outfile

    Ok(())
}
