use std::path::PathBuf;

use ahash::HashMap;
use anyhow::{Context, Result};
use epimetheus_io::io::{
    readers::bam::{BamReader, MethTags, TagKey, TagRecord},
    writers::bam::BamWriter,
};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::data::field::Tag;
use redb::{Database, ReadableDatabase, TableDefinition};

pub mod db;

pub struct BamMergeArgs {
    pub from_bam: PathBuf,
    pub to_bam: PathBuf,
    pub db_path: PathBuf,
    pub rename_tags: HashMap<TagKey, TagKey>,
    pub ignore_tags_from_bam: Vec<String>,
    pub keep_db: bool,
    pub keep_outfile: bool,
}

pub fn bam_tag_merge_service(args: &BamMergeArgs) -> Result<()> {
    let mut db_name = args.db_path.clone();
    db_name.push("tags.redb");
    let db = Database::create(&db_name).context("Could not create DB")?;
    let table_definition: TableDefinition<&str, (&str, &[u8])> = TableDefinition::new("tags");

    let mut rdr_from_bam = BamReader::new(&args.from_bam)?;
    let write_txn = db.begin_write()?;
    {
        let mut table = write_txn.open_table(table_definition)?;
        for result in rdr_from_bam.iter_tags()? {
            let tag_record = result?;
            if let (Some(mm), Some(ml)) = (&tag_record.mm_tags, &tag_record.ml_tag) {
                table.insert(&tag_record.read_id.as_str(), (mm.as_str(), ml.as_slice()))?;
            }
        }
    }
    write_txn.commit()?;

    // TODO: remove tags,

    let mut rdr_to_bam = BamReader::new(&args.to_bam)?;

    let temp_output = args.to_bam.with_extension("tmp.bam");
    let mut writer = BamWriter::new(&temp_output, rdr_to_bam.header.clone())?;

    let read_txn = db.begin_read()?;
    let table = read_txn.open_table(table_definition)?;

    for result in rdr_to_bam.reader.records() {
        let record = result?;

        let read_id_opt = record.name();
        if read_id_opt.is_none() {
            continue;
        }

        let mut record_buf = RecordBuf::try_from_alignment_record(&rdr_to_bam.header, &record)?;

        let mut all_meth_tags = MethTags::default();
        let tag_record = TagRecord::try_from(record.clone())?;

        if let (Some(mm), Some(ml)) = (tag_record.mm_tags, tag_record.ml_tag) {
            all_meth_tags.extend(mm, &ml)?;
        }

        let read_id_bstr = read_id_opt.unwrap();
        let read_id_string = read_id_bstr.to_string();
        let maybe_from_bam_record = table.get(read_id_string.as_str())?;
        if let Some(r) = maybe_from_bam_record {
            let (mm_str, ml) = r.value();
            let mut from_bam_meth_tags = MethTags::from_tags(mm_str.to_string(), ml)?;
            if !args.rename_tags.is_empty() {
                for (key, other_key) in &args.rename_tags.clone() {
                    println!(
                        "Key: {} | Other key: {}",
                        key.to_string(),
                        other_key.to_string()
                    );
                    from_bam_meth_tags.rename_tag(key, other_key.clone())?;
                }
            }

            all_meth_tags.extend_tags(from_bam_meth_tags)?;
        }

        for key in all_meth_tags.tags.keys() {
            println!("{}", key.to_string());
        }

        let (mm_value, ml_value) = all_meth_tags.to_sam_value();
        let data = record_buf.data_mut();

        let mm_tag = Tag::BASE_MODIFICATION_PROBABILITIES;
        let ml_tag = Tag::BASE_MODIFICATIONS;

        data.remove(&mm_tag);
        data.remove(&ml_tag);

        data.insert(mm_tag, mm_value);
        data.insert(ml_tag, ml_value);

        writer.write_record(record_buf)?;
    }
    drop(rdr_to_bam);
    drop(writer);

    Ok(())
}
