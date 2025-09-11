use anyhow::{anyhow, Result};
use log::{debug, info};
use methylome::{ModType, Strand};
use rust_htslib::tbx::{Read, Reader};
use std::{
    fmt,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use crate::bgzip::args::BgzipExtractArgs;

pub struct PileupRecordString(pub String);

impl PileupRecordString {
    pub fn new(_0: String) -> Self {
        Self(_0)
    }
}

pub struct PileupReader {
    reader: Reader,
    records: Vec<PileupRecordString>,
    file_path: PathBuf,
}

impl Clone for PileupReader {
    fn clone(&self) -> Self {
        Self::from_path(&self.file_path).unwrap()
    }
}

impl PileupReader {
    pub fn from_path(p: &Path) -> Result<Self> {
        let reader = Reader::from_path(p)
            .map_err(|e| anyhow!("Could not open file: {:?}. Error: {}", p, e.to_string()))?;

        Ok(Self {
            reader,
            records: Vec::with_capacity(500_000),
            file_path: p.to_path_buf(),
        })
    }

    pub fn query_contig(&mut self, contig: &String) -> Result<Vec<PileupRecordString>> {
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

    pub fn available_contigs(&self) -> Vec<String> {
        let contigs = self.reader.seqnames();

        contigs
    }
}

pub fn extract_from_pileup(args: &BgzipExtractArgs) -> Result<()> {
    let mut pileup_reader = PileupReader::from_path(Path::new(&args.input))?;
    let contigs = pileup_reader.available_contigs();

    if args.ls {
        for c in contigs {
            println!("{}", c);
        }
        return Ok(());
    }

    let requested_contigs = args.resolve_contigs()?;
    info!("Writing {} contigs.", &requested_contigs.len());

    let mut writer: Box<dyn Write> = match &args.output {
        Some(out) => {
            let file = File::create(out)?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    for c in &requested_contigs {
        let records = pileup_reader.query_contig(&c)?;

        for rec in records {
            let pileup_rec = PileupRecord::try_from(rec)?;
            writeln!(writer, "{}", pileup_rec)?;
        }
    }

    Ok(())
}

struct PileupRecord {
    pub contig: String,
    pub start: u32,
    pub end: u32,
    pub mod_type: ModType,
    pub score: u32,
    pub strand: Strand,
    pub start_pos: u32,
    pub end_pos: u32,
    pub color: String,
    pub n_valid_cov: u32,
    pub fraction_modified: f64,
    pub n_modified: u32,
    pub n_canonical: u32,
    pub n_other_mod: u32,
    pub n_delete: u32,
    pub n_fail: u32,
    pub n_diff: u32,
    pub n_no_call: u32,
}

impl TryFrom<PileupRecordString> for PileupRecord {
    type Error = anyhow::Error;

    fn try_from(value: PileupRecordString) -> std::result::Result<Self, Self::Error> {
        let fields: Vec<&str> = value.0.trim().split('\t').collect();

        Ok(Self {
            contig: fields[0].to_string(),
            start: fields[1].parse()?,
            end: fields[2].parse()?,
            mod_type: fields[3].parse()?,
            score: fields[4].parse()?,
            strand: fields[5].parse()?,
            start_pos: fields[6].parse()?,
            end_pos: fields[7].parse()?,
            color: fields[8].to_string(),
            n_valid_cov: fields[9].parse()?,
            fraction_modified: fields[10].parse()?,
            n_modified: fields[11].parse()?,
            n_canonical: fields[12].parse()?,
            n_other_mod: fields[13].parse()?,
            n_delete: fields[14].parse()?,
            n_fail: fields[15].parse()?,
            n_diff: fields[16].parse()?,
            n_no_call: fields[17].parse()?,
        })
    }
}

impl fmt::Display for PileupRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.contig,
            self.start,
            self.end,
            self.mod_type.to_pileup_code(),
            self.score,
            self.strand,
            self.start_pos,
            self.end_pos,
            self.color,
            self.n_valid_cov,
            self.fraction_modified,
            self.n_modified,
            self.n_canonical,
            self.n_other_mod,
            self.n_delete,
            self.n_fail,
            self.n_diff,
            self.n_no_call,
        )
    }
}
