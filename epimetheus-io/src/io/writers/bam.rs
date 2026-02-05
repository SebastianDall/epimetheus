use anyhow::{Result, anyhow};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::Record;
use noodles_sam::alignment::io::Write;
use noodles_sam::{self as sam, alignment::RecordBuf};
use std::{ffi::OsString, fs::File, io::BufWriter, path::PathBuf, str::FromStr};

pub struct BamWriter {
    wtr: bam::io::writer::Writer<bgzf::io::Writer<BufWriter<File>>>,
    header: sam::Header,
}

impl BamWriter {
    pub fn new(path: &PathBuf, header: sam::header::Header) -> Result<Self> {
        let ext = path.extension().ok_or(anyhow!(
            "Output file missing bam extension: {}",
            path.display()
        ))?;

        if ext != &OsString::from_str("bam").unwrap() {
            return Err(anyhow!("The output should have the bam extension"));
        }

        let file = File::create_new(path)?;
        let buf = BufWriter::new(file);

        let mut wtr = bam::io::Writer::new(buf);
        wtr.write_header(&header)?;

        Ok(Self { wtr, header })
    }

    pub fn write_record(&mut self, record_buf: RecordBuf) -> Result<()> {
        self.wtr
            .write_alignment_record(&self.header, &record_buf as &dyn Record)?;
        Ok(())
    }
}
