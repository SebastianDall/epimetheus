use anyhow::{Result, anyhow};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_bgzf::io::MultithreadedWriter;
use noodles_sam::alignment::Record;
use noodles_sam::alignment::io::Write as SamWrite;
use noodles_sam::{self as sam, alignment::RecordBuf};
use std::io::{StdoutLock, Write as StdWrite};
use std::num::NonZero;
use std::path::Path;
use std::{ffi::OsString, fs::File, str::FromStr};

pub struct BamWriter {
    wtr: bam::io::writer::Writer<MultithreadedWriter<File>>,
    header: sam::Header,
}

impl BamWriter {
    pub fn new(path: &Path, header: sam::header::Header) -> Result<Self> {
        let ext = path.extension().ok_or(anyhow!(
            "Output file missing bam extension: {}",
            path.display()
        ))?;

        if ext != &OsString::from_str("bam").unwrap() {
            return Err(anyhow!("The output should have the bam extension"));
        }

        let file = File::create_new(path)?;

        let bgzf_wtr =
            bgzf::io::MultithreadedWriter::with_worker_count(NonZero::new(6).unwrap(), file);
        let mut wtr = bam::io::Writer::from(bgzf_wtr);
        wtr.write_header(&header)?;

        Ok(Self { wtr, header })
    }

    pub fn write_record(&mut self, record_buf: RecordBuf) -> Result<()> {
        self.wtr
            .write_alignment_record(&self.header, &record_buf as &dyn Record)?;
        Ok(())
    }
}

pub struct BamStdOutWriter<'a> {
    pub wtr: bam::io::writer::Writer<bgzf::io::Writer<StdoutLock<'a>>>,
    header: sam::Header,
}

impl<'a> BamStdOutWriter<'a> {
    pub fn new(header: sam::header::Header) -> Result<Self> {
        let out = std::io::stdout().lock();
        let mut wtr = bam::io::Writer::new(out);
        wtr.write_header(&header)?;
        wtr.get_mut().flush()?;

        Ok(Self { wtr, header })
    }

    pub fn write_record(&mut self, record_buf: RecordBuf) -> Result<()> {
        self.wtr
            .write_alignment_record(&self.header, &record_buf as &dyn Record)?;
        Ok(())
    }
}
