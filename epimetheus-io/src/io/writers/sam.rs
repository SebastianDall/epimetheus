use anyhow::Result;
use noodles_sam::alignment::Record;
use noodles_sam::alignment::io::Write as SamWrite;
use noodles_sam::{self as sam, alignment::RecordBuf};
use std::io::{BufWriter, StdoutLock, Write};
pub struct SamStdOutWriter<'a> {
    pub wtr: sam::io::writer::Writer<BufWriter<StdoutLock<'a>>>,
    header: sam::Header,
}

impl<'a> SamStdOutWriter<'a> {
    pub fn new(header: sam::header::Header) -> Result<Self> {
        let stdout = std::io::stdout().lock();
        let out = BufWriter::new(stdout);
        let mut wtr = sam::io::Writer::new(out);
        wtr.write_header(&header)?;
        wtr.get_mut().flush()?;

        Ok(Self { wtr, header })
    }

    pub fn write_record(&mut self, record_buf: RecordBuf) -> Result<()> {
        self.wtr
            .write_alignment_record(&self.header, &record_buf as &dyn Record)?;
        Ok(())
    }

    pub fn try_finish(&mut self) -> Result<()> {
        self.wtr.finish(&self.header)?;
        Ok(())
    }
}
