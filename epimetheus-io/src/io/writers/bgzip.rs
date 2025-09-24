use anyhow::Result;
use epimetheus_core::models::pileup::{PileupRecord, PileupRecordString};
use noodles_bgzf::{self as bgzf};
use noodles_core::Position;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;
use std::{
    fs::File,
    io::{BufRead, BufWriter, Write},
    path::Path,
};

use crate::io::readers::bed::LineReader;

pub enum WriterType {
    File(Writer<File>),
    StdOut(Writer<BufWriter<std::io::Stdout>>),
}

impl WriterType {
    pub fn compress_from_reader<R: BufRead>(&mut self, reader: LineReader<R>) -> Result<()> {
        match self {
            WriterType::File(w) => w.compress_from_reader(reader),
            WriterType::StdOut(w) => w.compress_from_reader(reader),
        }
    }

    pub fn compress_from_lines(&mut self, lines: std::vec::IntoIter<String>) -> Result<()> {
        match self {
            WriterType::File(w) => w.compress_from_lines(lines),
            WriterType::StdOut(w) => w.compress_from_lines(lines),
        }
    }

    pub fn write_tabix(&mut self, path: &Path) -> Result<()> {
        match self {
            WriterType::File(w) => w.write_tabix(path),
            WriterType::StdOut(_) => Ok(()),
        }
    }

    pub fn finish(self) -> Result<()> {
        match self {
            WriterType::File(w) => w.finish(),
            WriterType::StdOut(w) => w.finish(),
        }
    }
}

pub struct Writer<W: Write> {
    writer: bgzf::io::Writer<W>,
    indexer: Option<tabix::index::Indexer>,
}

impl<W: Write> Writer<W> {
    pub fn write_pileup_record(&mut self, record: &PileupRecord) -> Result<()> {
        let line = format!("{}\n", record);
        let bytes = line.as_bytes();

        let start_position = self.writer.virtual_position();

        self.writer.write_all(bytes)?;
        let end_position = self.writer.virtual_position();

        if let Some(ref mut indexer) = self.indexer {
            let start_val = record.start as usize;
            let start = if start_val == 0 {
                Position::MIN
            } else {
                Position::try_from(start_val)?
            };

            let end_val = record.end as usize;
            let end = Position::try_from(end_val)?;

            let chunk = Chunk::new(start_position, end_position);

            indexer.add_record(&record.contig, start, end, chunk)?;
        }

        Ok(())
    }

    pub fn compress_from_reader<R: BufRead>(&mut self, mut reader: LineReader<R>) -> Result<()> {
        let mut line = String::new();

        while reader.read_line(&mut line)? > 0 {
            let record_string = PileupRecordString::new(line.clone());
            let record = PileupRecord::try_from(record_string)?;

            self.write_pileup_record(&record)?;
            line.clear();
        }

        Ok(())
    }

    pub fn compress_from_lines(&mut self, lines: std::vec::IntoIter<String>) -> Result<()> {
        for line in lines {
            let record_string = PileupRecordString::new(line);
            let record = PileupRecord::try_from(record_string)?;

            self.write_pileup_record(&record)?;
        }
        Ok(())
    }

    pub fn write_tabix(&mut self, path: &Path) -> Result<()> {
        assert_eq!(path.extension().unwrap(), "tbi");
        let mut tabix_writer = File::create(path).map(tabix::io::Writer::new)?;

        if let Some(indexer) = self.indexer.take() {
            let index = indexer.build();

            tabix_writer.write_index(&index)?;
        }

        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        self.writer.finish()?;
        Ok(())
    }
}

impl Writer<File> {
    pub fn from_path(output: &Path) -> Result<Self> {
        let writer = File::create(output).map(bgzf::io::Writer::new)?;
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::bed().build());

        Ok(Self {
            writer,
            indexer: Some(indexer),
        })
    }
}

impl Writer<BufWriter<std::io::Stdout>> {
    pub fn to_stdout() -> Result<Self> {
        let stdout = BufWriter::new(std::io::stdout());
        let writer = bgzf::io::Writer::new(stdout);

        Ok(Self {
            writer,
            indexer: None,
        })
    }
}
