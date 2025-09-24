use anyhow::Result;
use csv::ReaderBuilder;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    ops::{Deref, DerefMut},
};

pub enum InputReader {
    File(LineReader<BufReader<File>>),
    StdIn(LineReader<BufReader<std::io::Stdin>>),
    Lines(std::vec::IntoIter<String>),
}

pub struct BedReader<R: BufRead> {
    inner: csv::Reader<R>,
}

impl<R: BufRead> Deref for BedReader<R> {
    type Target = csv::Reader<R>;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<R: BufRead> DerefMut for BedReader<R> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

impl<R: BufRead> BedReader<R> {
    pub fn new(reader: R) -> Result<Self> {
        let csv_reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(reader);
        Ok(Self { inner: csv_reader })
    }
}

pub struct LineReader<R: BufRead> {
    inner: R,
}

impl<R: BufRead> LineReader<R> {
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    pub fn read_line(&mut self, buf: &mut String) -> Result<usize> {
        self.inner.read_line(buf).map_err(Into::into)
    }
}
