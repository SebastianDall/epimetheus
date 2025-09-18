use anyhow::Result;
use std::{fs::File, io::BufRead, path::Path};

pub struct Reader<R: BufRead> {
    reader: csv::Reader<R>,
}

impl Reader<File> {
    pub fn from_path(path: &Path) -> Result<Self> {
        let writer = File::open(path).map(csv::ReaderBuilder::new())?.
    }
}
