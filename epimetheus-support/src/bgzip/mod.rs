use std::{
    fs::File,
    io::{BufReader, Write},
};

use bgzip::{read::BGZFMultiThreadReader, BGZFWriter, Compression};

use crate::bgzip::args::BgzipArgs;

pub mod args;

pub fn zip_pileup(args: BgzipArgs) -> anyhow::Result<()> {
    let input_file = File::open(args.input)?;
    let mut reader = BufReader::new(input_file);

    let mut output_file = File::create(args.output)?;
    let mut writer = BGZFWriter::new(&mut output_file, Compression::default());

    std::io::copy(&mut reader, &mut writer)?;
    writer.close()?;

    Ok(())
}
