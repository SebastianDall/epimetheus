use bgzip::{BGZFWriter, Compression};
use log::{info, warn};
use noodles_csi::binning_index::index::Header;
use noodles_tabix as tabix;
use std::{
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use crate::bgzip::args::BgzipArgs;

pub mod args;

pub fn zip_pileup(args: BgzipArgs) -> anyhow::Result<()> {
    let input_file = File::open(&args.input)?;
    info!("Starting compression of {}", &args.input);
    if !&args.keep {
        warn!("Will remove uncompressed file after compression. Set --keep to change this.");
    }

    let mut reader = BufReader::new(input_file);

    let output_path = match args.output {
        Some(out) => {
            if !Path::new(&out).extension().map_or(false, |ext| ext == "gz") {
                anyhow::bail!("Output file must have .gz extension: {}", out);
            }
            info!("Writing to file: {}", &out);
            PathBuf::from(out)
        }
        None => {
            let mut new_out = PathBuf::from(&args.input);
            new_out.set_extension("bed.gz");
            info!("No output set. Writing to {:?}", &new_out);
            new_out
        }
    };
    let mut output_file = File::create(&output_path)?;

    let mut writer = BGZFWriter::new(&mut output_file, Compression::default());

    std::io::copy(&mut reader, &mut writer)?;
    writer.close()?;

    if !&args.keep {
        info!("Removing file: {}", &args.input);
        std::fs::remove_file(&args.input)?;
    }

    write_tabix(&Path::new(&output_path))?;

    Ok(())
}

fn write_tabix(file: &Path) -> anyhow::Result<()> {
    let outfile = format!("{}.tbi", file.display());

    let index = tabix::Index::builder()
        .set_header(Header::default())
        .build();

    tabix::fs::write(outfile, &index)?;

    Ok(())
}
