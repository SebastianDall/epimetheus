use anyhow::anyhow;
use log::{info, warn};
use noodles_csi::binning_index::index::Header;
use noodles_tabix as tabix;
use rust_htslib::bgzf::Writer as BGZFWriter;
use std::{
    fs::File,
    io::{BufReader, Write},
    path::{Path, PathBuf},
};

use crate::bgzip::args::BgzipWriterArgs;

pub fn zip_pileup(args: &BgzipWriterArgs) -> anyhow::Result<()> {
    let input_file = File::open(&args.input)?;
    info!("Starting compression of {}", &args.input);
    if !&args.keep {
        warn!("Will remove uncompressed file after compression. Set --keep to change this.");
    }

    let mut reader = BufReader::new(input_file);

    let output_path = match &args.output {
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

    let mut writer = BGZFWriter::from_path(&output_path)?;

    std::io::copy(&mut reader, &mut writer)?;
    writer.flush()?;

    if !&args.keep {
        info!("Removing file: {}", &args.input);
        std::fs::remove_file(&args.input)?;
    }

    info!("Writing tabix.");
    write_tabix(&Path::new(&output_path))
        .map_err(|e| anyhow!("Error writing tabix: {}", e.to_string()))?;

    Ok(())
}

fn write_tabix(file: &Path) -> anyhow::Result<()> {
    let outfile = format!("{}.tbi", file.display());

    let bed_head = Header::builder()
        .set_reference_sequence_name_index(0)
        .set_start_position_index(1)
        .set_end_position_index(Some(2))
        .build();

    let index = tabix::Index::builder().set_header(bed_head).build();

    tabix::fs::write(outfile, &index)?;

    Ok(())
}
