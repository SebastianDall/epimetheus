use log::{info, warn};
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;
use std::io::BufRead;
use std::{
    fs::File,
    io::{BufReader, Write},
    path::{Path, PathBuf},
};

use crate::bgzip::args::BgzipWriterArgs;

pub fn zip_pileup(args: &BgzipWriterArgs) -> anyhow::Result<()> {
    let input_file = Path::new(&args.input);
    info!("Starting compression of {}", &args.input);
    if !&args.keep {
        warn!("Will remove uncompressed file after compression. Set --keep to change this.");
    }

    let gz_output = match &args.output {
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

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::bed().build());

    let mut writer = File::create(&gz_output).map(bgzf::io::Writer::new)?;

    let reader = File::open(&input_file)?;
    let mut buf_reader = BufReader::new(reader);
    let mut line = String::new();

    let mut start_position = writer.virtual_position();

    while buf_reader.read_line(&mut line)? > 0 {
        let fields: Vec<&str> = line.trim().split('\t').collect();

        let reference = fields[0];

        let start_val = fields[1].parse::<usize>()?;
        let start = if start_val == 0 {
            Position::MIN
        } else {
            Position::try_from(start_val)?
        };

        let end_val = fields[2].parse::<usize>()?;
        let end = Position::try_from(end_val)?;

        writer.write_all(line.as_bytes())?;

        let end_position = writer.virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        indexer.add_record(reference, start, end, chunk)?;

        start_position = end_position;
        line.clear();
    }

    writer.finish()?;

    let index = indexer.build();

    let tab_outfile = format!("{}.tbi", gz_output.display());
    let mut writer = File::create(tab_outfile).map(tabix::io::Writer::new)?;
    writer.write_index(&index)?;

    if !&args.keep {
        info!("Removing file: {}", &args.input);
        std::fs::remove_file(&args.input)?;
    }

    Ok(())
}

// fn write_tabix(file: &Path) -> anyhow::Result<()> {
// let bed_head = Header::builder()
//     .set_format(noodles_csi::binning_index::index::header::Format::Generic(
//         noodles_csi::binning_index::index::header::format::CoordinateSystem::Bed,
//     ))
//     .build();
// }
