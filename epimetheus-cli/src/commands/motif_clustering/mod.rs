use anyhow::{Context, Result};
use epimetheus_core::algorithms::motif_processor::{collapse_child_motifs, create_motifs};
use log::info;
use std::{
    io::{BufWriter, Write},
    path::Path,
};

pub mod args;

pub use args::MotifClusteringArgs;

use crate::utils::create_output_file;

pub fn motif_clustering(args: &MotifClusteringArgs) -> Result<()> {
    let outpath = Path::new(&args.output);

    create_output_file(outpath)?;

    let motifs = match &args.motifs {
        Some(motifs) => {
            info!("Motifs loaded");
            create_motifs(motifs.clone()).context("Failed to parse motifs")?
        }
        _ => {
            anyhow::bail!("No motifs found");
        }
    };

    let motifs_with_no_childs = collapse_child_motifs(&motifs);

    let outfile = std::fs::File::create(outpath)
        .with_context(|| format!("Failed to create file at: {:?}", outpath))?;
    let mut writer = BufWriter::new(outfile);

    writeln!(writer, "motif\tmod_type\tmod_position")?;
    for m in motifs_with_no_childs {
        writeln!(
            writer,
            "{}\t{}\t{}",
            m.sequence_to_string(),
            m.mod_type.to_pileup_code(),
            m.mod_position
        )?;
    }

    Ok(())
}
