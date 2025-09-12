use anyhow::{Context, Result};
use std::{
    io::{BufWriter, Write},
    path::Path,
};

use crate::{
    algorithms::motif_processor::collapse_child_motifs,
    services::domain::motif_processor::create_motifs,
};

pub fn motif_clustering(output: &Path, motifs: &Vec<String>) -> Result<()> {
    let motifs = create_motifs(&motifs).context("Failed to parse motifs")?;
    let motifs_with_no_childs = collapse_child_motifs(&motifs);

    let outfile = std::fs::File::create(output).with_context(|| format!("{:#?}", output))?;
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
