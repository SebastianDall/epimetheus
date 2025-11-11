use anyhow::{Context, Result};
use epimetheus_io::io::readers::bam::BamReader;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use methylome::{Motif, find_motif_indices_in_sequence};
use rayon::prelude::*;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    sync::mpsc,
    thread,
};

pub fn extract_read_methylation_pattern(
    input_file: &Path,
    contigs_filter: Option<Vec<String>>,
    motifs: Vec<Motif>,
    output: &Path,
    threads: usize,
) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Could not initialize threadpool");

    let mut reader = BamReader::new(input_file)?;

    let contigs: Vec<String> = reader
        .query_contigs()?
        .into_iter()
        .filter(|c| {
            contigs_filter
                .as_ref()
                .map_or(true, |filter| filter.contains(c))
        })
        .collect();

    // multiprogressbar
    let multi = MultiProgress::new();
    let main_pb = multi.add(ProgressBar::new(contigs.len() as u64));
    main_pb.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} contigs")?
            .progress_chars("#>-"),
    );
    main_pb.set_message("Processing contigs");

    let writes_pb = multi.add(ProgressBar::new_spinner());
    writes_pb.set_style(
        ProgressStyle::default_spinner()
            .template("📝 {msg} [{elapsed_precise}] {spinner:.green} {pos} lines written")?,
    );
    writes_pb.set_message("Writing");

    let writes_pb_clone = writes_pb.clone();

    let (sender, receiver) = mpsc::channel();

    let output_path = output.to_path_buf();
    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut writer = BufWriter::new(File::create(&output_path)?);
        writeln!(
            writer,
            "contig_id\tread_id\tread_length\tmotif\tmod_type\tmod_position\tquality"
        )?;

        let mut lines_written = 0;
        while let Ok(line) = receiver.recv() {
            writeln!(writer, "{}", line)?;
            lines_written += 1;

            if lines_written % 1000 == 0 {
                writes_pb_clone.inc(1000);
                lines_written = 0;
            }
        }
        writer.flush()?;
        Ok(())
    });

    contigs.par_iter().try_for_each(|contig_id| -> Result<()> {
        main_pb.inc(1);
        let mut local_reader = BamReader::new(input_file)?;
        let reads = local_reader
            .query_contig_reads(contig_id)
            .with_context(|| format!("Reading contig: {}", contig_id))?;

        if reads.is_empty() {
            return Ok(());
        }

        for read in reads {
            let sequence = read.get_sequence();
            let modifications = read.get_modifications();
            // let read_length = read.get_sequence().len();

            for motif in &motifs {
                let indices = find_motif_indices_in_sequence(sequence, motif);
                // let mod_type = motif.mod_type;

                for &pos in &indices {
                    let quality = if let Some(meth_base) = modifications.0.get(&pos) {
                        meth_base.quality.0
                    } else {
                        0
                    };

                    let line = format! {
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        contig_id.clone(),
                        read.get_name().to_string(),
                        read.get_sequence().len(),
                        motif.sequence.to_string(),
                        motif.mod_type.to_pileup_code().to_string(),
                        motif.mod_position,
                        quality,
                    };

                    sender
                        .send(line)
                        .expect("Unable to send line to writer thread");
                }
            }
        }
        Ok(())
    })?;
    drop(sender);
    let _ = writer_handle.join().unwrap();
    Ok(())
}
