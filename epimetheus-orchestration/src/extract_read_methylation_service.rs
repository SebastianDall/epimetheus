use anyhow::{Context, Result};
use epimetheus_io::io::{
    readers::{bam::BamReader, fastq},
    traits::FastqReader,
};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use methylome::{
    Motif, find_motif_indices_in_sequence,
    read::{Alignment, MethBase},
};
use polars::{df, frame::DataFrame};
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
            "contig_id\tstart_contig\tstrand\tread_id\tread_length\tmapping_quality\tstart_read\tmotif\tmod_type\tmod_position\tquality\tmapping_status"
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
            let sequence_length = sequence.len();
            let modifications = read.get_modifications();
            let mapping = read.get_mapping().unwrap();

            let map_qual = mapping.get_mapping_quality();
            let strand = mapping.get_strand();

            // compute the read mapping from cigar string once.
            let read_mapping: Vec<Option<Alignment>> =
                mapping.build_full_position_map(sequence_length);
            for motif in &motifs {
                let motif_length = motif.sequence.len();
                let indices = find_motif_indices_in_sequence(sequence, &motif);
                for &pos in &indices {
                    let quality = if let Some(meth_base) = modifications.0.get(&pos) {
                        meth_base.quality.0
                    } else {
                        0
                    };

                    let original_pos = match strand {
                        methylome::Strand::Positive => pos,
                        methylome::Strand::Negative => sequence_length - pos - 1,
                    };

                    let genome_pos = match read_mapping.get(original_pos) {
                        Some(Some(Alignment::SequenceMatch(pos))) => *pos as i32,
                        Some(Some(Alignment::SequenceMismatch(pos))) => *pos as i32,
                        Some(Some(Alignment::AmbiguousMatch(pos))) => *pos as i32,
                        _ => -1,
                    };

                    let motif_start_in_bam_coords = match strand {
                        methylome::Strand::Positive => pos - motif.mod_position as usize,
                        methylome::Strand::Negative => original_pos - motif.mod_position as usize,
                    };

                    let alignments: Vec<Option<&Alignment>> = (0..motif_length)
                        .map(|offset| {
                            read_mapping
                                .get(motif_start_in_bam_coords + offset)
                                .and_then(|opt| opt.as_ref())
                        })
                        .collect();

                    let mapping_status = if genome_pos == -1 {
                        "unmapped"
                    } else if alignments
                        .iter()
                        .any(|a| a.is_none() || matches!(a, Some(Alignment::SoftClipped)))
                    {
                        "partial"
                    } else {
                        let positions: Vec<usize> = alignments
                            .iter()
                            .filter_map(|a| match a {
                                Some(Alignment::SequenceMatch(pos))
                                | Some(Alignment::SequenceMismatch(pos))
                                | Some(Alignment::AmbiguousMatch(pos)) => Some(*pos),
                                _ => None,
                            })
                            .collect();

                        if positions.len() != motif_length {
                            "partial"
                        } else if positions.windows(2).all(|w| w[1] == w[0] + 1) {
                            "complete"
                        } else {
                            "gapped"
                        }
                    };

                    let line = format! {
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        contig_id.clone(),
                        genome_pos,
                        strand.to_string(),
                        read.get_name().to_string(),
                        sequence_length,
                        map_qual,
                        pos,
                        motif.sequence.to_string(),
                        motif.mod_type.to_pileup_code().to_string(),
                        motif.mod_position,
                        quality,
                        mapping_status.to_string()
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

pub fn extract_read_methylation_pattern_fastq(
    input_file: &Path,
    read_ids_filter: Option<Vec<String>>,
    motifs: Vec<Motif>,
    threads: usize,
) -> Result<DataFrame> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Could not initialize threadpool");

    let reads = fastq::Reader::read_fastq(input_file, read_ids_filter)?;
    const BATCH_SIZE: usize = 1000;
    let batches: Vec<_> = reads.chunks(BATCH_SIZE).collect();

    // Process batches in parallel
    let results: Vec<(String, u32, u32, String, String, u32, u32)> = batches
        .into_par_iter()
        .map(|batch| {
            let mut batch_data = Vec::new();

            for read in batch {
                let sequence = read.get_sequence();
                let modifications = read.get_modifications();
                let read_length = read.get_sequence().len();

                for motif in &motifs {
                    // Find all motif positions in this read
                    let indices = find_motif_indices_in_sequence(sequence, motif);

                    if !indices.is_empty() {
                        let motif_sequence = motif
                            .sequence
                            .iter()
                            .map(|b| b.to_string())
                            .collect::<String>();

                        for pos in indices {
                            let quality = modifications
                                .0
                                .get(&pos)
                                .unwrap_or(&MethBase::new(
                                    motif.mod_type.clone(),
                                    methylome::read::MethQual(0),
                                ))
                                .clone();
                            let d = (
                                read.get_name().clone(),
                                pos as u32,
                                read_length as u32,
                                motif_sequence.clone(),
                                motif.mod_type.to_pileup_code().to_string(),
                                motif.mod_position as u32,
                                quality.quality.0 as u32,
                            );

                            batch_data.push(d);
                        }
                    }
                }
            }
            batch_data
        })
        .flatten()
        .collect();

    // Merge results from all batches
    // Convert results data to vectors for DataFrame
    let mut read_ids = Vec::with_capacity(results.len());
    let mut starts = Vec::with_capacity(results.len());
    let mut read_lengths = Vec::with_capacity(results.len());
    let mut motif_sequences = Vec::with_capacity(results.len());
    let mut mod_types = Vec::with_capacity(results.len());
    let mut mod_positions = Vec::with_capacity(results.len());
    let mut qualities = Vec::with_capacity(results.len());

    for (read_id, start, read_length, motif_seq, mod_type, mod_pos, quality) in results {
        read_ids.push(read_id);
        starts.push(start);
        read_lengths.push(read_length);
        motif_sequences.push(motif_seq);
        mod_types.push(mod_type);
        mod_positions.push(mod_pos);
        qualities.push(quality);
    }

    // Create DataFrame
    let df = df! [
        "read_id" => read_ids,
        "start" => starts,
        "read_length" => read_lengths,
        "motif_seq" => motif_sequences,
        "mod_type" => mod_types,
        "mod_pos" => mod_positions,
        "quality" => qualities,
    ]?;

    Ok(df)
}
