use ahash::AHashMap;
use anyhow::{Result, bail};

use crate::models::{contig::Contig, methylation::MethylationRecord};

pub struct GenomeWorkspaceBuilder {
    workspace: GenomeWorkspace,
}

impl GenomeWorkspaceBuilder {
    pub fn new() -> Self {
        Self {
            workspace: GenomeWorkspace::new(),
        }
    }

    pub fn add_contig(&mut self, contig: Contig) -> Result<&mut Self> {
        if self.workspace.contigs.contains_key(&contig.id) {
            bail!("Key error: '{}' already inserted", &contig.id)
        }

        self.workspace.contigs.insert(contig.id.clone(), contig);
        Ok(self)
    }

    #[allow(dead_code)]
    pub fn add_record(&mut self, record: MethylationRecord) -> Result<&mut Self> {
        if let Some(contig_entry) = self.workspace.get_mut_contig(&record.get_contig_id()) {
            contig_entry.add_methylation(
                record.position,
                record.strand,
                record.mod_type,
                record.methylation,
            )?;
        } else {
            bail!(
                "Warning: Contig: '{}' found in pileup, but not in assembly",
                record.contig
            );
        };
        Ok(self)
    }

    pub fn build(self) -> GenomeWorkspace {
        self.workspace
    }
}

pub struct GenomeWorkspace {
    contigs: AHashMap<String, Contig>,
}

impl GenomeWorkspace {
    fn new() -> Self {
        Self {
            contigs: AHashMap::new(),
        }
    }
    pub fn get_workspace(&self) -> AHashMap<String, Contig> {
        self.contigs.clone()
    }

    fn get_mut_contig(&mut self, id: &str) -> Option<&mut Contig> {
        self.contigs.get_mut(id)
    }

    pub fn is_empty(&self) -> bool {
        self.contigs.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use crate::models::methylation::MethylationCoverage;
    use crate::models::pileup::{PileupRecord, PileupRecordString};

    use super::*;
    use anyhow::Result;
    use methylome::{ModType, Strand};
    use std::io::BufRead;
    use std::str::FromStr;
    use std::{
        fs::File,
        io::{BufReader, Write},
    };
    use tempfile::NamedTempFile;

    #[test]
    fn test_strand_from_str() -> Result<()> {
        // Mock pileup data lines
        let pileup_data = vec![
            "contig_3\t0\t1\tm\t133\t-\t0\t1\t255,0,0\t133\t0.00\t0\t133\t0\t0\t6\t0\t0",
            "contig_3\t1\t2\ta\t174\t+\t1\t2\t255,0,0\t174\t1.72\t3\t171\t0\t0\t3\t0\t0",
        ];

        // Expected results for the strand column
        let expected_strands = vec![Strand::Negative, Strand::Positive];

        // Iterate through pileup data and validate strand parsing
        for (line, &expected_strand) in pileup_data.iter().zip(expected_strands.iter()) {
            let fields: Vec<&str> = line.split('\t').collect();
            let strand_field = fields[5]; // Extract strand field
            let strand = Strand::from_str(strand_field)?;

            // Assert that the parsed strand matches the expected value
            assert_eq!(strand, expected_strand);
        }

        Ok(())
    }

    #[test]
    fn test_populate_methylation() -> Result<()> {
        let mut workspace_builder = GenomeWorkspaceBuilder::new();

        // Add a mock contig to the workspace
        workspace_builder
            .add_contig(Contig::from_string("contig_3".to_string(), "ATCG".to_string()).unwrap())?;

        // Create a temporary pileup file
        let mut pileup_file = NamedTempFile::new()?;
        writeln!(
            pileup_file,
            "contig_3\t0\t1\tm\t133\t-\t0\t1\t255,0,0\t133\t0.00\t10\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t1\t2\ta\t174\t+\t1\t2\t255,0,0\t174\t1.72\t5\t169\t0\t0\t3\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t2\t3\ta\t172\t+\t2\t3\t255,0,0\t0\t0.00\t0\t0\t0\t0\t0\t0\t0" // Zero coverage, should be skipped
        )?;

        // Populate methylation data
        let file = File::open(pileup_file)?;
        let reader = BufReader::new(file);

        for res in reader.lines() {
            let record = res.unwrap();
            let pileup_record = PileupRecord::try_from(PileupRecordString::new(record)).unwrap();
            let meth_record = MethylationRecord::try_from_with_filters(pileup_record, 3, 0.8);

            let meth = match meth_record {
                Ok(Some(m)) => m,
                Ok(None) => continue,
                Err(e) => return Err(e),
            };
            workspace_builder.add_record(meth).unwrap();
        }

        let mut workspace = workspace_builder.build();

        // Get the contig
        let contig = workspace.get_mut_contig("contig_3").unwrap();

        // Check that methylated_positions are correctly populated
        assert_eq!(contig.methylated_positions.len(), 2);

        // Validate individual positions
        let pos_0 = contig
            .methylated_positions
            .get(&(0, Strand::Negative, ModType::FiveMC));

        let expected_methylation = MethylationCoverage::new(10, 133, 0).unwrap();

        assert!(pos_0.is_some());
        assert_eq!(pos_0, Some(&expected_methylation));

        let pos_1 = contig
            .methylated_positions
            .get(&(1, Strand::Positive, ModType::SixMA));
        let expected_methylation = MethylationCoverage::new(5, 174, 0).unwrap();
        assert!(pos_1.is_some());
        assert_eq!(pos_1, Some(&expected_methylation));

        // Ensure position with zero coverage is skipped
        let pos_2 = contig
            .methylated_positions
            .get(&(2, Strand::Positive, ModType::SixMA));
        assert!(pos_2.is_none());

        Ok(())
    }

    #[test]
    fn test_populate_methylation_missing_contig() {
        let mut workspace_builder = GenomeWorkspaceBuilder::new();
        // Create a temporary pileup file
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t0\t1\tm\t133\t-\t0\t1\t255,0,0\t133\t0.00\t10\t123\t0\t0\t6\t0\t0"
        )
        .unwrap();

        // Populate methylation data
        let file = File::open(pileup_file.path()).unwrap();
        let reader = BufReader::new(file);

        for res in reader.lines() {
            let record = res.unwrap();
            let pileup_record = PileupRecord::try_from(PileupRecordString::new(record)).unwrap();
            let meth_record = MethylationRecord::try_from_with_filters(pileup_record, 3, 0.8)
                .unwrap()
                .unwrap();

            let result = workspace_builder.add_record(meth_record);
            assert!(result.is_err());
        }
    }
}
