use ahash::HashMap;
use anyhow::anyhow;
use clap::Parser;
use epimetheus_io::io::modified_basecalls::descriptor::{ModCode, ModifiedBaseDescriptor};
use epimetheus_orchestration::bam_tag_merge_service::{BamMergeArgs, FromTags};
use std::{path::PathBuf, str::FromStr};

#[derive(Parser, Debug, Clone)]
pub struct BamMergeCliArgs {
    #[arg(short, long, help = "Path to bam to extract tags from")]
    pub from_bam: Option<PathBuf>,

    #[arg(
        short,
        long,
        help = "Path to db where tags have been extracted. This can not be used with '--from-bam'"
    )]
    pub from_db: Option<PathBuf>,

    #[arg(short, long, help = "path to bam to add tags to. file must exist.")]
    pub to_bam: PathBuf,

    #[arg(
        short,
        required = false,
        long,
        help = "path to construct db with tags for 'from_bams'."
    )]
    pub db_path: Option<PathBuf>,

    #[arg(required = false, num_args(1..), long, help = "Rename mod code in tag key in the 'from_bam'. <from:to>. ex C+21839.:m. Will change 21839 to m in the tag code. Has to be a lowercase string or a sequence of numbers.")]
    pub rename_tags_from_bam: Vec<String>,

    #[arg(required = false, num_args(1..), long, help = "Rename mod code in tag key in the 'in_bam'. <from:to>. ex C+21839.:m. Will change 21839 to m in the tag code. Has to be a lowercase string or a sequence of numbers.")]
    pub rename_tags_to_bam: Vec<String>,

    #[arg(required = false, num_args(1..), long, help = "Mod codes to be ignored in the from bam.")]
    pub ignore_tags_from_bam: Vec<String>,

    #[arg(
        long,
        default_value_t = false,
        help = "keep database to be reused. If user specified path '--db-path' this is ignored."
    )]
    pub keep_db: bool,
}

impl BamMergeCliArgs {
    pub fn validate_from_arguments(&self) -> anyhow::Result<FromTags> {
        match (&self.from_bam, &self.from_db) {
            (Some(b), None) => Ok(FromTags::Bam(b.to_path_buf())),
            (None, Some(db)) => Ok(FromTags::Db(db.to_path_buf())),
            (Some(_), Some(_)) => Err(anyhow!(
                "You must specify either '--from-bam' or '--from-db'. Not both"
            )),
            (None, None) => Err(anyhow!(
                "You must specify either '--from-bam' or '--from-db'."
            )),
        }
    }

    pub fn validate_rename_tags(
        tags: &Vec<String>,
    ) -> anyhow::Result<HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>> {
        tags.iter().map(|t_str| {
            let parts: Vec<&str> = t_str.split(":").collect();
            if parts.len() != 2 {
                return Err(anyhow!("Invalid modification code renaming. Should be <from_code>:<to_code> got: {}", t_str));
            }

            let from_key = ModifiedBaseDescriptor::from_str(parts[0])?;
            let new_mod_code = ModCode::new(parts[1].to_string())?;
            let to_key = from_key.mutate_mod_code(new_mod_code);

            Ok((from_key, to_key))
        }).collect::<Result<HashMap<ModifiedBaseDescriptor, ModifiedBaseDescriptor>, _>>()
    }
}

impl TryFrom<BamMergeCliArgs> for BamMergeArgs {
    type Error = anyhow::Error;

    fn try_from(cli: BamMergeCliArgs) -> Result<Self, Self::Error> {
        let from = cli.validate_from_arguments()?;

        let rename_tags_from = BamMergeCliArgs::validate_rename_tags(&cli.rename_tags_from_bam)?;
        let rename_tags_from_opt = if !rename_tags_from.is_empty() {
            Some(rename_tags_from)
        } else {
            None
        };

        let rename_tags_to = BamMergeCliArgs::validate_rename_tags(&cli.rename_tags_to_bam)?;
        let rename_tags_to_opt = if !rename_tags_to.is_empty() {
            Some(rename_tags_to)
        } else {
            None
        };

        let ignore_tags_from_bam = if !cli.ignore_tags_from_bam.is_empty() {
            let tag = cli
                .ignore_tags_from_bam
                .iter()
                .map(|s| ModifiedBaseDescriptor::from_str(s.as_str()))
                .collect::<Result<Vec<_>, _>>()?;

            Some(tag)
        } else {
            None
        };

        Ok(Self {
            from,
            to_bam: cli.to_bam,
            db_path: cli.db_path,
            rename_tags_to_bam: rename_tags_to_opt,
            rename_tags_from_bam: rename_tags_from_opt,
            ignore_tags_from_bam,
            keep_db: cli.keep_db,
        })
    }
}
