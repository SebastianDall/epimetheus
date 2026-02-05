use ahash::HashMap;
use anyhow::anyhow;
use clap::Parser;
use epimetheus_io::io::readers::bam::{ModCode, TagKey};
use epimetheus_orchestration::bam_tag_merge_service::BamMergeArgs;
use std::{path::PathBuf, str::FromStr};

#[derive(Parser, Debug, Clone)]
pub struct BamMergeCliArgs {
    #[arg(short, long, help = "paths to bams to extract tags from")]
    pub from_bam: PathBuf,

    #[arg(short, long, help = "path to bam to add tags to. file must exist.")]
    pub to_bam: PathBuf,

    #[arg(short, long, help = "path to construct db with tags for 'from_bams'.")]
    pub db_path: PathBuf,

    #[arg(required = false, num_args(1..), long, help = "Rename mod code in tag key. <from:to>. ex C+21839.:m. Will change 21839 to m in the tag code. Has to be a lowercase string or a sequence of numbers.")]
    pub rename_tags: Vec<String>,

    #[arg(required = false, num_args(1..), long, help = "Mod codes to be ignored in the from bam.")]
    pub ignore_tags_from_bam: Vec<String>,

    #[arg(long, default_value_t = false, help = "keep database to be reused.")]
    pub keep_db: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "leave original out bam untouched."
    )]
    pub keep_outfile: bool,
}

impl BamMergeCliArgs {
    pub fn construct_rename_tags(&self) -> anyhow::Result<HashMap<TagKey, TagKey>> {
        self.rename_tags.iter().map(|t_str| {
            let parts: Vec<&str> = t_str.split(":").collect();
            if parts.len() != 2 {
                return Err(anyhow!("Invalid modification code renaming. Should be <from_code>:<to_code> got: {}", t_str));
            }

            let from_key = TagKey::from_str(parts[0])?;
            let new_mod_code = ModCode::new(parts[1].to_string())?;
            let to_key = from_key.mutate_mod_code(new_mod_code);

            Ok((from_key, to_key))
        }).collect::<Result<HashMap<TagKey, TagKey>, _>>()
    }
}

impl TryFrom<BamMergeCliArgs> for BamMergeArgs {
    type Error = anyhow::Error;

    fn try_from(cli: BamMergeCliArgs) -> Result<Self, Self::Error> {
        let rename_tags = cli.construct_rename_tags()?;
        Ok(Self {
            from_bam: cli.from_bam,
            to_bam: cli.to_bam,
            db_path: cli.db_path,
            rename_tags: rename_tags,
            ignore_tags_from_bam: cli.ignore_tags_from_bam,
            keep_db: cli.keep_db,
            keep_outfile: cli.keep_outfile,
        })
    }
}
