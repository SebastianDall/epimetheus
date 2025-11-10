use anyhow::{Result, anyhow};
use std::ops::{Deref, DerefMut};

use crate::IupacBase;

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Sequence(pub Vec<IupacBase>);

impl Sequence {
    pub fn from_iupac(seq: Vec<IupacBase>) -> Self {
        Self(seq)
    }

    pub fn from_str(sequence_str: &str) -> anyhow::Result<Self> {
        let parsed_sequence = sequence_str
            .chars()
            .map(|b| {
                IupacBase::parse_char(b).map_err(|_| {
                    anyhow::anyhow!(
                        "Base '{}' in sequence '{}' is not a valid IUPAC code",
                        b,
                        sequence_str
                    )
                })
            })
            .collect::<Result<Vec<IupacBase>>>()?;

        Ok(Self::from_iupac(parsed_sequence))
    }

    pub fn from_u8(seq: &[u8]) -> Result<Self> {
        let parsed_sequence: Result<Vec<IupacBase>, anyhow::Error> = seq
            .iter()
            .filter(|&&byte| !byte.is_ascii_whitespace())
            .map(|&byte| {
                IupacBase::from_ascii(byte)
                    .ok_or_else(|| anyhow!("Invalid ascii byte: {}, '{}'", byte, byte as char))
            })
            .collect();

        Ok(Self(parsed_sequence?))
    }
}

impl ToString for Sequence {
    fn to_string(&self) -> String {
        self.into_iter().map(|b| b.to_string()).collect()
    }
}

impl Deref for Sequence {
    type Target = Vec<IupacBase>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Sequence {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl IntoIterator for Sequence {
    type Item = IupacBase;

    type IntoIter = std::vec::IntoIter<IupacBase>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a Sequence {
    type Item = &'a IupacBase;

    type IntoIter = std::slice::Iter<'a, IupacBase>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl FromIterator<IupacBase> for Sequence {
    fn from_iter<T: IntoIterator<Item = IupacBase>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}
