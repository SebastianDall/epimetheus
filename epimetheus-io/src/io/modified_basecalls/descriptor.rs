use std::str::FromStr;

use anyhow::{Result, anyhow};
use epimetheus_methylome::{IupacBase, Strand};

/// A tag key is like so:
/// [ATGCU][+-]([a-z]+[0-9]+)[.?]?
/// example: A+a.
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct ModifiedBaseDescriptor {
    fundamental_base: IupacBase,
    strand: Strand,
    base_modification_code: ModCode,
    skip_interpreter: Option<SkipInterpreter>,
    cached_string: String,
}

impl ModifiedBaseDescriptor {
    fn compute_cached_string(&self) -> String {
        let cached_string = format!(
            "{}{}{}{}",
            self.fundamental_base.to_string(),
            self.strand.to_string(),
            self.base_modification_code.0,
            self.skip_interpreter
                .map(|i| i.to_string())
                .unwrap_or_default()
        );
        cached_string
    }
    pub fn new(
        base: char,
        strand: char,
        base_modification_code: String,
        skip_interpreter: Option<char>,
    ) -> Result<Self> {
        let fundamental_base = IupacBase::parse_char(base)?;
        match fundamental_base {
            IupacBase::A | IupacBase::G | IupacBase::C | IupacBase::T => {
                // Valid continue
            }
            _ => {
                return Err(anyhow!(
                    "{} is not a valid fundamental base for MM tags.",
                    fundamental_base.to_string()
                ));
            }
        }

        let strand = Strand::try_from(strand)?;
        let base_modification_code = ModCode::new(base_modification_code)?;
        let skip_interpreter = skip_interpreter.and_then(|i| SkipInterpreter::try_from(i).ok());

        let cached_string = format!(
            "{}{}{}{}",
            fundamental_base.to_string(),
            strand.to_string(),
            base_modification_code.0,
            skip_interpreter.map(|i| i.to_string()).unwrap_or_default()
        );

        Ok(Self {
            fundamental_base,
            strand,
            base_modification_code,
            skip_interpreter,
            cached_string,
        })
    }
    pub fn mutate_mod_code(&self, code: ModCode) -> ModifiedBaseDescriptor {
        let mut new_key = self.clone();
        new_key.base_modification_code = code;
        new_key.cached_string = new_key.compute_cached_string();
        new_key
    }

    pub fn as_str(&self) -> &str {
        &self.cached_string
    }
}

impl ToString for ModifiedBaseDescriptor {
    fn to_string(&self) -> String {
        self.cached_string.clone()
    }
}

impl FromStr for ModifiedBaseDescriptor {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s.trim().is_empty() {
            return Err(anyhow!("Cannot parse empty tag key"));
        }

        let mut chars = s.chars();

        let base = chars
            .next()
            .ok_or_else(|| anyhow!("Missing fundamental base in tag key: '{}'", s))?;

        let strand = chars
            .next()
            .ok_or_else(|| anyhow!("Missing strand in tag key: '{}'", s))?;

        let remaining: String = chars.collect();
        let (code, interpreter) = if let Some(last_char) = remaining.chars().last() {
            if last_char == '?' || last_char == '.' {
                let mod_code = &remaining[..remaining.len() - last_char.len_utf8()];
                (mod_code, Some(last_char))
            } else {
                let mod_code = &remaining[..remaining.len()];
                (mod_code, None)
            }
        } else {
            return Err(anyhow!("Missing base modification code: {}", remaining));
        };

        let tag_key = ModifiedBaseDescriptor::new(base, strand, code.to_string(), interpreter)?;
        Ok(tag_key)
    }
}

/// Has to be either lowercase string ([a-z]+|[0-9]+)
#[derive(Hash, Debug, PartialEq, Eq, Clone)]
pub struct ModCode(pub String);
impl ModCode {
    pub fn new(code: String) -> Result<Self> {
        let all_lowercase = code.chars().all(|c| c.is_ascii_lowercase());
        let all_digits = code.chars().all(|c| c.is_numeric());
        if all_lowercase || all_digits {
            Ok(Self(code))
        } else {
            return Err(anyhow!(
                "Invalid mod base code. Must be lowercase string or all digits."
            ));
        }
    }
}

#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy)]
pub enum SkipInterpreter {
    Questionmark,
    Period,
}

impl FromStr for SkipInterpreter {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "?" => Ok(SkipInterpreter::Questionmark),
            "." => Ok(SkipInterpreter::Period),
            _ => return Err(anyhow!("Invalid skip interpreter: {}", s)),
        }
    }
}

impl TryFrom<char> for SkipInterpreter {
    type Error = anyhow::Error;

    fn try_from(value: char) -> std::result::Result<Self, Self::Error> {
        match value {
            '?' => Ok(SkipInterpreter::Questionmark),
            '.' => Ok(SkipInterpreter::Period),
            _ => return Err(anyhow!("Invalid skip interpreter: {}", value)),
        }
    }
}

impl ToString for SkipInterpreter {
    fn to_string(&self) -> String {
        match self {
            SkipInterpreter::Period => ".".to_string(),
            SkipInterpreter::Questionmark => "?".to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_modified_base_descriptor_from_str() {
        let desc = ModifiedBaseDescriptor::from_str("C+m").unwrap();
        assert_eq!(desc.as_str(), "C+m");

        let desc_with_period = ModifiedBaseDescriptor::from_str("A+a.").unwrap();
        assert_eq!(desc_with_period.as_str(), "A+a.");

        let desc_with_question = ModifiedBaseDescriptor::from_str("T+g?").unwrap();
        assert_eq!(desc_with_question.as_str(), "T+g?");

        let desc_numeric = ModifiedBaseDescriptor::from_str("C+21839").unwrap();
        assert_eq!(desc_numeric.as_str(), "C+21839");
    }

    #[test]
    fn test_modified_base_descriptor_invalid() {
        assert!(ModifiedBaseDescriptor::from_str("").is_err());
        assert!(ModifiedBaseDescriptor::from_str("C").is_err());
        assert!(ModifiedBaseDescriptor::from_str("C+").is_err());
        assert!(ModifiedBaseDescriptor::from_str("X+m").is_err()); // Invalid base
    }

    #[test]
    fn test_mod_code_valid() {
        assert!(ModCode::new("m".to_string()).is_ok());
        assert!(ModCode::new("abc".to_string()).is_ok());
        assert!(ModCode::new("123".to_string()).is_ok());
        assert!(ModCode::new("21839".to_string()).is_ok());
    }

    #[test]
    fn test_mod_code_invalid() {
        assert!(ModCode::new("M".to_string()).is_err()); // Uppercase
        assert!(ModCode::new("m1".to_string()).is_err()); // Mixed
        assert!(ModCode::new("a_b".to_string()).is_err()); // Special char
    }
}
