use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

/// The standard taxonomic ranks used in kraken2 reports.
///
/// Each variant maps to a single-character rank code used in kraken2 output.
/// Non-standard intermediate ranks (e.g. "G2") are represented by pairing
/// a `Rank` with a depth in [`TaxonomicRank`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Rank {
    /// Unclassified sequences (U)
    Unclassified,
    /// Root of the taxonomy tree (R)
    Root,
    /// Domain / superkingdom (D)
    Domain,
    /// Kingdom (K)
    Kingdom,
    /// Phylum (P)
    Phylum,
    /// Class (C)
    Class,
    /// Order (O)
    Order,
    /// Family (F)
    Family,
    /// Genus (G)
    Genus,
    /// Species (S)
    Species,
}

impl Rank {
    /// Returns the single-character code for this rank.
    #[must_use]
    pub fn code(self) -> char {
        match self {
            Self::Unclassified => 'U',
            Self::Root => 'R',
            Self::Domain => 'D',
            Self::Kingdom => 'K',
            Self::Phylum => 'P',
            Self::Class => 'C',
            Self::Order => 'O',
            Self::Family => 'F',
            Self::Genus => 'G',
            Self::Species => 'S',
        }
    }
}

impl fmt::Display for Rank {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.code())
    }
}

impl FromStr for Rank {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "U" => Ok(Self::Unclassified),
            "R" => Ok(Self::Root),
            "D" => Ok(Self::Domain),
            "K" => Ok(Self::Kingdom),
            "P" => Ok(Self::Phylum),
            "C" => Ok(Self::Class),
            "O" => Ok(Self::Order),
            "F" => Ok(Self::Family),
            "G" => Ok(Self::Genus),
            "S" => Ok(Self::Species),
            _ => anyhow::bail!("invalid rank code: {s:?}"),
        }
    }
}

/// A taxonomic rank with an optional depth suffix for non-standard intermediate ranks.
///
/// Standard ranks (e.g. Genus) have `depth: None` and display as a single character ("G").
/// Non-standard ranks (e.g. two levels below Genus) have `depth: Some(2)` and display
/// as the rank code followed by the depth ("G2").
///
/// Kraken2 increments the depth suffix during DFS traversal for each non-standard rank
/// encountered between two standard ranks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct TaxonomicRank {
    rank: Rank,
    depth: Option<u32>,
}

impl TaxonomicRank {
    /// Returns the base rank.
    #[must_use]
    pub fn rank(&self) -> Rank {
        self.rank
    }

    /// Returns the depth suffix, or `None` for standard ranks.
    #[must_use]
    pub fn depth(&self) -> Option<u32> {
        self.depth
    }

    /// Returns `true` if this is a standard rank (no depth suffix).
    #[must_use]
    pub fn is_standard(&self) -> bool {
        self.depth.is_none()
    }
}

impl fmt::Display for TaxonomicRank {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.rank.code())?;
        if let Some(d) = self.depth {
            write!(f, "{d}")?;
        }
        Ok(())
    }
}

impl FromStr for TaxonomicRank {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        anyhow::ensure!(!s.is_empty(), "empty rank string");

        let code = &s[..1];
        let rank = Rank::from_str(code)?;
        let suffix = &s[1..];

        if suffix.is_empty() {
            return Ok(Self { rank, depth: None });
        }

        let depth: u32 =
            suffix.parse().map_err(|_| anyhow::anyhow!("invalid rank depth suffix: {s:?}"))?;
        anyhow::ensure!(depth > 0, "rank depth suffix must be > 0, got: {s:?}");
        Ok(Self { rank, depth: Some(depth) })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rank_round_trip_all_standard() {
        let cases = [
            (Rank::Unclassified, "U"),
            (Rank::Root, "R"),
            (Rank::Domain, "D"),
            (Rank::Kingdom, "K"),
            (Rank::Phylum, "P"),
            (Rank::Class, "C"),
            (Rank::Order, "O"),
            (Rank::Family, "F"),
            (Rank::Genus, "G"),
            (Rank::Species, "S"),
        ];

        for (rank, code) in cases {
            assert_eq!(rank.to_string(), code, "Display failed for {rank:?}");
            assert_eq!(Rank::from_str(code).unwrap(), rank, "FromStr failed for {code}");
            assert_eq!(rank.code(), code.chars().next().unwrap());
        }
    }

    #[test]
    fn test_rank_invalid() {
        assert!(Rank::from_str("X").is_err());
        assert!(Rank::from_str("").is_err());
        assert!(Rank::from_str("GG").is_err());
    }

    #[test]
    fn test_taxonomic_rank_standard_round_trip() {
        for code in ["U", "R", "D", "K", "P", "C", "O", "F", "G", "S"] {
            let tr = TaxonomicRank::from_str(code).unwrap();
            assert!(tr.is_standard());
            assert_eq!(tr.depth(), None);
            assert_eq!(tr.to_string(), code);
        }
    }

    #[test]
    fn test_taxonomic_rank_non_standard_parse() {
        let tr = TaxonomicRank::from_str("G2").unwrap();
        assert_eq!(tr.rank(), Rank::Genus);
        assert_eq!(tr.depth(), Some(2));
        assert!(!tr.is_standard());
        assert_eq!(tr.to_string(), "G2");
    }

    #[test]
    fn test_taxonomic_rank_domain_depth() {
        let tr = TaxonomicRank::from_str("D1").unwrap();
        assert_eq!(tr.rank(), Rank::Domain);
        assert_eq!(tr.depth(), Some(1));
    }

    #[test]
    fn test_taxonomic_rank_display_non_standard() {
        let tr = TaxonomicRank::from_str("S3").unwrap();
        assert_eq!(tr.rank(), Rank::Species);
        assert_eq!(tr.depth(), Some(3));
        assert_eq!(tr.to_string(), "S3");
    }

    #[test]
    fn test_taxonomic_rank_invalid_code() {
        assert!(TaxonomicRank::from_str("X").is_err());
        assert!(TaxonomicRank::from_str("X1").is_err());
    }

    #[test]
    fn test_taxonomic_rank_zero_depth() {
        assert!(TaxonomicRank::from_str("G0").is_err());
    }

    #[test]
    fn test_taxonomic_rank_empty() {
        assert!(TaxonomicRank::from_str("").is_err());
    }

    #[test]
    fn test_taxonomic_rank_non_numeric_suffix() {
        assert!(TaxonomicRank::from_str("Gx").is_err());
    }
}
