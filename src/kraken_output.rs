use std::io::BufRead;

use anyhow::{Context, Result};

/// A single record from kraken2's per-read classification output.
///
/// Each line of the kraken2 output file has five tab-delimited columns:
/// 1. Classification status (`C` = classified, `U` = unclassified)
/// 2. Read name
/// 3. Taxon ID (0 if unclassified)
/// 4. Sequence length(s)
/// 5. K-mer hit list
///
/// Only columns 1-3 are parsed; 4 and 5 are skipped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KrakenRecord {
    classified: bool,
    read_name: String,
    taxon_id: u64,
}

impl KrakenRecord {
    /// Whether this read was classified by kraken2.
    #[must_use]
    pub fn classified(&self) -> bool {
        self.classified
    }

    /// The read name from the kraken2 output.
    #[must_use]
    pub fn read_name(&self) -> &str {
        &self.read_name
    }

    /// The taxon ID assigned to the read (0 if unclassified).
    #[must_use]
    pub fn taxon_id(&self) -> u64 {
        self.taxon_id
    }
}

/// An iterator-like reader over kraken2 per-read classification output.
///
/// Wraps a `BufRead` and yields one [`KrakenRecord`] per line. Reuses an
/// internal line buffer to minimize allocations.
pub struct KrakenOutputReader<R: BufRead> {
    reader: R,
    line_buf: String,
    line_number: u64,
}

impl<R: BufRead> KrakenOutputReader<R> {
    /// Creates a new reader wrapping the given buffered input.
    pub fn new(reader: R) -> Self {
        Self { reader, line_buf: String::new(), line_number: 0 }
    }

    /// Returns the next record, or `None` at end-of-input.
    ///
    /// # Errors
    /// Returns an error if the line cannot be read or parsed.
    pub fn next_record(&mut self) -> Result<Option<KrakenRecord>> {
        self.line_buf.clear();
        let bytes_read = self
            .reader
            .read_line(&mut self.line_buf)
            .context("failed to read kraken output line")?;
        if bytes_read == 0 {
            return Ok(None);
        }
        self.line_number += 1;
        let line = self.line_buf.trim_end_matches(['\n', '\r']);
        parse_kraken_line(line)
            .map(Some)
            .with_context(|| format!("kraken output line {}", self.line_number))
    }

    /// Returns the number of lines read so far.
    #[must_use]
    pub fn line_number(&self) -> u64 {
        self.line_number
    }
}

impl<R: BufRead> Iterator for KrakenOutputReader<R> {
    type Item = Result<KrakenRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(Some(rec)) => Some(Ok(rec)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Parses a single kraken2 output line into a `KrakenRecord`.
fn parse_kraken_line(line: &str) -> Result<KrakenRecord> {
    let mut fields = line.split('\t');

    let status = fields.next().context("missing classification status field")?;
    let classified = match status {
        "C" => true,
        "U" => false,
        other => anyhow::bail!("invalid classification status: {other:?} (expected 'C' or 'U')"),
    };

    let read_name = fields.next().context("missing read name field")?.to_owned();
    anyhow::ensure!(!read_name.is_empty(), "read name is empty");

    let taxon_str = fields.next().context("missing taxon ID field")?;
    let taxon_id: u64 =
        taxon_str.parse().with_context(|| format!("invalid taxon ID: {taxon_str:?}"))?;

    Ok(KrakenRecord { classified, read_name, taxon_id })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_classified_line() {
        let line = "C\tread1\t9606\t150\tA:1 9606:50";
        let rec = parse_kraken_line(line).unwrap();
        assert!(rec.classified());
        assert_eq!(rec.read_name(), "read1");
        assert_eq!(rec.taxon_id(), 9606);
    }

    #[test]
    fn test_parse_unclassified_line() {
        let line = "U\tread2\t0\t150\tA:150";
        let rec = parse_kraken_line(line).unwrap();
        assert!(!rec.classified());
        assert_eq!(rec.read_name(), "read2");
        assert_eq!(rec.taxon_id(), 0);
    }

    #[test]
    fn test_missing_fields() {
        let line = "C\tread1";
        assert!(parse_kraken_line(line).is_err());
    }

    #[test]
    fn test_invalid_status() {
        let line = "X\tread1\t9606\t150\tA:1";
        assert!(parse_kraken_line(line).is_err());
    }

    #[test]
    fn test_invalid_taxon_id() {
        let line = "C\tread1\tnope\t150\tA:1";
        assert!(parse_kraken_line(line).is_err());
    }

    #[test]
    fn test_empty_read_name() {
        let line = "C\t\t9606\t150\tA:1";
        assert!(parse_kraken_line(line).is_err());
    }

    #[test]
    fn test_reader_multiple_records() {
        let input = "C\tread1\t9606\t150\tA:1\nU\tread2\t0\t150\tA:150\n";
        let mut reader = KrakenOutputReader::new(input.as_bytes());

        let r1 = reader.next_record().unwrap().unwrap();
        assert_eq!(r1.read_name(), "read1");
        assert_eq!(r1.taxon_id(), 9606);

        let r2 = reader.next_record().unwrap().unwrap();
        assert_eq!(r2.read_name(), "read2");
        assert_eq!(r2.taxon_id(), 0);

        assert!(reader.next_record().unwrap().is_none());
        assert_eq!(reader.line_number(), 2);
    }

    #[test]
    fn test_reader_empty_input() {
        let mut reader = KrakenOutputReader::new("".as_bytes());
        assert!(reader.next_record().unwrap().is_none());
        assert_eq!(reader.line_number(), 0);
    }
}
