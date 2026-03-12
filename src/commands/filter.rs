use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use fgoxide::io::Io;
use fgoxide::iter::IntoChunkedReadAheadIterator;
use pooled_writer::bgzf::BgzfCompressor;
use pooled_writer::{PoolBuilder, PooledWriter};
use seq_io::fastq::{Error as FastqError, OwnedRecord, Reader as FastqReader, Record};

use crate::commands::command::Command;
use crate::kraken_output::{KrakenOutputReader, KrakenRecord};
use crate::progress::{ProgressLogger, format_count};
use crate::report::KrakenReport;

/// Number of records per chunk sent through the read-ahead channel.
const READ_AHEAD_CHUNK_SIZE: usize = 1024;

/// Number of buffered chunks in the read-ahead channel.
const READ_AHEAD_NUM_CHUNKS: usize = 1024;

/// Buffer size used when opening input files for reading.
const IO_BUFFER_SIZE: usize = 512 * 1024;

/// Filter reads from FASTQ files based on kraken2 classification results.
///
/// Extracts reads classified to one or more taxon IDs from FASTQ files, using the
/// kraken2 report (taxonomy tree) and per-read classification output. Supports both
/// single-end and paired-end reads, and writes bgzf-compressed output.
///
/// # Required inputs
///
/// The command needs three pieces of data that must all come from the same kraken2 run:
///
/// - **`--kraken-report`** (`-r`): The kraken2 report file containing the taxonomy tree
///   and per-taxon read counts. This is used to resolve taxon IDs, expand descendants,
///   and estimate the expected number of matching reads.
/// - **`--kraken-output`** (`-k`): The per-read classification output from kraken2
///   (generated with `--output`). Each line maps a read name to a taxon ID.
/// - **`--input`** (`-i`): One FASTQ file for single-end data, or two for paired-end.
///   Gzip and bgzf compressed inputs are detected and handled automatically.
///
/// The kraken output and FASTQ file(s) must contain the same reads in the same order.
/// The command verifies read name agreement and will error if the files are mismatched
/// or have different numbers of records.
///
/// # Taxon selection
///
/// At least one of `--taxon-ids` or `--include-unclassified` must be specified.
///
/// - **`--taxon-ids`** (`-t`): One or more NCBI taxon IDs to extract. By default, only
///   reads classified directly to these exact taxon IDs are included.
/// - **`--include-descendants`** (`-d`): Expand each taxon ID to include all of its
///   descendants in the taxonomy tree. For example, specifying a genus-level taxon ID
///   with `-d` will also extract reads classified to any species or strain within that
///   genus.
/// - **`--include-unclassified`** (`-u`): Include reads that kraken2 could not classify
///   (taxon ID 0). Can be combined with `--taxon-ids` to extract both classified and
///   unclassified reads in a single pass.
///
/// # Output
///
/// - **`--output`** (`-o`): Output FASTQ file path(s). Must provide the same number of
///   output files as input files (one for single-end, two for paired-end). Outputs are
///   always bgzf-compressed regardless of file extension.
/// - **`--threads`**: Number of threads used for bgzf compression (default: 4).
/// - **`--compression-level`**: Bgzf compression level from 0 (fastest) to 9 (smallest),
///   default 5.
///
/// # Examples
///
/// Extract all reads classified as _E. coli_ (taxon 562):
///
/// ```bash
/// k2tools filter -r report.txt -k output.txt -i reads.fq.gz -o ecoli.fq.gz -t 562
/// ```
///
/// Extract all Enterobacteriaceae (taxon 543) including every species and strain beneath
/// it in the taxonomy:
///
/// ```bash
/// k2tools filter -r report.txt -k output.txt \
///     -i reads.fq.gz -o entero.fq.gz -t 543 -d
/// ```
///
/// Extract unclassified reads from a paired-end run:
///
/// ```bash
/// k2tools filter -r report.txt -k output.txt \
///     -i r1.fq.gz r2.fq.gz -o unclass_r1.fq.gz unclass_r2.fq.gz -u
/// ```
///
/// Extract human reads plus unclassified in a single pass:
///
/// ```bash
/// k2tools filter -r report.txt -k output.txt \
///     -i reads.fq.gz -o host_and_unclass.fq.gz -t 9606 -d -u
/// ```
#[derive(clap::Args)]
pub struct Filter {
    /// Path to the kraken2 report file.
    #[arg(short = 'r', long)]
    kraken_report: PathBuf,

    /// Path to the kraken2 per-read classification output.
    #[arg(short = 'k', long)]
    kraken_output: PathBuf,

    /// Input FASTQ file(s). One for single-end, two for paired-end.
    /// Supports gzip/bgzf compressed inputs.
    #[arg(short, long, num_args = 1..=2, required = true)]
    input: Vec<PathBuf>,

    /// Output FASTQ file(s). Must match the number of inputs.
    /// Written with bgzf compression.
    #[arg(short, long, num_args = 1..=2, required = true)]
    output: Vec<PathBuf>,

    /// Taxon ID(s) to extract reads for. At least one taxon ID or
    /// --include-unclassified must be specified.
    #[arg(short, long, num_args = 1..)]
    taxon_ids: Vec<u64>,

    /// Include reads assigned to any descendant of the specified taxa.
    #[arg(short = 'd', long, default_value_t = false)]
    include_descendants: bool,

    /// Include unclassified reads (taxon ID 0) in the output.
    #[arg(short = 'u', long, default_value_t = false)]
    include_unclassified: bool,

    /// Number of threads for bgzf compression.
    #[arg(long, default_value_t = 4)]
    threads: usize,

    /// Bgzf compression level (0-9).
    #[arg(long, default_value_t = 5)]
    compression_level: u8,
}

impl Command for Filter {
    fn execute(&self) -> Result<()> {
        self.validate_args()?;

        let report = KrakenReport::from_path(&self.kraken_report)?;
        if report.is_empty() {
            return self.handle_empty_inputs();
        }

        let (taxon_set, expected) = build_taxon_set_and_expected_count(
            &report,
            &self.taxon_ids,
            self.include_descendants,
            self.include_unclassified,
        )?;
        log::info!(
            "Filtering for {} taxa; expecting approximately {} reads",
            format_count(taxon_set.len() as u64),
            format_count(expected),
        );

        let (total, kept) = self.run_filter_pipeline(&taxon_set).map_err(|e| {
            let banner = "#".repeat(72);
            let output_paths: Vec<_> =
                self.output.iter().map(|p| format!("  {}", p.display())).collect();
            eprintln!(
                "\n{banner}\n\
                 # ERROR: invalid inputs detected\n\
                 #\n\
                 # {e}\n\
                 #\n\
                 # WARNING: partial/invalid output files may have been written to:\n\
                 # {}\n\
                 {banner}\n",
                output_paths.join("\n"),
            );
            e
        })?;

        #[allow(clippy::cast_precision_loss)]
        let pct = if total > 0 { kept as f64 / total as f64 * 100.0 } else { 0.0 };
        log::info!(
            "Kept {} / {} reads ({pct:.2}%), expected {}.",
            format_count(kept),
            format_count(total),
            format_count(expected),
        );

        Ok(())
    }
}

impl Filter {
    /// Handles the case where kraken2 was run on empty FASTQ files, producing an
    /// empty report and no kraken output file. Verifies that all FASTQ inputs are
    /// truly empty, then writes valid empty bgzf output files.
    fn handle_empty_inputs(&self) -> Result<()> {
        let io = Io::new(u32::from(self.compression_level), IO_BUFFER_SIZE);
        for path in &self.input {
            let reader = io
                .new_reader(path)
                .with_context(|| format!("failed to open FASTQ: {}", path.display()))?;
            let mut fq = FastqReader::new(reader);
            if fq.next().is_some() {
                anyhow::bail!(
                    "kraken2 report is empty but FASTQ input {} contains records; \
                     inputs are inconsistent",
                    path.display()
                );
            }
        }

        let (mut pool, writers) = self.build_writer_pool()?;
        for w in writers {
            w.close()?;
        }
        pool.stop_pool()?;

        log::info!("Report is empty; all inputs are empty. Wrote empty output files.");
        Ok(())
    }

    /// Validates command-line arguments beyond what clap enforces.
    fn validate_args(&self) -> Result<()> {
        anyhow::ensure!(
            self.input.len() == self.output.len(),
            "number of input files ({}) must match number of output files ({})",
            self.input.len(),
            self.output.len()
        );
        anyhow::ensure!(self.threads >= 1, "threads must be at least 1");
        anyhow::ensure!(self.compression_level <= 9, "compression level must be 0-9");
        anyhow::ensure!(
            !self.taxon_ids.is_empty() || self.include_unclassified,
            "at least one --taxon-ids value or --include-unclassified must be specified"
        );
        Ok(())
    }

    /// Opens all inputs, creates writers, runs the main filter loop, and closes
    /// everything down. Returns (total_reads, kept_reads).
    fn run_filter_pipeline(&self, taxon_set: &HashSet<u64>) -> Result<(u64, u64)> {
        let io = Io::new(u32::from(self.compression_level), IO_BUFFER_SIZE);
        let kraken_reader = io.new_reader(&self.kraken_output).with_context(|| {
            format!("failed to open kraken output: {}", self.kraken_output.display())
        })?;
        let mut kraken_iter = KrakenOutputReader::new(kraken_reader)
            .read_ahead(READ_AHEAD_CHUNK_SIZE, READ_AHEAD_NUM_CHUNKS);

        let is_paired = self.input.len() == 2;
        let mut fq_iter1 = FastqReader::new(
            io.new_reader(&self.input[0])
                .with_context(|| format!("failed to open FASTQ: {}", self.input[0].display()))?,
        )
        .into_records()
        .read_ahead(READ_AHEAD_CHUNK_SIZE, READ_AHEAD_NUM_CHUNKS);

        let mut fq_iter2 = if is_paired {
            Some(
                FastqReader::new(io.new_reader(&self.input[1]).with_context(|| {
                    format!("failed to open FASTQ: {}", self.input[1].display())
                })?)
                .into_records()
                .read_ahead(READ_AHEAD_CHUNK_SIZE, READ_AHEAD_NUM_CHUNKS),
            )
        } else {
            None
        };

        let (mut pool, mut writers) = self.build_writer_pool()?;
        let mut progress = ProgressLogger::new("k2tools::filter", "reads", 5_000_000);

        // Run the filter and verification, capturing any error so we can
        // shut down the pool cleanly before propagating it (avoids panics
        // in PooledWriter::drop when writers outlive the pool).
        let result = filter_reads(
            &mut kraken_iter,
            &mut fq_iter1,
            fq_iter2.as_mut(),
            taxon_set,
            &mut writers,
            &mut progress,
        )
        .and_then(|(total, kept)| {
            verify_fastq_exhausted(&mut fq_iter1, fq_iter2.as_mut(), total)?;
            Ok((total, kept))
        });

        progress.finish();

        // Always close writers before stopping the pool
        for w in writers {
            w.close()?;
        }
        pool.stop_pool()?;

        result
    }

    /// Constructs the bgzf writer pool and exchanges output files into pooled writers.
    /// Returns (pool, writers) so that destructuring as `let (pool, writers) = ...`
    /// ensures writers are dropped before the pool (reverse declaration order).
    fn build_writer_pool(&self) -> Result<(pooled_writer::Pool, Vec<PooledWriter>)> {
        let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
            .threads(self.threads)
            .queue_size(self.threads * 50)
            .compression_level(self.compression_level)?;

        let mut writers: Vec<PooledWriter> = Vec::new();
        for path in &self.output {
            let file = File::create(path)
                .with_context(|| format!("failed to create output: {}", path.display()))?;
            writers.push(pool_builder.exchange(BufWriter::new(file)));
        }
        let pool = pool_builder.build()?;
        Ok((pool, writers))
    }
}

/// Runs the main filter loop: co-iterates kraken output and FASTQ iterator(s) in
/// lockstep, writing matching records to the output writers.
///
/// Returns (total_reads_processed, reads_kept).
fn filter_reads(
    kraken_iter: &mut impl Iterator<Item = Result<KrakenRecord>>,
    fq_iter1: &mut impl Iterator<Item = Result<OwnedRecord, FastqError>>,
    mut fq_iter2: Option<&mut impl Iterator<Item = Result<OwnedRecord, FastqError>>>,
    taxon_set: &HashSet<u64>,
    writers: &mut [PooledWriter],
    progress: &mut ProgressLogger,
) -> Result<(u64, u64)> {
    let mut total: u64 = 0;
    let mut kept: u64 = 0;

    for kraken_result in kraken_iter {
        let kraken_rec = kraken_result?;
        total += 1;
        progress.record();

        let fq_rec1 = fq_iter1
            .next()
            .context("FASTQ input ended before kraken output")?
            .with_context(|| format!("failed to read FASTQ record at kraken line {total}"))?;

        let fq_rec2: Option<OwnedRecord> = if let Some(ref mut iter2) = fq_iter2 {
            Some(
                iter2
                    .next()
                    .context("second FASTQ input ended before kraken output")?
                    .with_context(|| {
                        format!("failed to read FASTQ R2 record at kraken line {total}")
                    })?,
            )
        } else {
            None
        };

        if taxon_set.contains(&kraken_rec.taxon_id()) {
            // Validate read names only for matching reads to avoid overhead
            validate_read_name(kraken_rec.read_name(), fq_rec1.head(), total)?;
            if let Some(ref rec2) = fq_rec2 {
                validate_read_name(kraken_rec.read_name(), rec2.head(), total)?;
            }

            write_fastq_record(&mut writers[0], &fq_rec1)?;
            if let Some(ref rec2) = fq_rec2 {
                write_fastq_record(&mut writers[1], rec2)?;
            }
            kept += 1;
        }
    }

    Ok((total, kept))
}

/// Verifies that the FASTQ streams are exhausted after the kraken output ends.
fn verify_fastq_exhausted(
    fq_iter1: &mut impl Iterator<Item = Result<OwnedRecord, FastqError>>,
    fq_iter2: Option<&mut impl Iterator<Item = Result<OwnedRecord, FastqError>>>,
    total: u64,
) -> Result<()> {
    if fq_iter1.next().is_some() {
        anyhow::bail!("FASTQ input has more records than kraken output ({total} kraken records)");
    }
    if let Some(iter2) = fq_iter2 {
        if iter2.next().is_some() {
            anyhow::bail!(
                "second FASTQ input has more records than kraken output ({total} kraken records)"
            );
        }
    }
    Ok(())
}

/// Builds the set of taxon IDs to filter for and computes the expected number of
/// matching reads from the report's count fields.
///
/// If `include_descendants` is true, expands each taxon ID to include all its
/// descendants in the report taxonomy tree. If `include_unclassified` is true,
/// adds taxon ID 0. The expected count uses `clade_count` when descendants are
/// included, `direct_count` otherwise.
///
/// Returns `(taxon_id_set, expected_read_count)`.
fn build_taxon_set_and_expected_count(
    report: &KrakenReport,
    taxon_ids: &[u64],
    include_descendants: bool,
    include_unclassified: bool,
) -> Result<(HashSet<u64>, u64)> {
    let mut set = HashSet::new();
    let mut expected: u64 = 0;

    for &tid in taxon_ids {
        let idx = report
            .index_of_taxon_id(tid)
            .with_context(|| format!("taxon ID {tid} not found in report"))?;
        let row = report.row(idx);
        set.insert(tid);

        if include_descendants {
            expected += row.clade_count();
            for desc_idx in report.descendants(idx) {
                set.insert(report.row(desc_idx).taxon_id());
            }
        } else {
            expected += row.direct_count();
        }
    }

    if include_unclassified {
        set.insert(0);
        if let Some(row) = report.get_by_taxon_id(0) {
            expected += row.clade_count();
        }
    }

    Ok((set, expected))
}

/// Validates that a kraken read name matches a FASTQ record header.
///
/// Expects the FASTQ header to start with the kraken read name (byte-for-byte),
/// optionally followed by `/1` or `/2` (paired-end suffix) and/or whitespace
/// plus a comment. Avoids scanning the full header — only checks the prefix
/// at the kraken name length boundary.
fn validate_read_name(kraken_name: &str, fastq_head: &[u8], line_number: u64) -> Result<()> {
    let k = kraken_name.as_bytes();
    let f = fastq_head;

    if f.len() >= k.len() && f[..k.len()] == *k {
        let rest = &f[k.len()..];
        if rest.is_empty()
            || rest[0] == b' '
            || rest[0] == b'\t'
            || (rest.len() >= 2
                && rest[0] == b'/'
                && (rest[1] == b'1' || rest[1] == b'2')
                && (rest.len() == 2 || rest[2] == b' ' || rest[2] == b'\t'))
        {
            return Ok(());
        }
    }

    // Build a readable FASTQ name for the error message only on failure
    let name_end = f.iter().position(|&b| b == b' ' || b == b'\t').unwrap_or(f.len());
    anyhow::bail!(
        "read name mismatch at kraken line {line_number}: \
         kraken={kraken_name:?}, FASTQ={:?}",
        String::from_utf8_lossy(&f[..name_end])
    );
}

/// Writes a single FASTQ record to a writer.
fn write_fastq_record<W: Write>(writer: &mut W, rec: &impl Record) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(rec.head())?;
    writer.write_all(b"\n")?;
    writer.write_all(rec.seq())?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(rec.qual())?;
    writer.write_all(b"\n")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_report() -> KrakenReport {
        // unclassified(0), root(1), Bacteria(2), E.coli(3), Eukaryota(4), Human(5)
        let lines = [
            " 10.00\t100\t100\tU\t0\tunclassified",
            " 90.00\t900\t5\tR\t1\troot",
            " 60.00\t600\t10\tD\t2\t  Bacteria",
            " 50.00\t500\t500\tS\t3\t    Escherichia coli",
            " 30.00\t300\t10\tD\t4\t  Eukaryota",
            " 20.00\t200\t200\tS\t5\t    Homo sapiens",
        ]
        .join("\n");
        KrakenReport::from_reader(lines.as_bytes()).unwrap()
    }

    #[test]
    fn test_build_taxon_set_exact() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[3], false, false).unwrap();
        assert_eq!(set, HashSet::from([3]));
        assert_eq!(expected, 500);
    }

    #[test]
    fn test_build_taxon_set_with_descendants() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[2], true, false).unwrap();
        assert_eq!(set, HashSet::from([2, 3]));
        assert_eq!(expected, 600);
    }

    #[test]
    fn test_build_taxon_set_with_descendants_root() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[1], true, false).unwrap();
        assert_eq!(set, HashSet::from([1, 2, 3, 4, 5]));
        assert_eq!(expected, 900);
    }

    #[test]
    fn test_build_taxon_set_unknown_taxon() {
        let report = make_report();
        let result = build_taxon_set_and_expected_count(&report, &[99999], false, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_taxon_set_include_unclassified() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[3], false, true).unwrap();
        assert_eq!(set, HashSet::from([0, 3]));
        assert_eq!(expected, 600);
    }

    #[test]
    fn test_build_taxon_set_only_unclassified() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[], false, true).unwrap();
        assert_eq!(set, HashSet::from([0]));
        assert_eq!(expected, 100);
    }

    #[test]
    fn test_expected_count_with_descendants() {
        let report = make_report();
        let (_, expected) = build_taxon_set_and_expected_count(&report, &[2], true, false).unwrap();
        assert_eq!(expected, 600);
    }

    #[test]
    fn test_expected_count_without_descendants() {
        let report = make_report();
        let (_, expected) =
            build_taxon_set_and_expected_count(&report, &[2], false, false).unwrap();
        assert_eq!(expected, 10);
    }

    #[test]
    fn test_expected_count_with_unclassified() {
        let report = make_report();
        let (_, expected) = build_taxon_set_and_expected_count(&report, &[3], false, true).unwrap();
        assert_eq!(expected, 600);
    }

    #[test]
    fn test_validate_read_name_match() {
        assert!(validate_read_name("read1", b"read1", 1).is_ok());
    }

    #[test]
    fn test_validate_read_name_mismatch() {
        assert!(validate_read_name("read1", b"read2", 1).is_err());
    }

    #[test]
    fn test_validate_read_name_strip_suffix_1() {
        assert!(validate_read_name("read1", b"read1/1", 1).is_ok());
    }

    #[test]
    fn test_validate_read_name_strip_suffix_2() {
        assert!(validate_read_name("read1", b"read1/2", 1).is_ok());
    }

    #[test]
    fn test_validate_read_name_with_comment() {
        assert!(validate_read_name("read1", b"read1 length=150", 1).is_ok());
    }

    #[test]
    fn test_validate_read_name_suffix_and_comment() {
        assert!(validate_read_name("read1", b"read1/1 length=150", 1).is_ok());
    }

    #[test]
    fn test_validate_args_mismatched_counts() {
        let filter = Filter {
            kraken_report: PathBuf::from("r.txt"),
            kraken_output: PathBuf::from("k.txt"),
            input: vec![PathBuf::from("a.fq"), PathBuf::from("b.fq")],
            output: vec![PathBuf::from("c.fq")],
            taxon_ids: vec![1],
            include_descendants: false,
            include_unclassified: false,
            threads: 4,
            compression_level: 6,
        };
        assert!(filter.validate_args().is_err());
    }

    #[test]
    fn test_validate_args_no_taxa_or_unclassified() {
        let filter = Filter {
            kraken_report: PathBuf::from("r.txt"),
            kraken_output: PathBuf::from("k.txt"),
            input: vec![PathBuf::from("a.fq")],
            output: vec![PathBuf::from("b.fq")],
            taxon_ids: vec![],
            include_descendants: false,
            include_unclassified: false,
            threads: 4,
            compression_level: 6,
        };
        assert!(filter.validate_args().is_err());
    }
}
