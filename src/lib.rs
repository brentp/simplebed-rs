//! A Rust library for reading BED files, supporting BGZF compression and performance.
//!
//! This library provides `BedReader` and `BedWriter` structs to efficiently parse and write BED files,
//! handling both plain text and BGZF compressed formats using `rust-htslib`.
//! It represents BED records as `BedRecord` structs, providing access to
//! chrom, start, end, name, score, and optional fields with type-safe access.
//!
//! ## Usage
//!
//! ```rust
//! use simplebed::{BedReader, BedWriter, BedRecord, BedValue};
//! use std::path::Path;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Reading a BED file
//!     let bed_path = Path::new("tests/test.bed");
//!     let mut bed_reader = BedReader::new(bed_path)?;
//!
//!     // Writing to a BGZF compressed file
//!     let output_path = Path::new("output.bed.gz"); // .gz extension enables BGZF compression
//!     let mut bed_writer = BedWriter::new(output_path)?;
//!
//!     // Process records and write them to the compressed file
//!     for record_result in bed_reader.records() {
//!         let record = record_result?;
//!         println!(
//!             "Processing: {}, Start: {}, End: {}",
//!             record.chrom(),
//!             record.start(),
//!             record.end()
//!         );
//!         bed_writer.write_record(&record)?;
//!     }
//!
//!     // Ensure all data is written
//!     bed_writer.flush()?;
//!
//!     // Create and write a new record
//!     let new_record = BedRecord::new(
//!         "chr1".to_string(),
//!         1000,
//!         2000,
//!         Some("gene_name".to_string()),
//!         Some(100.0),
//!         vec![BedValue::String("additional_field".to_string())],
//!     );
//!     bed_writer.write_record(&new_record)?;
//!
//!     Ok(())
//! }
//! ```
//!
//! ## Features
//!
//! * **BED format parsing and writing:** Handles standard BED format with chrom, start, end,
//!   name, score, and optional columns.
//! * **BGZF support:** Uses `rust-htslib` to read and write BGZF compressed BED files efficiently.
//! * **Type-safe optional fields:** Represents optional columns as a vector of `BedValue` enum,
//!   supporting String, Integer, and Float types.
//! * **Performance-oriented:** Uses buffered reading and efficient string parsing for speed.
//! * **Error handling:** Provides a custom `BedError` enum for robust error management.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use rust_htslib::bgzf;

/// Represents the possible data types for optional BED fields.
#[derive(Debug, Clone, PartialEq)]
pub enum BedValue {
    String(String),
    Int(i64),
    Float(f64),
}

impl std::fmt::Display for BedValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BedValue::String(s) => write!(f, "{}", s),
            BedValue::Int(i) => write!(f, "{}", i),
            BedValue::Float(fl) => write!(f, "{:.4}", fl),
        }
    }
}

impl BedValue {
    /// Attempts to parse a string slice into a `BedValue`.
    /// Tries parsing as Int, then Float, then falls back to String.
    fn parse(s: &str) -> Self {
        if let Ok(int_val) = s.parse::<i64>() {
            BedValue::Int(int_val)
        } else if let Ok(float_val) = s.parse::<f64>() {
            BedValue::Float(float_val)
        } else {
            BedValue::String(s.to_string())
        }
    }
}

/// Represents a BED record.
#[derive(Debug, Clone, PartialEq)]
pub struct BedRecord {
    chrom: String,
    start: u64,
    end: u64,
    name: Option<String>,
    score: Option<f64>,
    other_fields: Vec<BedValue>,
}

impl BedRecord {
    /// Creates a new `BedRecord`.
    pub fn new(
        chrom: String,
        start: u64,
        end: u64,
        name: Option<String>,
        score: Option<f64>,
        other_fields: Vec<BedValue>,
    ) -> Self {
        BedRecord {
            chrom,
            start,
            end,
            name,
            score,
            other_fields,
        }
    }

    /// Returns the chromosome.
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// Returns the start coordinate (0-based).
    pub fn start(&self) -> u64 {
        self.start
    }

    /// Returns the end coordinate (exclusive).
    pub fn end(&self) -> u64 {
        self.end
    }

    /// Returns the name, if present.
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Returns the score, if present.
    pub fn score(&self) -> Option<f64> {
        self.score
    }

    /// Returns a slice of the other fields.
    pub fn other_fields(&self) -> &[BedValue] {
        &self.other_fields
    }
}

/// Errors that can occur during BED file reading.
#[derive(Debug, thiserror::Error)]
pub enum BedError {
    #[error("IO error: {0}")]
    IoError(#[from] io::Error),
    #[error("BGZF error: {0}")]
    BgzfError(#[from] rust_htslib::errors::Error),
    #[error("Invalid BED format: {0}")]
    InvalidFormat(String),
    #[error("Parse error: {0}")]
    ParseError(String),
}

/// A reader for BED files, supporting plain text and BGZF compression.
pub struct BedReader<R: BufRead> {
    reader: R,
}

impl BedReader<BufReader<Box<dyn Read>>> {
    /// Creates a new `BedReader` from a file path.
    /// Automatically detects BGZF compression based on file extension (.gz or .bgz).
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, BedError> {
        let path = path.as_ref();
        let reader: BufReader<Box<dyn Read>> = if path.extension().and_then(|s| s.to_str())
            == Some("gz")
            || path.extension().and_then(|s| s.to_str()) == Some("bgz")
        {
            let bgzf_reader = bgzf::Reader::from_path(path)?;
            BufReader::new(Box::new(bgzf_reader))
        } else {
            let file = File::open(path)?;
            BufReader::new(Box::new(file))
        };
        Ok(BedReader { reader })
    }
}

impl<R: BufRead> BedReader<R> {
    /// Creates a new `BedReader` from a `BufRead` source.
    pub fn from_buf_read(reader: R) -> Self {
        BedReader { reader }
    }

    /// Returns an iterator over `BedRecord`s in the BED file.
    pub fn records(&mut self) -> Records<'_, R> {
        Records { reader: self }
    }

    /// Reads a single `BedRecord` from the reader.
    fn read_record(&mut self) -> Result<Option<BedRecord>, BedError> {
        let mut line = String::new();
        let bytes_read = self.reader.read_line(&mut line)?;
        if bytes_read == 0 {
            return Ok(None); // EOF
        }

        let line = line.trim_end(); // Remove trailing newline

        if line.starts_with('#') || line.is_empty() {
            return self.read_record(); // Skip comments and empty lines
        }

        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 3 {
            return Err(BedError::InvalidFormat(
                "BED format requires at least 3 columns (chrom, start, end)".to_string(),
            ));
        }

        let chrom = fields[0].to_string();
        let start = fields[1].parse::<u64>().map_err(|e| {
            BedError::ParseError(format!(
                "Invalid start coordinate: {}. Error: {}",
                fields[1], e
            ))
        })?;
        let end = fields[2].parse::<u64>().map_err(|e| {
            BedError::ParseError(format!(
                "Invalid end coordinate: {}. Error: {}",
                fields[2], e
            ))
        })?;

        let mut name = None;
        if fields.len() > 3 && !fields[3].is_empty() {
            name = Some(fields[3].to_string());
        }

        let mut score = None;
        if fields.len() > 4 && !fields[4].is_empty() {
            score = fields[4].parse::<f64>().ok(); // Ignore parse errors for score, as it's optional and can be "."
        }

        let mut other_fields = Vec::new();
        for field in fields.iter().skip(5) {
            other_fields.push(BedValue::parse(field));
        }

        Ok(Some(BedRecord::new(
            chrom,
            start,
            end,
            name,
            score,
            other_fields,
        )))
    }
}

/// An iterator over `BedRecord`s.
pub struct Records<'reader, R: BufRead> {
    reader: &'reader mut BedReader<R>,
}

impl<'reader, R: BufRead> Iterator for Records<'reader, R> {
    type Item = Result<BedRecord, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None, // EOF
            Err(e) => Some(Err(e)),
        }
    }
}

/// A writer for BED files, supporting plain text and BGZF compression.
pub struct BedWriter<W: Write> {
    writer: W,
}

impl BedWriter<Box<dyn Write>> {
    /// Creates a new `BedWriter` from a file path.
    /// Automatically uses BGZF compression if the path ends with .gz or .bgz.
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, BedError> {
        let path = path.as_ref();
        let writer: Box<dyn Write> = if path.extension().and_then(|s| s.to_str()) == Some("gz")
            || path.extension().and_then(|s| s.to_str()) == Some("bgz")
        {
            Box::new(bgzf::Writer::from_path(path)?)
        } else {
            Box::new(File::create(path)?)
        };
        Ok(BedWriter { writer })
    }
}

impl<W: Write> BedWriter<W> {
    /// Creates a new `BedWriter` from any type implementing `Write`.
    pub fn from_writer(writer: W) -> Self {
        BedWriter { writer }
    }

    /// Writes a single BED record.
    pub fn write_record(&mut self, record: &BedRecord) -> Result<(), BedError> {
        // Write required fields
        write!(
            self.writer,
            "{}\t{}\t{}",
            record.chrom(),
            record.start(),
            record.end()
        )?;

        // Write optional name
        if let Some(name) = record.name() {
            write!(self.writer, "\t{}", name)?;
        } else if !record.other_fields().is_empty() || record.score().is_some() {
            write!(self.writer, "\t.")?;
        }

        // Write optional score
        if let Some(score) = record.score() {
            write!(self.writer, "\t{}", score)?;
        } else if !record.other_fields().is_empty() {
            write!(self.writer, "\t.")?;
        }

        // Write other fields
        for field in record.other_fields() {
            write!(self.writer, "\t{}", field)?;
        }

        writeln!(self.writer)?;
        Ok(())
    }

    /// Flushes any buffered data to the underlying writer.
    pub fn flush(&mut self) -> Result<(), BedError> {
        self.writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use tempfile::TempDir;

    // Test fixture struct to manage test files
    struct TestBed {
        _temp_dir: TempDir, // Keep TempDir alive while we need the files
        bed_path: PathBuf,
        bgzf_path: PathBuf,
    }

    impl TestBed {
        fn new() -> Result<Self, BedError> {
            let temp_dir = TempDir::new().unwrap();
            let bed_path = temp_dir.path().join("test.bed");
            let bgzf_path = temp_dir.path().join("test.bed.gz");

            // Create test records
            let records = vec![
                BedRecord::new(
                    "chr1".to_string(),
                    100,
                    200,
                    Some("gene1".to_string()),
                    Some(1.5),
                    vec![BedValue::String("field1".to_string()), BedValue::Int(42)],
                ),
                BedRecord::new(
                    "chr2".to_string(),
                    300,
                    400,
                    Some("gene2".to_string()),
                    None,
                    vec![BedValue::String("field2".to_string())],
                ),
                BedRecord::new("chr3".to_string(), 500, 600, None, None, vec![]),
            ];

            // Write plain text BED file
            let mut writer = BedWriter::new(&bed_path)?;
            for record in &records {
                writer.write_record(record)?;
            }
            writer.flush()?;

            // Write BGZF compressed BED file
            let mut bgzf_writer = BedWriter::new(&bgzf_path)?;
            for record in &records {
                bgzf_writer.write_record(record)?;
            }
            bgzf_writer.flush()?;

            Ok(TestBed {
                _temp_dir: temp_dir,
                bed_path,
                bgzf_path,
            })
        }
    }

    #[test]
    fn test_parse_bed_value() {
        assert_eq!(
            BedValue::parse("test"),
            BedValue::String("test".to_string())
        );
        assert_eq!(BedValue::parse("123"), BedValue::Int(123));
        assert_eq!(BedValue::parse("-456"), BedValue::Int(-456));
        assert_eq!(BedValue::parse("1.23"), BedValue::Float(1.23));
        assert_eq!(BedValue::parse("-3.14"), BedValue::Float(-3.14));
    }

    #[test]
    fn test_read_bed_record() -> Result<(), BedError> {
        let test_bed = TestBed::new()?;
        let mut reader = BedReader::new(&test_bed.bed_path)?;

        let record1 = reader.read_record()?.unwrap();
        assert_eq!(record1.chrom(), "chr1");
        assert_eq!(record1.start(), 100);
        assert_eq!(record1.end(), 200);
        assert_eq!(record1.name(), Some("gene1"));
        assert_eq!(record1.score(), Some(1.5));
        assert_eq!(
            record1.other_fields(),
            &[BedValue::String("field1".to_string()), BedValue::Int(42)]
        );

        let record2 = reader.read_record()?.unwrap();
        assert_eq!(record2.chrom(), "chr2");
        assert_eq!(record2.start(), 300);
        assert_eq!(record2.end(), 400);
        assert_eq!(record2.name(), Some("gene2"));
        assert_eq!(record2.score(), None);
        assert_eq!(
            record2.other_fields(),
            &[BedValue::String("field2".to_string())]
        );

        let record3 = reader.read_record()?.unwrap();
        assert_eq!(record3.chrom(), "chr3");
        assert_eq!(record3.start(), 500);
        assert_eq!(record3.end(), 600);
        assert_eq!(record3.name(), None);
        assert_eq!(record3.score(), None);
        assert!(record3.other_fields().is_empty());

        assert!(reader.read_record()?.is_none()); // EOF

        Ok(())
    }

    #[test]
    fn test_read_bgzf_bed() -> Result<(), BedError> {
        let test_bed = TestBed::new()?;
        let mut reader = BedReader::new(&test_bed.bgzf_path)?;

        // Read and verify first record
        let record = reader.read_record()?.unwrap();
        assert_eq!(record.chrom(), "chr1");
        assert_eq!(record.start(), 100);
        assert_eq!(record.end(), 200);
        assert_eq!(record.name(), Some("gene1"));
        assert_eq!(record.score(), Some(1.5));

        // Skip remaining records
        let count = reader.records().count();
        assert_eq!(count, 2); // Two more records

        Ok(())
    }

    #[test]
    fn test_bed_writer() -> Result<(), BedError> {
        let mut output = Vec::new();
        let mut writer = BedWriter::from_writer(&mut output);

        let record = BedRecord::new(
            "chr1".to_string(),
            100,
            200,
            Some("gene1".to_string()),
            Some(1.5),
            vec![BedValue::String("field1".to_string()), BedValue::Int(42)],
        );

        writer.write_record(&record)?;
        writer.flush()?;

        let written = String::from_utf8(output).unwrap();
        assert_eq!(written, "chr1\t100\t200\tgene1\t1.5\tfield1\t42\n");

        Ok(())
    }

    #[test]
    fn test_bed_writer_minimal() -> Result<(), BedError> {
        let mut output = Vec::new();
        let mut writer = BedWriter::from_writer(&mut output);

        let record = BedRecord::new("chr1".to_string(), 100, 200, None, None, vec![]);

        writer.write_record(&record)?;
        writer.flush()?;

        let written = String::from_utf8(output).unwrap();
        assert_eq!(written, "chr1\t100\t200\n");

        Ok(())
    }
}
