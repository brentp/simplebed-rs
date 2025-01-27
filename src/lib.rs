//! A Rust library for reading BED files, supporting BGZF compression and performance.
//!
//! This library provides `BedReader` and `BedWriter` structs to efficiently parse and write BED files,
//! handling both plain text and BGZF compressed formats using `rust-htslib`.
//! It represents BED records as `BedRecord` structs, providing access to
//! chrom, start, end, name, score, and optional fields with type-safe access.
//!
//! # Example
//!
//! ```
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     use simplebed::{BedRecord,BedValue,BedReader,BedWriter};
//!     // Reading a BED file
//!     let bed_path = std::path::Path::new("tests/test.bed");
//!     let mut bed_reader = BedReader::new(bed_path)?;
//!
//!     // Writing to a BGZF compressed file
//!     let output_path = std::path::Path::new("output.bed.gz"); // .gz extension enables BGZF compression
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
//! * **BGZF support:** Uses "noodles" to read and write BGZF compressed BED files efficiently.
//! * **Type-safe optional fields:** Represents optional columns as a vector of "BedValue" enum,
//!   supporting String, Integer, and Float types.
//! * **Performance-oriented:** Uses buffered reading and efficient string parsing for speed.

use flate2::read::GzDecoder;
use noodles::bgzf;

use noodles::tabix;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use noodles::bgzf::VirtualPosition;
use noodles::core::{Position, Region};
use noodles::csi::BinningIndex;
use noodles::csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use std::ops::Bound;

pub mod value;
pub use value::BedValue;

pub mod writer;
pub use writer::BedWriter;

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

#[derive(Debug, thiserror::Error)]
pub enum BedError {
    #[error("Invalid BED format: {0}")]
    InvalidFormat(String),
    #[error("Parse error: {0}")]
    ParseError(String),
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

/// Represents the type of index available for the BED file
#[derive(Debug)]
enum BedIndex {
    Tabix(tabix::Index),
    Csi(csi::Index),
    None,
}

impl BedIndex {
    fn query(&self, tid: usize, region: &Region) -> io::Result<Vec<Chunk>> {
        match self {
            BedIndex::Tabix(index) => index.query(tid, region.interval()),
            BedIndex::Csi(index) => index.query(tid, region.interval()),
            _ => unimplemented!(),
        }
    }
}

/// A reader for BED files, supporting plain text, BGZF compression, and tabix indexing.
pub struct BedReader {
    reader: Box<dyn BufRead>,
    index: BedIndex,
}

impl BedReader {
    /// Creates a new `BedReader` from a file path.
    /// Automatically detects BGZF compression based on file extension (.gz or .bgz).
    pub fn new<P: AsRef<Path>>(path: P) -> Result<BedReader, Box<dyn std::error::Error>> {
        let mut file = BufReader::new(File::open(path.as_ref())?);
        let compression = detect_compression(&mut file)?;

        let reader = match compression {
            Compression::GZ => {
                let rdr: Box<dyn BufRead> = Box::new(BufReader::new(GzDecoder::new(file)));
                rdr
            }
            Compression::BGZF => {
                let rdr: Box<dyn BufRead> = Box::new(bgzf::Reader::new(file));
                rdr
            }
            _ => Box::new(file),
        };

        BedReader::from_reader(reader, path)
    }

    pub fn from_reader<P: AsRef<Path>>(
        reader: Box<dyn BufRead>,
        path: P,
    ) -> Result<BedReader, Box<dyn std::error::Error>> {
        let path = path.as_ref();
        let index = if let Ok(csi_file) = File::open(format!("{}.csi", path.display())) {
            let csi_reader = BufReader::new(Box::new(csi_file));
            let mut csi = csi::io::Reader::new(csi_reader);
            BedIndex::Csi(csi.read_index()?)
        } else if let Ok(tbx_file) = File::open(format!("{}.tbi", path.display())) {
            let tbx_reader = BufReader::new(Box::new(tbx_file));
            let mut tbx = tabix::io::Reader::new(tbx_reader);
            BedIndex::Tabix(tbx.read_index()?)
        } else {
            BedIndex::None
        };

        Ok(BedReader { reader, index })
    }

    /// Reads a single `BedRecord` from the reader.
    pub fn read_record(&mut self) -> Result<Option<BedRecord>, BedError> {
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
            score = fields[4].parse::<f64>().ok();
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

    /// Returns an iterator over the records in the BED file.
    pub fn records(&mut self) -> Records<'_> {
        Records { reader: self }
    }

    /// Query the BED file for records overlapping the given region.
    ///
    /// # Arguments
    ///
    /// * `tid` - Target ID (reference sequence index)
    /// * `chrom` - Chromosome name
    /// * `start` - Start position (0-based)
    /// * `stop` - Stop position (exclusive)
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use simplebed::BedReader;
    /// use std::path::Path;
    ///
    /// let bed_path = Path::new("test.bed.gz");
    /// let mut reader = BedReader::new(bed_path)?;
    ///
    /// // Query chromosome 1 starting at position 1000 and ending at position 2000
    /// for record in reader.query(0, "chr1", 1000, 2000)? {
    ///     let record = record?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query(
        &mut self,
        tid: usize,
        chrom: &str,
        start: usize,
        stop: usize,
    ) -> Result<QueryRecords<'_>, BedError> {
        // Check if we have an index
        if matches!(self.index, BedIndex::None) {
            return Err(BedError::InvalidFormat(
                "No index found. File must be indexed with tabix or CSI.".to_string(),
            ));
        }

        // Create a region query
        let start = Position::try_from(start).map_err(|e| {
            BedError::ParseError(format!("Invalid start coordinate: {}. Error: {}", start, e))
        })?;
        let stop = Position::try_from(stop).map_err(|e| {
            BedError::ParseError(format!("Invalid stop coordinate: {}. Error: {}", stop, e))
        })?;
        let region = Region::new(chrom.to_string(), start..=stop);

        // Get chunks that overlap this region
        let chunks = self
            .index
            .query(tid, &region)
            .map_err(|e| BedError::IoError(e))?;

        Ok(QueryRecords {
            reader: self,
            chunks,
            region: region.clone(),
        })
    }
}

/// An iterator over the records in a BED file.
pub struct Records<'a> {
    reader: &'a mut BedReader,
}

impl<'a> Iterator for Records<'a> {
    type Item = Result<BedRecord, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.reader.read_record() {
                Ok(Some(record)) => return Some(Ok(record)),
                Ok(None) => return None, // EOF
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

/// An iterator over query results from a BED file.
pub struct QueryRecords<'a> {
    reader: &'a mut BedReader,
    chunks: Vec<Chunk>,
    region: Region,
}

impl<'a> Iterator for QueryRecords<'a> {
    type Item = Result<BedRecord, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.reader.read_record() {
                Ok(Some(record)) => {
                    // Check if record overlaps query region
                    if record.chrom() == self.region.name()
                        && record.end()
                            > match self.region.start() {
                                Bound::Included(pos) => pos.get() as u64,
                                Bound::Excluded(pos) => pos.get() as u64,
                                Bound::Unbounded => 0,
                            }
                        && record.start()
                            < match self.region.end() {
                                Bound::Included(pos) => pos.get() as u64,
                                Bound::Excluded(pos) => pos.get() as u64,
                                Bound::Unbounded => u64::MAX,
                            }
                    {
                        return Some(Ok(record));
                    }
                    // Keep reading until we find an overlapping record
                    continue;
                }
                Ok(None) => return None,
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

#[derive(Debug, PartialEq)]
enum Compression {
    BGZF,
    RAZF,
    GZ,
    None,
}

fn detect_compression<R: BufRead>(
    reader: &mut R,
) -> Result<Compression, Box<dyn std::error::Error>> {
    let buf = reader.fill_buf()?;
    let mut dec_buf = vec![0u8; buf.len()];

    let is_gzipped = &buf[0..2] == b"\x1f\x8b";

    if is_gzipped && buf[3] & 4 != 0 && buf.len() >= 18 {
        let c = match &buf[12..16] {
            b"BC\x02\x00" => Compression::BGZF,
            b"RAZF" => Compression::RAZF,
            _ => Compression::GZ,
        };
        let mut gz = GzDecoder::new(buf);
        match gz.read_exact(&mut dec_buf) {
            Ok(_) => {}
            Err(e) => {
                if e.kind() != io::ErrorKind::UnexpectedEof {
                    return Err(e.into());
                }
            }
        }
        Ok(c)
    } else if is_gzipped {
        Ok(Compression::GZ)
    } else {
        Ok(Compression::None)
    }
}

#[cfg(test)]
mod tests {
    use crate::{BedReader, BedRecord, BedValue, BufRead};
    use std::error::Error;
    use std::path::Path;

    #[test]
    fn test_read_bed_file() -> Result<(), Box<dyn Error>> {
        let bed_path = Path::new("tests/test.bed");
        let mut bed_reader = BedReader::new(bed_path)?;

        let record1 = bed_reader.read_record()?.unwrap();
        assert_eq!(record1.chrom(), "chr1");
        assert_eq!(record1.start(), 22);
        assert_eq!(record1.end(), 44);
        assert_eq!(record1.name(), Some("TTN"));
        assert_eq!(record1.score(), Some(34.5));
        assert_eq!(
            record1.other_fields(),
            &[
                BedValue::String("other_field".to_string()),
                BedValue::Integer(23),
                BedValue::Float(65.4)
            ]
        );

        let record2 = bed_reader.read_record()?.unwrap();
        assert_eq!(record2.chrom(), "chr1");
        assert_eq!(record2.start(), 33);
        assert_eq!(record2.end(), 54);
        assert_eq!(record2.name(), Some("TTN"));
        assert_eq!(record2.score(), None);
        assert_eq!(record2.other_fields(), &[]);

        let record3 = bed_reader.read_record()?.unwrap();
        assert_eq!(record3.chrom(), "chr1");
        assert_eq!(record3.start(), 53);
        assert_eq!(record3.end(), 94);
        assert_eq!(record3.name(), None);
        assert_eq!(record3.score(), None);
        assert_eq!(record3.other_fields(), &[]);

        let record_none = bed_reader.read_record()?;
        assert!(record_none.is_none());

        Ok(())
    }
}
