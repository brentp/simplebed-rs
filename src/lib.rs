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
//!     let mut bed_reader = BedReader::from_path(bed_path)?;
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

use flate2::bufread::GzDecoder;
use noodles::bgzf;

//use noodles::bgzf::io::{BufRead, Read, Seek};
use noodles::core::{Position, Region};
use noodles::csi::BinningIndex;
use noodles::csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles::tabix;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek};
use std::marker::PhantomData;
use std::path::Path;
pub mod value;
pub use value::BedValue;

pub mod record;
pub use record::BedRecord;

pub mod writer;
pub use writer::BedWriter;

#[derive(Debug, thiserror::Error)]
pub enum BedError {
    #[error("Invalid BED format: {0}")]
    InvalidFormat(String),
    #[error("Parse error: {0}")]
    ParseError(String),
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("No chromosome order found.")]
    NoChromosomeOrder,
}

/// Represents the type of index available for the BED file
#[derive(Debug)]
enum BedIndex {
    Tabix(tabix::Index),
    Csi(csi::Index),
    None,
}

impl BedIndex {
    fn query(&self, tid: usize, region: &Region) -> Option<io::Result<Vec<Chunk>>> {
        match self {
            BedIndex::Tabix(index) => Some(index.query(tid, region.interval())),
            BedIndex::Csi(index) => Some(index.query(tid, region.interval())),
            BedIndex::None => None,
        }
    }

    /// Return the csi header which contains the reference sequence names
    fn header(&self) -> Option<&csi::binning_index::index::Header> {
        match self {
            BedIndex::Tabix(index) => index.header(),
            BedIndex::Csi(index) => index.header(),
            BedIndex::None => None,
        }
    }
}

/// Represents different types of readers for BED files
pub enum BedReaderType<R: Read + Seek> {
    Plain(BufReader<R>),
    Gzip(BufReader<GzDecoder<BufReader<R>>>), // Wrap GzDecoder in BufReader
    Bgzf(bgzf::Reader<BufReader<R>>),
}

impl std::fmt::Debug for BedReaderType<File> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BedReaderType::Plain(_reader) => write!(f, "Plain(reader)"),
            BedReaderType::Gzip(_reader) => write!(f, "Gzip(reader)"),
            BedReaderType::Bgzf(_reader) => write!(f, "Bgzf(reader)"),
        }
    }
}

impl BedReaderType<File> {
    fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        match self {
            BedReaderType::Plain(reader) => reader.read_line(buf),
            BedReaderType::Gzip(reader) => reader.read_line(buf), // Now works with BufReader
            BedReaderType::Bgzf(reader) => reader.read_line(buf),
        }
    }
}

/// A reader for BED files, supporting plain text, BGZF compression, and tabix indexing.
#[derive(Debug)]
pub struct BedReader {
    reader: BedReaderType<File>,
    index: BedIndex,
    current_region: Option<Region>,
    chromosome_order: Option<HashMap<String, usize>>,
}

impl BedReader {
    /// Creates a new `BedReader` from a BufReader.
    pub fn new<P: AsRef<Path>>(
        reader: File,
        path: P,
    ) -> Result<BedReader, Box<dyn std::error::Error>> {
        let mut reader = BufReader::new(reader);
        let compression = detect_compression(&mut reader)?;

        let reader = match compression {
            Compression::GZ => BedReaderType::Gzip(BufReader::new(GzDecoder::new(reader))),
            Compression::BGZF => BedReaderType::Bgzf(bgzf::Reader::new(reader)),
            _ => BedReaderType::Plain(reader),
        };

        BedReader::from_reader(reader, path)
    }

    // Add a new constructor for path-based creation
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<BedReader, Box<dyn std::error::Error>> {
        let reader = File::open(path.as_ref())?;
        Self::new(reader, path)
    }

    pub fn from_reader<P: AsRef<Path>>(
        reader: BedReaderType<File>,
        path: P,
    ) -> Result<BedReader, Box<dyn std::error::Error>> {
        let path = path.as_ref();
        let index = if let Ok(csi_file) = File::open(format!("{}.csi", path.display())) {
            let csi_reader = BufReader::new(csi_file);
            let mut csi = csi::io::Reader::new(csi_reader);
            BedIndex::Csi(csi.read_index()?)
        } else if let Ok(tbx_file) = File::open(format!("{}.tbi", path.display())) {
            let tbx_reader = BufReader::new(tbx_file);
            let mut tbx = tabix::io::Reader::new(tbx_reader);
            BedIndex::Tabix(tbx.read_index()?)
        } else {
            BedIndex::None
        };

        Ok(BedReader {
            reader,
            index,
            current_region: None,
            chromosome_order: None,
        })
    }

    /// Set the chromosome order for the reader. Used for query without index.
    pub fn set_chromosome_order(&mut self, chromosome_order: HashMap<String, usize>) {
        self.chromosome_order = Some(chromosome_order);
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
        Self::parse_line(line)
    }

    fn parse_line(line: &str) -> Result<Option<BedRecord>, BedError> {
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
    pub fn records(&mut self) -> Records {
        Records { reader: self }
    }
}

impl BedReader {
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
    /// let mut reader = BedReader::from_path(bed_path)?;
    ///
    /// // Query chromosome 1 starting at position 1000 and ending at position 2000
    /// for record in reader.query(0, "chr1", 1000, 2000)? {
    ///     let record = record?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r>(
        &'r mut self,
        tid: usize,
        chrom: &str,
        start: usize,
        stop: usize,
    ) -> Result<Box<dyn Iterator<Item = Result<BedRecord, BedError>> + 'r>, BedError> {
        /*
        // Check if we have an index, otherwise fallback to linear scan
        if !matches!(self.index, BedIndex::None) {
            return Err(BedError::InvalidFormat(
                "Index found but not used in this implementation.".to_string(),
            ));
        }
        */

        let target_chrom_order = if let Some(chrom_order) = &self.chromosome_order {
            match chrom_order.get(chrom) {
                Some(order) => Ok(*order),
                None => {
                    return Err(BedError::InvalidFormat(format!(
                        "Chromosome {} not found in chromosome order.",
                        chrom
                    )))
                }
            }
        } else {
            Err(BedError::NoChromosomeOrder)
        };

        // Create a region query
        let start_pos = Position::try_from(start).map_err(|e| {
            BedError::ParseError(format!("Invalid start coordinate: {}. Error: {}", start, e))
        })?;
        let stop_pos = Position::try_from(stop).map_err(|e| {
            BedError::ParseError(format!("Invalid stop coordinate: {}. Error: {}", stop, e))
        })?;
        self.current_region = Some(Region::new(chrom.to_string(), start_pos..=stop_pos));

        // Get chunks that overlap this region
        let qchunks = self.index.query(tid, self.current_region.as_ref().unwrap());

        let iter: Box<dyn Iterator<Item = Result<BedRecord, BedError>>> =
            match (qchunks, &mut self.reader) {
                (None, BedReaderType::Plain(reader)) => Box::new(LinearScanIterator::new(
                    reader,
                    target_chrom_order?,
                    start as u64,
                    stop as u64,
                    self.chromosome_order
                        .as_ref()
                        .ok_or(BedError::NoChromosomeOrder)?,
                )),
                (None, BedReaderType::Gzip(reader)) => Box::new(LinearScanIterator::new(
                    reader,
                    target_chrom_order?,
                    start as u64,
                    stop as u64,
                    self.chromosome_order
                        .as_ref()
                        .ok_or(BedError::NoChromosomeOrder)?,
                )),
                (None, BedReaderType::Bgzf(reader)) => Box::new(LinearScanIterator::new(
                    reader,
                    target_chrom_order?,
                    start as u64,
                    stop as u64,
                    self.chromosome_order
                        .as_ref()
                        .ok_or(BedError::NoChromosomeOrder)?,
                )),
                (Some(Ok(chunks)), BedReaderType::Bgzf(ref mut reader)) => {
                    let q = csi::io::Query::new(reader, chunks);
                    let header = self.index.header().expect("No header found");
                    let q = q
                        .indexed_records(header)
                        .filter_by_region(self.current_region.as_ref().unwrap());

                    // Create an iterator that converts Record to BedRecord
                    let iter = BedRecordIterator {
                        inner: q,
                        _phantom: PhantomData,
                        //region: &self.current_region.as_ref().unwrap(),
                    };
                    Box::new(iter)
                }
                _ => unimplemented!(),
            };

        Ok(iter)
    }
}

struct LinearScanIterator<'a, R> {
    reader: R,
    target_chrom_order: usize,
    start_pos: u64,
    stop_pos: u64,
    chromosome_order: &'a HashMap<String, usize>,
    buffer: String,
    finished: bool,
}

impl<'a, R: BufRead> LinearScanIterator<'a, R> {
    fn new(
        reader: R,
        target_chrom_order: usize,
        start_pos: u64,
        stop_pos: u64,
        chromosome_order: &'a HashMap<String, usize>,
    ) -> Self {
        LinearScanIterator {
            reader,
            target_chrom_order,
            start_pos,
            stop_pos,
            chromosome_order,
            buffer: String::new(),
            finished: false,
        }
    }
}

impl<R: BufRead> Iterator for LinearScanIterator<'_, R> {
    type Item = Result<BedRecord, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        loop {
            self.buffer.clear();
            let bytes_read = match self.reader.read_line(&mut self.buffer) {
                Ok(bytes) => bytes,
                Err(e) => return Some(Err(BedError::IoError(e))),
            };

            if bytes_read == 0 {
                self.finished = true;
                return None; // EOF
            }

            let line = self.buffer.trim_end(); // Remove trailing newline

            if line.starts_with('#') || line.is_empty() {
                continue; // Skip comments and empty lines
            }

            let record = match BedReader::parse_line(line) {
                Ok(Some(rec)) => rec,
                Ok(None) => continue,
                Err(e) => return Some(Err(e)),
            };

            let current_chrom_order = match self.chromosome_order.get(record.chrom()) {
                Some(order) => *order,
                None => {
                    return Some(Err(BedError::InvalidFormat(format!(
                        "Chromosome {} not found in chromosome order.",
                        record.chrom()
                    ))))
                }
            };

            // If we've moved past the target chromosome, we're done
            if current_chrom_order > self.target_chrom_order {
                self.finished = true;
                return None;
            }

            // Skip until we reach the target chromosome
            if current_chrom_order < self.target_chrom_order {
                continue;
            }

            // Now we're on the target chromosome, check positions
            if record.start() >= self.stop_pos {
                self.finished = true;
                return None; // Past the target interval
            }

            if record.end() > self.start_pos {
                return Some(Ok(record)); // Record overlaps with the target interval
            }
        }
    }
}

/// An iterator over the records in a BED file.
pub struct Records<'a> {
    reader: &'a mut BedReader,
}

impl Iterator for Records<'_> {
    type Item = Result<BedRecord, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None, // EOF
            Err(e) => Some(Err(e)),
        }
    }
}

#[derive(Debug, PartialEq)]
#[allow(clippy::upper_case_acronyms)]
enum Compression {
    BGZF,
    RAZF,
    GZ,
    None,
}

fn detect_compression<R: std::io::BufRead>(
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
        if matches!(c, Compression::BGZF) {
            let mut gz = GzDecoder::new(buf);
            match gz.read_exact(&mut dec_buf) {
                Ok(_) => {}
                Err(e) => {
                    if e.kind() != io::ErrorKind::UnexpectedEof {
                        return Err(e.into());
                    }
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

// Add this new struct to convert Record to BedRecord
struct BedRecordIterator<I>
where
    I: Iterator<Item = std::io::Result<Record>>,
{
    inner: I,
    _phantom: PhantomData<I>,
}

use noodles::csi::io::indexed_records::Record;

impl<I> Iterator for BedRecordIterator<I>
where
    I: Iterator<Item = std::io::Result<Record>>,
{
    type Item = Result<BedRecord, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|result| {
            result
                .map(|record| {
                    let buf = record.as_ref();
                    BedReader::parse_line(buf)
                        .expect("Failed to parse line")
                        .expect("Failed to parse line")
                })
                .map_err(BedError::IoError)
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::{BedReader, BedRecord, BedValue, HashMap};
    use std::error::Error;
    use std::path::Path;

    #[test]
    fn test_read_bed_file() -> Result<(), Box<dyn Error>> {
        let bed_path = Path::new("tests/test.bed");
        let mut bed_reader = BedReader::from_path(bed_path)?;

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

    #[test]
    fn test_query_bed_file() -> Result<(), Box<dyn Error>> {
        let bed_path = Path::new("tests/compr.bed.gz");
        let mut bed_reader = BedReader::from_path(bed_path)?;
        bed_reader.set_chromosome_order(HashMap::from([
            ("chr1".to_string(), 0),
            ("chr2".to_string(), 1),
        ]));

        let records: Vec<BedRecord> = bed_reader
            .query(0, "chr1", 22, 34)?
            .collect::<Result<Vec<_>, _>>()?;

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom(), "chr1");
        assert_eq!(records[0].start(), 22);
        assert_eq!(records[0].end(), 44);
        assert_eq!(records[1].chrom(), "chr1");
        assert_eq!(records[1].start(), 33);
        assert_eq!(records[1].end(), 54);

        Ok(())
    }

    #[test]
    fn test_query_first_interval() -> Result<(), Box<dyn Error>> {
        let bed_path = Path::new("tests/compr.bed.gz");
        let mut bed_reader = BedReader::from_path(bed_path)?;
        bed_reader.set_chromosome_order(HashMap::from([
            ("chr1".to_string(), 0),
            ("chr2".to_string(), 1),
        ]));

        let records: Vec<BedRecord> = bed_reader
            .query(0, "chr1", 1, 3)?
            .collect::<Result<Vec<_>, _>>()?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chrom(), "chr1");
        assert_eq!(records[0].start(), 1);
        assert_eq!(records[0].end(), 10);
        Ok(())
    }

    #[test]
    fn test_get_big_interval_query() -> Result<(), Box<dyn Error>> {
        let bed_path = Path::new("tests/compr.bed.gz");
        eprintln!("starting");
        let mut bed_reader = BedReader::from_path(bed_path)?;
        eprintln!("{:?}", bed_reader);
        bed_reader.set_chromosome_order(HashMap::from([
            ("chr1".to_string(), 0),
            ("chr2".to_string(), 1),
        ]));
        eprintln!("{:?}", bed_reader.chromosome_order);

        let records: Vec<BedRecord> = bed_reader
            .query(0, "chr1", 999998, 999999)?
            .collect::<Result<Vec<_>, _>>()?;
        assert_eq!(records.len(), 1);
        Ok(())
    }
}
