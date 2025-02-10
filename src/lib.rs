#![warn(missing_docs)]
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
//!     use std::fs::File;
//!     use std::collections::HashMap;
//!     // Reading a BED file
//!     let bed_path = std::path::Path::new("tests/test.bed");
//!     let mut bed_reader = BedReader::<File>::from_path(bed_path)?;
//!
//!     // the chromosome order makes sure the chromsome order is as expected when querying without an index.
//!     let chromosome_order = HashMap::from([("chr1".to_string(), 0), ("chr2".to_string(), 1)]);
//!     bed_reader.set_chromosome_order(chromosome_order);
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
use noodles::csi;
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

pub mod index;
use index::{BedIndex, Skip, SkipIndex};

/// Errors available in this library.
#[derive(Debug, thiserror::Error)]
pub enum BedError {
    /// Invalid BED format.
    #[error("Invalid BED format: {0}")]
    InvalidFormat(String),
    /// Parse error.
    #[error("Parse error: {0}")]
    ParseError(String),
    /// IO error.
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    /// Ng no chromosome order found when querying.
    #[error("No chromosome order found.")]
    NoChromosomeOrder,
    /// Bed File chromosome order does not match genome file.
    #[error("Bed File chromosome order does not match genome file.")]
    ChromosomeOrderMismatch,
    /// Chromosome not found in chromosome order.
    #[error("Chromosome not found in chromosome order.")]
    ChromosomeNotFound,
}

/// Represents different types of readers for BED files
pub enum BedReaderType<R: Read> {
    /// A plain text reader.
    Plain(BufReader<R>),
    /// A gzip compressed reader.
    Gzip(BufReader<GzDecoder<BufReader<R>>>), // Wrap GzDecoder in BufReader
    /// A bgzf compressed reader.
    Bgzf(bgzf::Reader<BufReader<R>>),
}

impl<R: Read> std::fmt::Debug for BedReaderType<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BedReaderType::Plain(_reader) => write!(f, "Plain(reader)"),
            BedReaderType::Gzip(_reader) => write!(f, "Gzip(reader)"),
            BedReaderType::Bgzf(_reader) => write!(f, "Bgzf(reader)"),
        }
    }
}

impl<R: Read> BedReaderType<R> {
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
pub struct BedReader<R: Read> {
    reader: BedReaderType<R>,
    index: BedIndex,
    current_region: Option<Region>,
    chromosome_order: Option<HashMap<String, usize>>,
    skip_index: Option<SkipIndex>,
    buf: String,
}

impl<R: Read> BedReader<R> {
    /// Set the chromosome order for the reader. Used for query without index.
    pub fn set_chromosome_order(&mut self, chromosome_order: HashMap<String, usize>) {
        self.chromosome_order = Some(chromosome_order);
    }

    /// Reads a single `BedRecord` from the reader.
    pub fn read_record(&mut self) -> Result<Option<BedRecord>, BedError> {
        self.buf.clear();
        let bytes_read = self.reader.read_line(&mut self.buf)?;
        if bytes_read == 0 {
            return Ok(None); // EOF
        }

        let line = self.buf.trim_end(); // Remove trailing newline

        if line.starts_with('#') || line.is_empty() {
            return self.read_record(); // Skip comments and empty lines
        }
        BedRecord::parse_line(line)
    }

    /// Returns an iterator over the records in the BED file.
    pub fn records(&mut self) -> Records<R> {
        Records { reader: self }
    }

    /// Creates a new `BedReader` from a BufReader.
    pub fn new<P: AsRef<Path>>(
        reader: R,
        path: P,
    ) -> Result<BedReader<R>, Box<dyn std::error::Error>> {
        let mut reader = BufReader::new(reader);
        let compression = detect_compression(&mut reader)?;

        let reader = match compression {
            Compression::GZ => BedReaderType::Gzip(BufReader::new(GzDecoder::new(reader))),
            Compression::BGZF => BedReaderType::Bgzf(bgzf::Reader::new(reader)),
            _ => BedReaderType::Plain(reader),
        };

        BedReader::from_reader(reader, path)
    }

    /// Creates a new `BedReader` from a file path.
    pub fn from_path<P: AsRef<Path>>(
        path: P,
    ) -> Result<BedReader<File>, Box<dyn std::error::Error>> {
        let reader = File::open(path.as_ref())?;
        BedReader::<File>::new(reader, path)
    }

    /// Creates a new `BedReader` from a reader.
    pub fn from_reader<P: AsRef<Path>>(
        reader: BedReaderType<R>,
        path: P,
    ) -> Result<BedReader<R>, Box<dyn std::error::Error>> {
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

        let mut br = BedReader {
            reader,
            index,
            current_region: None,
            chromosome_order: None,
            skip_index: None,
            buf: String::new(),
        };
        br.chromosome_order = br.chrom_order_from_index();
        Ok(br)
    }

    fn chrom_order_from_index(&self) -> Option<HashMap<String, usize>> {
        match &self.index {
            BedIndex::Csi(_) | BedIndex::Tabix(_) => {
                let header = self.index.header().expect("No header found");
                let chrom_order = header
                    .reference_sequence_names()
                    .iter()
                    .enumerate()
                    .map(|(i, name)| (name.to_string(), i))
                    .collect();
                Some(chrom_order)
            }
            BedIndex::None => None,
        }
    }
}

const SKIP_CHUNKS: u64 = 1000;

impl<R: Read + Seek> BedReader<R> {
    /// Query the BED file for records overlapping the given region.
    ///
    /// # Arguments
    ///
    /// * `chrom` - Chromosome name
    /// * `start` - Start position (0-based)
    /// * `stop` - Stop position (exclusive)
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use simplebed::BedReader;
    /// use std::{path::Path, fs::File};
    ///
    /// let bed_path = Path::new("test.bed.gz");
    /// let mut reader = BedReader::<File>::from_path(bed_path)?;
    ///
    /// // Query chromosome 1 starting at position 1000 and ending at position 2000
    /// for record in reader.query("chr1", 1000, 2000)? {
    ///     let record = record?;
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r>(
        &'r mut self,
        chrom: &str,
        start: usize,
        stop: usize,
    ) -> Result<Box<dyn Iterator<Item = Result<BedRecord, BedError>> + 'r>, BedError> {
        let tid = if let Some(chrom_order) = &self.chromosome_order {
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

        // Continue with existing query logic
        let start_pos = Position::try_from(start + 1).map_err(|e| {
            BedError::ParseError(format!("Invalid start coordinate: {}. Error: {}", start, e))
        })?;
        let stop_pos = Position::try_from(stop).map_err(|e| {
            BedError::ParseError(format!("Invalid stop coordinate: {}. Error: {}", stop, e))
        })?;

        // If we have a skip index and we're querying a different chromosome,
        // or we are on same chrom but the new start is before the current end,
        // try to seek to the approximate position
        if let Some(skip_index) = &self.skip_index {
            let tid = tid.as_ref().unwrap(); // we know we have tid since we checked it above
            if let Some(seek_pos) = skip_index.find_seek_position(*tid) {
                match &mut self.reader {
                    BedReaderType::Plain(reader) => {
                        reader.seek(std::io::SeekFrom::Start(seek_pos))?;
                    }
                    _ => unreachable!(),
                }
            }
        }
        self.current_region = Some(Region::new(chrom.to_string(), start_pos..=stop_pos));

        // Get chunks that overlap this region
        let qchunks = if let Ok(tid) = tid {
            self.index.query(tid, self.current_region.as_ref().unwrap())
        } else {
            None
        };

        let iter: Box<dyn Iterator<Item = Result<BedRecord, BedError>>> =
            match (qchunks, &mut self.reader) {
                (None, BedReaderType::Plain(reader)) => Box::new(LinearScanIterator::new(
                    reader,
                    tid?,
                    start as u64,
                    stop as u64,
                    self.chromosome_order
                        .as_ref()
                        .ok_or(BedError::NoChromosomeOrder)?,
                )),
                (None, BedReaderType::Gzip(reader)) => Box::new(LinearScanIterator::new(
                    reader,
                    tid?,
                    start as u64,
                    stop as u64,
                    self.chromosome_order
                        .as_ref()
                        .ok_or(BedError::NoChromosomeOrder)?,
                )),
                (None, BedReaderType::Bgzf(reader)) => Box::new(LinearScanIterator::new(
                    reader,
                    tid?,
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

    /// Build a skip index for faster chromosome access
    pub fn build_skip_index(&mut self) -> Result<(), BedError> {
        // Only build if we have chromosome order
        let chromosome_order = self
            .chromosome_order
            .as_ref()
            .ok_or(BedError::NoChromosomeOrder)?;

        // Get current position to restore later
        let start_pos = match &mut self.reader {
            BedReaderType::Plain(reader) => reader.stream_position()?,
            BedReaderType::Gzip(_) => {
                return Err(BedError::InvalidFormat(
                    "Cannot build skip index for gzip files".to_string(),
                ))
            }
            BedReaderType::Bgzf(_) => {
                return Err(BedError::InvalidFormat(
                    "Bgzf skip index not yet implemented".to_string(),
                ))
            }
        };

        // Get file size
        let file_size = match &mut self.reader {
            BedReaderType::Plain(reader) => {
                reader.seek(std::io::SeekFrom::End(0))?;
                let size = reader.stream_position()?;
                reader.seek(std::io::SeekFrom::Start(start_pos))?;
                size
            }
            _ => unreachable!(),
        };

        let chunk_size = file_size / (SKIP_CHUNKS + 1);
        let mut skip_index = SkipIndex::new();
        let mut buf = String::new();

        for i in 1..=SKIP_CHUNKS {
            let target_pos = i * chunk_size;

            // Seek to approximate position
            match &mut self.reader {
                BedReaderType::Plain(reader) => {
                    reader.seek(std::io::SeekFrom::Start(target_pos))?;

                    // Read until newline
                    if target_pos > 0 {
                        buf.clear();
                        reader.read_line(&mut buf)?;
                    }

                    // Get position after newline
                    let pos = reader.stream_position()?;

                    // Read next line to get chromosome
                    buf.clear();
                    if reader.read_line(&mut buf)? == 0 {
                        break; // EOF
                    }

                    // Parse chromosome from line
                    let fields: Vec<&str> = buf.trim_end().split('\t').collect();
                    if !fields.is_empty() {
                        let chrom = fields[0].to_string();
                        if let Some(&tid) = chromosome_order.get(&chrom) {
                            skip_index.entries.push(Skip::new(chrom, tid, pos));
                        } else {
                            return Err(BedError::ChromosomeNotFound);
                        }
                    }
                }
                _ => unreachable!(),
            }
        }

        // Restore original position
        match &mut self.reader {
            BedReaderType::Plain(reader) => {
                reader.seek(std::io::SeekFrom::Start(start_pos))?;
            }
            _ => unreachable!(),
        }

        skip_index.optimize();
        // check that entries are sorted by tid
        if !skip_index.entries.is_sorted_by_key(|skip| skip.tid) {
            return Err(BedError::ChromosomeOrderMismatch);
        }
        assert!(skip_index.entries.is_sorted_by_key(|skip| skip.seek_pos));

        // Store the index
        self.skip_index = Some(skip_index);

        Ok(())
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

            let record = match BedRecord::parse_line(line) {
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
pub struct Records<'a, R: Read> {
    reader: &'a mut BedReader<R>,
}

impl<R: Read> Iterator for Records<'_, R> {
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
                    BedRecord::parse_line(buf)
                        .expect("Failed to parse line")
                        .expect("Empty line to parse line")
                })
                .map_err(BedError::IoError)
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::{BedReader, BedRecord, BedValue, HashMap};
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::error::Error;
    use std::fs::File;
    use std::io::Cursor;
    use std::path::Path;

    fn setup_test_bed() -> Result<(), Box<dyn Error>> {
        use std::fs::File;
        use std::io::Write;

        let test_dir = Path::new("tests");
        if !test_dir.exists() {
            std::fs::create_dir(test_dir)?;
        }

        let path = test_dir.join("generated.bed");
        let mut file = File::create(path)?;

        // Write 10,000 entries for each chromosome from chr1 to chr22
        for chrom in 1..=22 {
            let chrom_name = format!("chr{}", chrom);

            for i in 0..5000 {
                let start = i * 100;
                let end = start + 50; // Each entry is 50 bases long
                let score = (i % 100) as f64; // Cycle scores from 0 to 99
                let name = format!("feat_{}_{}", chrom, i);

                writeln!(
                    file,
                    "{}\t{}\t{}\t{}\t{}",
                    chrom_name, start, end, name, score,
                )?;
            }
        }

        file.flush()?;
        Ok(())
    }

    #[test]
    fn test_skip_index() -> Result<(), Box<dyn Error>> {
        setup_test_bed()?;

        let path = Path::new("tests/generated.bed");
        let mut bed_reader = BedReader::<File>::from_path(path)?;
        bed_reader.set_chromosome_order(
            (1..=22)
                .map(|i| (format!("chr{}", i), i - 1))
                .collect::<HashMap<_, _>>(),
        );
        bed_reader.build_skip_index()?;
        eprintln!("{:?}", bed_reader.skip_index);
        /* this is no longer true because we optimize the skip index
        assert_eq!(
            bed_reader.skip_index.as_ref().unwrap().entries.len(),
            crate::SKIP_CHUNKS as usize
        );
        */

        let records: Vec<BedRecord> = bed_reader
            .query("chr19", 1, 102)?
            .collect::<Result<Vec<_>, _>>()?;
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom(), "chr19");
        assert_eq!(records[1].chrom(), "chr19");

        let records: Vec<BedRecord> = bed_reader
            .query("chr19", 1, 10)?
            .collect::<Result<Vec<_>, _>>()?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chrom(), "chr19");

        // Time 1000 queries and calculate queries per second
        let start = std::time::Instant::now();
        let iterations = 1000;

        for _ in 0..iterations {
            let records: Vec<BedRecord> = bed_reader
                .query("chr19", 1, 10)?
                .collect::<Result<Vec<_>, _>>()?;
            assert_eq!(records.len(), 1);
            assert_eq!(records[0].chrom(), "chr19");
        }

        let duration = start.elapsed();
        let qps = iterations as f64 / duration.as_secs_f64();
        println!(
            "Performed {} queries in {:?} ({:.2} queries/sec)",
            iterations, duration, qps
        );

        Ok(())
    }

    #[test]
    fn test_read_bed_file() -> Result<(), Box<dyn Error>> {
        let bed_path = Path::new("tests/test.bed");
        let mut bed_reader = BedReader::<File>::from_path(bed_path)?;

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
        let mut bed_reader = BedReader::<File>::from_path(bed_path)?;
        bed_reader.set_chromosome_order(
            (1..=22)
                .map(|i| (format!("chr{}", i), i - 1))
                .collect::<HashMap<_, _>>(),
        );

        let records: Vec<BedRecord> = bed_reader
            .query("chr1", 22, 34)?
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
        let mut bed_reader = BedReader::<File>::from_path(bed_path)?;
        bed_reader.set_chromosome_order(
            (1..=22)
                .map(|i| (format!("chr{}", i), i - 1))
                .collect::<HashMap<_, _>>(),
        );

        let records: Vec<BedRecord> = bed_reader
            .query("chr1", 1, 3)?
            .collect::<Result<Vec<_>, _>>()?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chrom(), "chr1");
        assert_eq!(records[0].start(), 1);
        assert_eq!(records[0].end(), 10);
        Ok(())
    }

    #[test]
    fn test_query_0_start() -> Result<(), Box<dyn Error>> {
        let bed_path = Path::new("tests/compr.bed.gz");
        let mut bed_reader = BedReader::<File>::from_path(bed_path)?;
        bed_reader.set_chromosome_order(
            (1..=22)
                .map(|i| (format!("chr{}", i), i - 1))
                .collect::<HashMap<_, _>>(),
        );

        let records: Vec<BedRecord> = bed_reader
            .query("chr1", 0, 10)?
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
        let mut bed_reader = BedReader::<File>::from_path(bed_path)?;
        eprintln!("{:?}", bed_reader);
        bed_reader.set_chromosome_order(
            (1..=22)
                .map(|i| (format!("chr{}", i), i - 1))
                .collect::<HashMap<_, _>>(),
        );
        eprintln!("{:?}", bed_reader.chromosome_order);

        let records: Vec<BedRecord> = bed_reader
            .query("chr1", 999998, 999999)?
            .collect::<Result<Vec<_>, _>>()?;
        assert_eq!(records.len(), 1);
        Ok(())
    }

    #[test]
    fn test_query_first_interval_multiple_formats() -> Result<(), Box<dyn Error>> {
        // Test data
        let test_data = b"chr1\t1\t10\tgene1\n\
                         chr1\t20\t30\tgene2\n\
                         chr2\t5\t15\tgene3\n";

        use crate::BedReaderType;
        use crate::GzDecoder;
        use std::io::BufReader;
        use std::io::Write;

        // Test uncompressed data using Cursor
        let cursor = Cursor::new(test_data.to_vec());
        let reader = BufReader::new(cursor);
        let mut bed_reader = BedReader::from_reader(BedReaderType::Plain(reader), "memory")?;
        bed_reader.set_chromosome_order(
            (1..=22)
                .map(|i| (format!("chr{}", i), i - 1))
                .collect::<HashMap<_, _>>(),
        );

        let records: Vec<BedRecord> = bed_reader
            .query("chr1", 1, 3)?
            .collect::<Result<Vec<_>, _>>()?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chrom(), "chr1");
        assert_eq!(records[0].start(), 1);
        assert_eq!(records[0].end(), 10);

        // Test gzipped data
        let mut gz_buf = Vec::new();
        {
            let mut encoder = GzEncoder::new(&mut gz_buf, Compression::default());
            encoder.write_all(test_data)?;
        }
        let cursor = Cursor::new(gz_buf.to_vec());
        let reader = BufReader::new(GzDecoder::new(BufReader::new(cursor)));
        let mut bed_reader = BedReader::from_reader(BedReaderType::Gzip(reader), "gz_memory")?;
        bed_reader.set_chromosome_order(
            (1..=22)
                .map(|i| (format!("chr{}", i), i - 1))
                .collect::<HashMap<_, _>>(),
        );

        let records: Vec<BedRecord> = bed_reader
            .query("chr1", 1, 3)?
            .collect::<Result<Vec<_>, _>>()?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chrom(), "chr1");
        assert_eq!(records[0].start(), 1);
        assert_eq!(records[0].end(), 10);

        Ok(())
    }

    #[test]
    fn test_chrom_order_from_index() -> Result<(), Box<dyn Error>> {
        // Test with indexed file (using the existing compressed test file)
        let bed_path = Path::new("tests/compr.bed.gz");
        let bed_reader = BedReader::<File>::from_path(bed_path)?;

        // Get chromosome order from index
        let chrom_order = bed_reader.chrom_order_from_index();

        // Verify that we got a chromosome order (since compr.bed.gz should have an index)
        assert!(chrom_order.is_some());

        if let Some(order) = chrom_order {
            // Verify some expected chromosomes are present with correct ordering
            assert_eq!(order.get("chr1"), Some(&0));
            // Add more assertions based on your test data
        }

        // Test with non-indexed file
        let bed_path = Path::new("tests/test.bed");
        let bed_reader = BedReader::<File>::from_path(bed_path)?;

        // Get chromosome order from non-indexed file
        let chrom_order = bed_reader.chrom_order_from_index();

        // Verify that we get None for a non-indexed file
        assert!(chrom_order.is_none());

        Ok(())
    }
}
