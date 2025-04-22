//! A writer for BED files.
use crate::BedError;
use crate::BedRecord;
use noodles::bgzf;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// A writer for BED files, supporting plain text and BGZF compression.
pub struct BedWriter {
    writer: Box<dyn Write>,
}

impl BedWriter {
    /// Creates a new `BedWriter` from a file path.
    /// Automatically detects BGZF compression based on file extension (.gz or .bgz).
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, BedError> {
        let path_ref = path.as_ref();
        let file = File::create(path_ref)?;
        let writer: Box<dyn Write> = if path_ref.extension().and_then(|s| s.to_str()) == Some("gz")
        {
            let bgzf_writer = bgzf::Writer::new(file);
            Box::new(bgzf_writer) // Box the Writer<bgzf::Writer>
        } else {
            Box::new(BufWriter::new(file)) // Box the Writer<File>
        };
        Self::from_writer(writer)
    }

    /// Creates a new `BedWriter` from a writer.
    pub fn from_writer(writer: Box<dyn Write>) -> Result<Self, BedError> {
        Ok(BedWriter { writer })
    }

    /// Writes a `BedRecord` to the writer.
    pub fn write_record(&mut self, record: &BedRecord) -> Result<(), BedError> {
        write!(self.writer, "{}", record.chrom())?;
        write!(self.writer, "\t{}", record.start())?;
        write!(self.writer, "\t{}", record.end())?;

        if let Some(name) = record.name() {
            write!(self.writer, "\t{}", name)?;
        }

        if let Some(score) = record.score() {
            write!(self.writer, "\t{}", score)?;
        }

        for field in record.other_fields() {
            write!(self.writer, "\t{}", field)?;
        }

        writeln!(self.writer)?; // Newline after each record
        Ok(())
    }

    /// Flushes the writer to ensure all data is written to the underlying file.
    pub fn flush(&mut self) -> Result<(), BedError> {
        self.writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BedReader, BedRecord, BedValue};
    use std::error::Error;
    use std::fs;

    #[test]
    fn test_write_bed_file() -> Result<(), Box<dyn Error>> {
        let output_path = Path::new("output_test.bed");
        let mut bed_writer = BedWriter::new(output_path)?;

        let record1 = BedRecord::new(
            "chr1".to_string(),
            100,
            200,
            Some("test_record".to_string()),
            Some(99.0),
            vec![BedValue::String("extra_field".to_string())],
        );
        bed_writer.write_record(&record1)?;

        let record2 = BedRecord::new("chr2".to_string(), 300, 400, None, None, vec![]);
        bed_writer.write_record(&record2)?;
        bed_writer.flush()?;

        let mut bed_reader = BedReader::<File>::from_path(output_path)?;
        let read_record1 = bed_reader.read_record()?.unwrap();
        assert_eq!(read_record1, record1);
        let read_record2 = bed_reader.read_record()?.unwrap();
        assert_eq!(read_record2, record2);

        fs::remove_file(output_path)?; // Cleanup test file
        Ok(())
    }

    #[test]
    fn test_write_bgzf_bed_file() -> Result<(), Box<dyn Error>> {
        let output_path = Path::new("output_test.bed.gz");
        let mut bed_writer = BedWriter::new(output_path)?;

        let record1 = BedRecord::new(
            "chr1".to_string(),
            100,
            200,
            Some("test_record".to_string()),
            Some(99.0),
            vec![BedValue::String("extra_field".to_string())],
        );
        bed_writer.write_record(&record1)?;

        let record2 = BedRecord::new("chr2".to_string(), 300, 400, None, None, vec![]);
        bed_writer.write_record(&record2)?;
        bed_writer.flush()?;

        let mut bed_reader = BedReader::<File>::from_path(output_path)?;
        let read_record1 = bed_reader.read_record()?.unwrap();
        assert_eq!(read_record1, record1);
        let read_record2 = bed_reader.read_record()?.unwrap();
        assert_eq!(read_record2, record2);

        fs::remove_file(output_path)?; // Cleanup test file
        Ok(())
    }
}
