//! Represents a BED record.
use crate::value::BedValue;
use crate::BedError;

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
    #[inline]
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// Returns the start coordinate (0-based).
    #[inline]
    pub fn start(&self) -> u64 {
        self.start
    }

    /// Returns the end coordinate (exclusive).
    #[inline]
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

    /// Parse a BED record from a string.
    pub fn parse_line(line: &str) -> Result<Option<BedRecord>, BedError> {
        let mut fields = line.split('\t');

        // Get required fields
        let chrom = match fields.next() {
            Some(c) => c.to_string(),
            None => {
                return Err(BedError::InvalidFormat(
                    "BED format requires at least 3 columns (chrom, start, end)".to_string(),
                ))
            }
        };

        let start = match fields.next() {
            Some(s) => s.parse::<u64>().map_err(|e| {
                BedError::ParseError(format!("Invalid start coordinate: {}. Error: {}", s, e))
            })?,
            None => {
                return Err(BedError::InvalidFormat(
                    "BED format requires at least 3 columns (chrom, start, end)".to_string(),
                ))
            }
        };

        let end = match fields.next() {
            Some(e) => e.parse::<u64>().map_err(|e| {
                BedError::ParseError(format!("Invalid end coordinate: {}. Error: {}", e, e))
            })?,
            None => {
                return Err(BedError::InvalidFormat(
                    "BED format requires at least 3 columns (chrom, start, end)".to_string(),
                ))
            }
        };

        // Optional fields
        let name = fields.next().filter(|f| !f.is_empty()).map(String::from);

        let score = fields
            .next()
            .filter(|f| !f.is_empty())
            .and_then(|s| s.parse::<f64>().ok());

        // Remaining fields
        let other_fields: Vec<BedValue> = fields.map(BedValue::parse).collect();

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
