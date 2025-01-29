use crate::value::BedValue;

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
