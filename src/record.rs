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

    /// Sets the start coordinate (0-based).
    pub fn set_start(&mut self, start: u64) {
        self.start = start;
    }

    /// Sets the end coordinate (exclusive).
    pub fn set_end(&mut self, end: u64) {
        self.end = end;
    }

    /// Returns the name, if present.
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Returns the score, if present.
    pub fn score(&self) -> Option<f64> {
        self.score
    }

    /// Sets the score
    pub fn set_score(&mut self, score: f64) {
        self.score = Some(score);
    }

    /// Sets the name
    pub fn set_name(&mut self, name: String) {
        self.name = Some(name);
    }

    /// Returns a slice of the other fields.
    pub fn other_fields(&self) -> &[BedValue] {
        &self.other_fields
    }

    /// Adds a new field to other_fields
    pub fn push_field(&mut self, value: BedValue) {
        self.other_fields.push(value);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_start_end() {
        let mut record = BedRecord::new(
            "chr1".to_string(),
            1000,
            2000,
            Some("feature1".to_string()),
            Some(100.0),
            vec![],
        );

        // Test set_start
        record.set_start(1500);
        assert_eq!(record.start(), 1500);

        // Test set_end
        record.set_end(2500);
        assert_eq!(record.end(), 2500);

        // Test multiple sets
        record.set_start(1600);
        record.set_end(2600);
        assert_eq!(record.start(), 1600);
        assert_eq!(record.end(), 2600);
    }

    #[test]
    fn test_push_field() {
        let mut record = BedRecord::new("chr1".to_string(), 1000, 2000, None, None, vec![]);

        record.push_field(BedValue::String("test_value".to_string()));
        assert_eq!(record.other_fields().len(), 1);
        assert_eq!(
            record.other_fields()[0],
            BedValue::String("test_value".to_string())
        );
    }

    #[test]
    fn test_set_name() {
        let mut record = BedRecord::new("chr1".to_string(), 1000, 2000, None, None, vec![]);

        // Test setting name on a record with no name
        record.set_name("feature1".to_string());
        assert_eq!(record.name(), Some("feature1"));

        // Test updating an existing name
        record.set_name("feature2".to_string());
        assert_eq!(record.name(), Some("feature2"));
    }

    #[test]
    fn test_set_score() {
        let mut record = BedRecord::new("chr1".to_string(), 1000, 2000, None, None, vec![]);

        // Test setting score on a record with no score
        record.set_score(45.5);
        assert_eq!(record.score(), Some(45.5));

        // Test updating an existing score
        record.set_score(100.0);
        assert_eq!(record.score(), Some(100.0));
    }
}
