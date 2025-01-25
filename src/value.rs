use std::fmt;

/// Represents the possible data types for optional BED fields.
#[derive(Debug, Clone, PartialEq)]
pub enum BedValue {
    String(String),
    Int(i64),
    Float(f64),
}

impl BedValue {
    /// Attempts to parse a string slice into a `BedValue`.
    /// Tries parsing as Int, then Float, then falls back to String.
    pub fn parse(s: &str) -> Self {
        if let Ok(int_val) = s.parse::<i64>() {
            BedValue::Int(int_val)
        } else if let Ok(float_val) = s.parse::<f64>() {
            BedValue::Float(float_val)
        } else {
            BedValue::String(s.to_string())
        }
    }
}

impl fmt::Display for BedValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BedValue::String(s) => write!(f, "{}", s),
            BedValue::Int(i) => write!(f, "{}", i),
            BedValue::Float(fl) => write!(f, "{:.4}", fl),
        }
    }
}
