//! Indexing for BED files
use noodles::{
    core::Region,
    csi::{self, binning_index::index::reference_sequence::bin::Chunk, BinningIndex},
    tabix,
};
use std::io;

/// Represents the type of index available for the BED file
#[derive(Debug)]
pub enum BedIndex {
    /// Tabix index
    Tabix(tabix::Index),
    /// CSI index
    Csi(csi::Index),
    /// No index
    None,
}

impl BedIndex {
    /// Query the index for the given chromosome and region
    pub fn query(&self, tid: usize, region: &Region) -> Option<io::Result<Vec<Chunk>>> {
        match self {
            BedIndex::Tabix(index) => Some(index.query(tid, region.interval())),
            BedIndex::Csi(index) => Some(index.query(tid, region.interval())),
            BedIndex::None => None,
        }
    }

    /// Return the csi header which contains the reference sequence names
    pub fn header(&self) -> Option<&csi::binning_index::index::Header> {
        match self {
            BedIndex::Tabix(index) => index.header(),
            BedIndex::Csi(index) => index.header(),
            BedIndex::None => None,
        }
    }
}

/// A skip index entry storing chromosome and file position information
#[allow(dead_code)] // chromosome field is not used but useful for debugging
#[derive(Debug, PartialEq, Eq)]
pub struct Skip {
    chromosome: String,
    pub(crate) tid: usize,
    pub(crate) seek_pos: u64,
}

impl Skip {
    pub(crate) fn new(chromosome: String, tid: usize, seek_pos: u64) -> Self {
        Skip {
            chromosome,
            tid,
            seek_pos,
        }
    }
}

/// A simple index that allows skipping to approximate chromosome positions
#[derive(Debug)]
pub struct SkipIndex {
    pub(crate) entries: Vec<Skip>,
}

impl SkipIndex {
    pub(crate) fn new() -> Self {
        SkipIndex {
            entries: Vec::new(),
        }
    }

    /// Optimize the index by keeping only the first and last entries for each tid
    pub(crate) fn optimize(&mut self) {
        // Early return if we have 2 or fewer entries
        if self.entries.len() <= 2 {
            return;
        }

        // Group entries by tid
        let mut current_tid = self.entries[0].tid;
        let mut first_idx = 0;
        let mut to_remove = Vec::new();

        for i in 1..self.entries.len() {
            if self.entries[i].tid != current_tid {
                // We've found a new tid. Mark middle entries of previous tid for removal
                if i - first_idx > 2 {
                    // Add all indices except first and last for this tid group
                    to_remove.extend((first_idx + 1)..(i - 1));
                }
                // Update for new tid group
                current_tid = self.entries[i].tid;
                first_idx = i;
            }
        }

        // Handle the last group
        let last_idx = self.entries.len() - 1;
        if last_idx - first_idx > 2 {
            to_remove.extend((first_idx + 1)..(last_idx));
        }

        // Remove marked entries in reverse order to maintain valid indices
        for idx in to_remove.iter().rev() {
            self.entries.remove(*idx);
        }
    }

    /// Binary search for the closest position before the target tid
    pub fn find_seek_position(&self, target_tid: usize) -> Option<u64> {
        let idx = self.entries.partition_point(|skip| skip.tid < target_tid);
        if idx == 0 {
            None
        } else {
            eprintln!(
                "found seek position for tid {} with entry: {:?}",
                target_tid,
                self.entries[idx - 1]
            );
            Some(self.entries[idx - 1].seek_pos)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_skip_index_optimize() {
        let mut index = SkipIndex {
            entries: vec![
                Skip::new("chr1".to_string(), 0, 100),
                Skip::new("chr1".to_string(), 0, 200),
                Skip::new("chr1".to_string(), 0, 300),
                Skip::new("chr1".to_string(), 0, 400),
                Skip::new("chr2".to_string(), 1, 500),
                Skip::new("chr2".to_string(), 1, 600),
                Skip::new("chr2".to_string(), 1, 700),
                Skip::new("chr3".to_string(), 2, 800),
                Skip::new("chr4".to_string(), 3, 900),
                Skip::new("chr4".to_string(), 3, 1000),
                Skip::new("chr4".to_string(), 3, 1100),
                Skip::new("chr4".to_string(), 3, 1200),
            ],
        };

        index.optimize();
        eprintln!("optimized index: {:?}", index);

        // Check that we kept only first and last entries for each tid
        assert_eq!(index.entries.len(), 7);

        assert_eq!(index.entries[0], Skip::new("chr1".to_string(), 0, 100));
        assert_eq!(index.entries[1], Skip::new("chr1".to_string(), 0, 400));
        assert_eq!(index.entries[2], Skip::new("chr2".to_string(), 1, 500));
        assert_eq!(index.entries[3], Skip::new("chr2".to_string(), 1, 700));
        assert_eq!(index.entries[4], Skip::new("chr3".to_string(), 2, 800));
        assert_eq!(index.entries[5], Skip::new("chr4".to_string(), 3, 900));
        assert_eq!(index.entries[6], Skip::new("chr4".to_string(), 3, 1200));
    }

    #[test]
    fn test_skip_index_optimize_edge_cases() {
        // Test case with single entry per tid
        let mut index = SkipIndex {
            entries: vec![
                Skip::new("chr1".to_string(), 0, 100),
                Skip::new("chr2".to_string(), 1, 200),
                Skip::new("chr3".to_string(), 2, 300),
            ],
        };

        index.optimize();
        assert_eq!(index.entries.len(), 3);
        assert_eq!(index.entries[0], Skip::new("chr1".to_string(), 0, 100));
        assert_eq!(index.entries[1], Skip::new("chr2".to_string(), 1, 200));
        assert_eq!(index.entries[2], Skip::new("chr3".to_string(), 2, 300));

        // Test case with exactly two entries per tid
        let mut index = SkipIndex {
            entries: vec![
                Skip::new("chr1".to_string(), 0, 100),
                Skip::new("chr1".to_string(), 0, 200),
                Skip::new("chr2".to_string(), 1, 300),
                Skip::new("chr2".to_string(), 1, 400),
            ],
        };

        index.optimize();
        assert_eq!(index.entries.len(), 4);
        assert_eq!(index.entries[0], Skip::new("chr1".to_string(), 0, 100));
        assert_eq!(index.entries[1], Skip::new("chr1".to_string(), 0, 200));
        assert_eq!(index.entries[2], Skip::new("chr2".to_string(), 1, 300));
        assert_eq!(index.entries[3], Skip::new("chr2".to_string(), 1, 400));

        // Test empty and single-entry cases
        let mut empty_index = SkipIndex::new();
        empty_index.optimize();
        assert_eq!(empty_index.entries.len(), 0);

        let mut single_entry = SkipIndex {
            entries: vec![Skip::new("chr1".to_string(), 0, 100)],
        };
        single_entry.optimize();
        assert_eq!(single_entry.entries.len(), 1);
        assert_eq!(
            single_entry.entries[0],
            Skip::new("chr1".to_string(), 0, 100)
        );
    }
}
