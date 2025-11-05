//! Basic building blocks of an alignment.

use crate::utils;

use std::ops::Range;

//-----------------------------------------------------------------------------

/// An alignment between a substring of a query sequence and a single node using a single edit operation.
///
/// This can also represent unaligned insertions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mapping {
    // Aligned interval of the query sequence.
    seq_interval: Range<usize>,
    // Handle of the target node, or `gbwt::ENDMARKER` for no target node.
    handle: usize,
    // Length of the target node.
    node_len: usize,
    // Aligned interval of the target node.
    node_interval: Range<usize>,
    // Edit operation.
    edit: Difference,
}

impl Mapping {
    /// Creates a new mapping starting at the given position.
    ///
    /// # Panics
    ///
    /// Will panic if the edit extends past the end of the node.
    /// Will panic if `handle` is [`gbwt::ENDMARKER`].
    pub fn new(seq_offset: usize, handle: usize, node_len: usize, node_offset: usize, edit: Difference) -> Self {
        assert!(node_offset + edit.target_len() <= node_len, "The edit extends past the end of the node");
        assert!(handle != gbwt::ENDMARKER, "A mapping cannot have {} as the target handle", gbwt::ENDMARKER);
        Self {
            seq_interval: seq_offset..seq_offset + edit.query_len(),
            handle,
            node_len,
            node_interval: node_offset..node_offset + edit.target_len(),
            edit,
        }
    }

    /// Creates a mapping for an unaligned insertion.
    ///
    /// # Panics
    ///
    /// Will panic if the edit is not an insertion.
    pub fn unaligned(seq_offset: usize, edit: Difference) -> Self {
        assert!(matches!(edit, Difference::Insertion(_)), "An unaligned mapping must have an insertion edit");
        Self {
            seq_interval: seq_offset..seq_offset + edit.query_len(),
            handle: gbwt::ENDMARKER,
            node_len: 0,
            node_interval: 0..0,
            edit,
        }
    }

    /// Returns the aligned interval of the query sequence.
    #[inline]
    pub fn seq_interval(&self) -> &Range<usize> {
        &self.seq_interval
    }

    /// Returns the length of the query interval.
    #[inline]
    pub fn query_len(&self) -> usize {
        self.seq_interval.len()
    }

    /// Returns the handle of the target node.
    ///
    /// Returns [`gbwt::ENDMARKER`] for an unaligned insertion.
    #[inline]
    pub fn handle(&self) -> usize {
        self.handle
    }

    /// Returns `true` if this mapping is an unaligned insertion.
    #[inline]
    pub fn is_unaligned(&self) -> bool {
        self.handle == gbwt::ENDMARKER
    }

    /// Returns the length of the target node.
    #[inline]
    pub fn node_len(&self) -> usize {
        self.node_len
    }

    /// Returns the aligned interval of the target node.
    #[inline]
    pub fn node_interval(&self) -> &Range<usize> {
        &self.node_interval
    }

    /// Returns the length of the target interval.
    #[inline]
    pub fn target_len(&self) -> usize {
        self.node_interval.len()
    }

    /// Returns the edit operation.
    #[inline]
    pub fn edit(&self) -> &Difference {
        &self.edit
    }

    /// Returns `true` this mapping follows the given one in the same node.
    #[inline]
    pub fn follows(&self, other: &Self) -> bool {
        self.handle == other.handle && self.node_interval.start == other.node_interval.end
    }

    /// Returns `true` if this mapping is at the start of the target node.
    #[inline]
    pub fn is_at_start(&self) -> bool {
        self.node_interval.start == 0
    }

    /// Returns `true` if this mapping is at the end of the target node.
    #[inline]
    pub fn is_at_end(&self) -> bool {
        self.node_interval.end == self.node_len
    }
}

//-----------------------------------------------------------------------------

/// An operation in a difference string describing an alignment between a query sequence and a target sequence.
///
/// This implementation supports the following operations:
///
/// * `=`: A match given as the matching sequence.
/// * `:`: A match given as the match length.
/// * `*`: A mismatch given as the target base and the query base.
/// * `+`: An insertion given as the inserted sequence.
/// * `-`: A deletion given as the deleted sequence.
///
/// The operations do not store target bases, as the query sequence can be reconstructed without that information.
/// Operation `~` (intron length and splice signal) is not supported yet.
/// Parsing is based on bytes rather than characters to avoid unnecessary UTF-8 validation.
///
/// # Examples
///
/// ```
/// use gbz_base::Difference;
///
/// let with_gaps = b":48-CAT:44+GATTACA:51";
/// let ops = Difference::parse(with_gaps);
/// assert!(ops.is_ok());
/// let ops = ops.unwrap();
/// assert_eq!(ops.len(), 5);
/// assert_eq!(ops[0], Difference::Match(48));
/// assert_eq!(ops[1], Difference::Deletion(3));
/// assert_eq!(ops[2], Difference::Match(44));
/// assert_eq!(ops[3], Difference::Insertion(b"GATTACA".to_vec()));
/// assert_eq!(ops[4], Difference::Match(51));
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Difference {
    /// A match of the given length.
    Match(usize),
    /// Mismatch represented as the query base.
    Mismatch(u8),
    /// Insertion to the reference represented as the inserted sequence.
    Insertion(Vec<u8>),
    /// Deletion from the reference represented as deletion length.
    Deletion(usize),
    /// Marker for the end of a sequence when concatenating multiple difference strings.
    End,
}

impl Difference {
    /// Number of supported operation types.
    pub const NUM_TYPES: usize = 5;

    /// Symbol used as a substitute for an unknown base.
    pub const UNKNOWN_BASE: u8 = b'X';

    // TODO: This does not support `~` (intron length and splice signal) yet.
    const OPS: &'static [u8] = b"=:*+-";

    fn base_to_upper(c: u8) -> u8 {
        if c.is_ascii_lowercase() {
            c - b'a' + b'A'
        } else {
            c
        }
    }

    fn seq_to_upper(seq: &[u8]) -> Vec<u8> {
        seq.iter().map(|&c| Self::base_to_upper(c)).collect()
    }

    fn matching_sequence(value: &[u8]) -> Option<Self> {
        Some(Self::Match(value.len()))
    }

    fn match_length(value: &[u8]) -> Option<Self> {
        let len = str::from_utf8(value).ok()?;
        let len = len.parse::<usize>().ok()?;
        Some(Self::Match(len))
    }

    fn mismatch(value: &[u8]) -> Option<Self> {
        if value.len() != 2 {
            return None;
        }
        Some(Self::Mismatch(Self::base_to_upper(value[1])))
    }

    fn insertion(value: &[u8]) -> Option<Self> {
        Some(Self::Insertion(Self::seq_to_upper(value)))
    }

    fn deletion(value: &[u8]) -> Option<Self> {
        Some(Self::Deletion(value.len()))
    }

    /// Parses a difference string and returns it as a vector of operations.
    ///
    /// Returns an error if the difference string is invalid.
    pub fn parse(difference_string: &[u8]) -> Result<Vec<Self>, String> {
        let mut result: Vec<Self> = Vec::new();
        if difference_string.is_empty() {
            return Ok(result);
        }
        if !Self::OPS.contains(&difference_string[0]) {
            return Err(format!("Invalid difference string operation: {}", difference_string[0] as char));
        }

        let mut start = 0;
        while start < difference_string.len() {
            let mut end = start + 1;
            while end < difference_string.len() && !Self::OPS.contains(&difference_string[end]) {
                end += 1;
            }
            let value = &difference_string[start + 1..end];
            let op = match difference_string[start] {
                b'=' => Self::matching_sequence(value),
                b':' => Self::match_length(value),
                b'*' => Self::mismatch(value),
                b'+' => Self::insertion(value),
                b'-' => Self::deletion(value),
                _ => return Err(format!("Invalid difference string operation: {}", difference_string[start] as char)),
            }.ok_or(format!("Invalid difference string field: {}", String::from_utf8_lossy(&difference_string[start..end])))?;
            result.push(op);
            start = end;
        }

        Ok(result)
    }

    /// Parses a difference string and returns it as a normalized vector of operations.
    ///
    /// The operations are merged and empty operations are removed.
    /// Returns an error if the difference string is invalid.
    pub fn parse_normalized(difference_string: &[u8]) -> Result<Vec<Self>, String> {
        let ops = Self::parse(difference_string)?;
        Ok(Self::normalize(ops))
    }

    /// Calculates various statistics from a sequence of operations.
    ///
    /// The return value is (query length, target length, matches, edits).
    pub fn stats(difference_string: &[Self]) -> (usize, usize, usize, usize) {
        let mut query_len = 0;
        let mut target_len = 0;
        let mut matches = 0;
        let mut edits = 0;
        for diff in difference_string.iter() {
            match diff {
                Self::Match(len) => {
                    query_len += len;
                    target_len += len;
                    matches += len;
                },
                Self::Mismatch(_) => {
                    query_len += 1;
                    target_len += 1;
                    edits += 1;
                }
                Self::Insertion(seq) => {
                    query_len += seq.len();
                    edits += seq.len();
                }
                Self::Deletion(len) => {
                    target_len += len;
                    edits += len;
                },
                Self::End => {},
            }
        }
        (query_len, target_len, matches, edits)
    }

    /// Returns the length of the operation.
    pub fn len(&self) -> usize {
        match self {
            Self::Match(len) => *len,
            Self::Mismatch(_) => 1,
            Self::Insertion(seq) => seq.len(),
            Self::Deletion(len) => *len,
            Self::End => 0,
        }
    }

    /// Returns `true` if the operation is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the length of the operation in the target sequence.
    pub fn target_len(&self) -> usize {
        match self {
            Self::Match(len) => *len,
            Self::Mismatch(_) => 1,
            Self::Insertion(_) => 0,
            Self::Deletion(len) => *len,
            Self::End => 0,
        }
    }

    /// Returns the length of the operation in the query sequence.
    pub fn query_len(&self) -> usize {
        match self {
            Self::Match(len) => *len,
            Self::Mismatch(_) => 1,
            Self::Insertion(seq) => seq.len(),
            Self::Deletion(_) => 0,
            Self::End => 0,
        }
    }

    /// Merges the given operation into this operation if they can be merged.
    ///
    /// Returns `true` if the operations were merged.
    pub fn try_merge(&mut self, op: &Self) -> bool {
        match (self, op) {
            (Self::Match(len1), Self::Match(len2)) => {
                *len1 += len2;
                true
            },
            (Self::Insertion(seq1), Self::Insertion(seq2)) => {
                seq1.extend_from_slice(seq2);
                true
            },
            (Self::Deletion(len1), Self::Deletion(len2)) => {
                *len1 += len2;
                true
            },
            _ => false,
        }
    }

    /// Normalizes the sequence of operations.
    ///
    /// This merges adjacent matches and insertions and removes empty operations.
    pub fn normalize(ops: Vec<Self>) -> Vec<Self> {
        let mut result = ops;
        let mut tail = 0;
        for i in 0..result.len() {
            if result[i].is_empty() {
                continue;
            }
            if tail > 0 {
                // We need to be careful around the borrow checker.
                let (left, right) = result.split_at_mut(i);
                if left[tail - 1].try_merge(&right[0]) {
                    continue;
                }
            }
            result.swap(tail, i);
            tail += 1;
        }

        result.truncate(tail);
        result
    }

    /// Writes a difference string as a `Vec<u8>` string.
    pub fn to_bytes(ops: &[Difference], target_sequence: &[u8]) -> Vec<u8> {
        let mut result = Vec::new();
        let mut target_offset = 0;
        for op in ops.iter() {
            match op {
                Self::Match(len) => {
                    result.push(b':');
                    utils::append_usize(&mut result, *len);
                    target_offset += *len;
                },
                Self::Mismatch(base) => {
                    result.push(b'*');
                    result.push(target_sequence[target_offset]);
                    result.push(*base);
                    target_offset += 1;
                },
                Self::Insertion(seq) => {
                    result.push(b'+');
                    result.extend_from_slice(seq);
                },
                Self::Deletion(len) => {
                    result.push(b'-');
                    result.extend_from_slice(&target_sequence[target_offset..target_offset + *len]);
                    target_offset += *len;
                },
                Self::End => {},
            }
        }
        result
    }
}

//-----------------------------------------------------------------------------

#[cfg(test)]
mod tests {

use super::*;

// Tests for `Difference`.

fn check_difference(seq: &[u8], truth: &[Difference], name: &str, normalize: bool) {
    let result = Difference::parse(seq);
    assert!(result.is_ok(), "Failed to parse the difference string for {}: {}", name, result.err().unwrap());
    let mut result = result.unwrap();
    if normalize {
        result = Difference::normalize(result);
    }
    assert_eq!(result.len(), truth.len(), "Wrong number of operations for {}", name);
    for i in 0..truth.len() {
        assert_eq!(result[i], truth[i], "Wrong operation at position {} for {}", i, name);
    }
}

fn check_to_bytes(ops: &[Difference], truth: &[u8], target_sequence: &[u8], name: &str) {
    let result = Difference::to_bytes(ops, target_sequence);
    assert_eq!(result, truth, "Wrong byte representation for {}", name);
}

#[test]
fn difference_empty() {
    let name = "empty";
    let seq = b"";
    let truth = [];
    check_difference(seq, &truth, name, false);
    let target_sequence = b"";
    check_to_bytes(&truth, seq, target_sequence, name);
}

#[test]
fn difference_single() {
    {
        let name = "matching sequence";
        let seq = b"=ACGT";
        let truth = [ Difference::Match(4) ];
        check_difference(seq, &truth, name, false);
        let expected = b":4"; // Converted to match length.
        let target_sequence = b"ACGT";
        check_to_bytes(&truth, expected, target_sequence, name);
    }
    {
        let name = "match length";
        let seq = b":5";
        let truth = [ Difference::Match(5) ];
        check_difference(seq, &truth, name, false);
        let target_sequence = b"GATTA";
        check_to_bytes(&truth, seq, target_sequence, name);
    }
    {
        let name = "mismatch";
        let seq = b"*AC";
        let truth = [ Difference::Mismatch(b'C') ];
        check_difference(seq, &truth, name, false);
        let target_sequence = b"A";
        check_to_bytes(&truth, seq, target_sequence, name);
    }
    {
        let name = "insertion";
        let seq = b"+ACGT";
        let truth = [ Difference::Insertion(b"ACGT".to_vec()) ];
        check_difference(seq, &truth, name, false);
        let target_sequence = b"";
        check_to_bytes(&truth, seq, target_sequence, name);
    }
    {
        let name = "deletion";
        let seq = b"-ACGT";
        let truth = [ Difference::Deletion(4) ];
        check_difference(seq, &truth, name, false);
        let target_sequence = b"ACGT";
        check_to_bytes(&truth, seq, target_sequence, name);
    }
}

#[test]
fn difference_multi() {
    // All of these are from real reads mapped with Giraffe.
    // Target sequences only contain the relevant bases.
    {
        let name = "ERR3239454.129477191 fragment 2";
        let seq = b":24*CG:78-C:35*TG:2*TG:1*TG:6";
        let truth = [
            Difference::Match(24),
            Difference::Mismatch(b'G'),
            Difference::Match(78),
            Difference::Deletion(1),
            Difference::Match(35),
            Difference::Mismatch(b'G'),
            Difference::Match(2),
            Difference::Mismatch(b'G'),
            Difference::Match(1),
            Difference::Mismatch(b'G'),
            Difference::Match(6)
        ];
        check_difference(seq, &truth, name, false);
        let mut target_sequence = vec![b'-'; 151];
        target_sequence[24] = b'C';
        target_sequence[103] = b'C';
        target_sequence[139] = b'T';
        target_sequence[142] = b'T';
        target_sequence[144] = b'T';
        check_to_bytes(&truth, seq, &target_sequence, name);
    }
    {
        let name = "ERR3239454.11251898 fragment 1";
        let seq = b":129+TTGAGGGGGTATAAGAGGTCG";
        let truth = [
            Difference::Match(129),
            Difference::Insertion(b"TTGAGGGGGTATAAGAGGTCG".to_vec())
        ];
        check_difference(seq, &truth, name, false);
        let target_sequence = vec![b'-', 129];
        check_to_bytes(&truth, seq, &target_sequence, name);
    }
    {
        let name = "ERR3239454.97848632 fragment 1";
        let seq = b":111*GT:20*CA:6*GT:2*GT:3*GT:3";
        let truth = [
            Difference::Match(111),
            Difference::Mismatch(b'T'),
            Difference::Match(20),
            Difference::Mismatch(b'A'),
            Difference::Match(6),
            Difference::Mismatch(b'T'),
            Difference::Match(2),
            Difference::Mismatch(b'T'),
            Difference::Match(3),
            Difference::Mismatch(b'T'),
            Difference::Match(3)
        ];
        check_difference(seq, &truth, name, false);
        let mut target_sequence = vec![b'-'; 150];
        target_sequence[111] = b'G';
        target_sequence[132] = b'C';
        target_sequence[139] = b'G';
        target_sequence[142] = b'G';
        target_sequence[146] = b'G';
        check_to_bytes(&truth, seq, &target_sequence, name);
    }
}

#[test]
fn difference_case() {
    {
        // Lower case mismatch.
        let seq = b"*ac";
        let truth = [ Difference::Mismatch(b'C') ];
        check_difference(seq, &truth, "lower case mismatch", false);
    }
    {
        // Lower case insertion.
        let seq = b"+acgt";
        let truth = [ Difference::Insertion(b"ACGT".to_vec()) ];
        check_difference(seq, &truth, "lower case insertion", false);
    }
}

fn invalid_difference(seq: &[u8], name: &str) {
    let result = Difference::parse(seq);
    assert!(result.is_err(), "Parsed an invalid difference string: {}", name);
}

#[test]
fn difference_invalid() {
    invalid_difference(b"ABC", "invalid operation at start");
    invalid_difference(b"+AC:", "empty operation at end");
    invalid_difference(b"~AC30AC", "intron / splice operation");
    invalid_difference(b":123A", "match length is not an integer");
    invalid_difference(b"+GATTACA:", "empty match length");
    invalid_difference(b"*ACT", "mismatch is too long");
    invalid_difference(b"*A", "mismatch is too short");
}

#[test]
fn difference_len() {
    {
        let diff = Difference::Match(42);
        assert_eq!(diff.len(), 42, "Wrong length for a match");
    }
    {
        let diff = Difference::Mismatch(b'A');
        assert_eq!(diff.len(), 1, "Wrong length for a mismatch");
    }
    {
        let diff = Difference::Insertion(b"GATTACA".to_vec());
        assert_eq!(diff.len(), 7, "Wrong length for an insertion");
    }
    {
        let diff = Difference::Deletion(42);
        assert_eq!(diff.len(), 42, "Wrong length for a deletion");
    }
}

#[test]
fn difference_normalize() {
    {
        let seq = b"=ACGT:4*AC*GT:2:3=ACT*AC-ACGT-ACGT+GATTACA+CAT";
        let truth = [
            Difference::Match(8),
            Difference::Mismatch(b'C'),
            Difference::Mismatch(b'T'),
            Difference::Match(8),
            Difference::Mismatch(b'C'),
            Difference::Deletion(8),
            Difference::Insertion(b"GATTACACAT".to_vec()),
        ];
        check_difference(seq, &truth, "normalized", true);
    }
    {
        let seq = b"=:0+-*AC:30=GATTA+-:0=";
        let truth = [
            Difference::Mismatch(b'C'),
            Difference::Match(35),
        ];
        check_difference(seq, &truth, "normalized with empty operations", true);
    }
}

}

//-----------------------------------------------------------------------------
