//! Utility functions and structures.

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::io::{BufRead, BufReader, Read};

use flate2::read::MultiGzDecoder;

//-----------------------------------------------------------------------------

/// Returns the full file name for a specific test file.
pub fn get_test_data(filename: &'static str) -> PathBuf {
    let mut buf = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    buf.push("test-data");
    buf.push(filename);
    buf
}

//-----------------------------------------------------------------------------

// Utilities for working with files.

const SIZE_UNITS: [(f64, &str); 6] = [
    (1.0, "B"),
    (1024.0, "KiB"),
    (1024.0 * 1024.0, "MiB"),
    (1024.0 * 1024.0 * 1024.0, "GiB"),
    (1024.0 * 1024.0 * 1024.0 * 1024.0, "TiB"),
    (1024.0 * 1024.0 * 1024.0 * 1024.0 * 1024.0, "PiB"),
];

/// Returns a human-readable representation of the given number of bytes.
pub fn human_readable_size(bytes: usize) -> String {
    let mut unit = 0;
    let value = bytes as f64;
    while unit + 1 < SIZE_UNITS.len() && value >= SIZE_UNITS[unit + 1].0 {
        unit += 1;
    }
    format!("{:.3} {}", value / SIZE_UNITS[unit].0, SIZE_UNITS[unit].1)
}

/// Returns a human-readable size of the file.
pub fn file_size<P: AsRef<Path>>(filename: P) -> Option<String> {
    let metadata = fs::metadata(filename).map_err(|x| x.to_string());
    if metadata.is_err() {
        return None;
    }
    Some(human_readable_size(metadata.unwrap().len() as usize))
}

/// Returns `true` if the file exists.
pub fn file_exists<P: AsRef<Path>>(filename: P) -> bool {
    fs::metadata(filename).is_ok()
}

/// Returns `true` if the file appears to be gzip-compressed.
pub fn is_gzipped<P: AsRef<Path>>(filename: P) -> bool {
    let file = File::open(filename).ok();
    if file.is_none() {
        return false;
    }
    let mut reader = BufReader::new(file.unwrap());
    let mut magic = [0; 2];
    let len = reader.read(&mut magic).ok();
    len == Some(2) && magic == [0x1F, 0x8B]
}

/// Returns a buffered reader for the file, which may be gzip-compressed.
pub fn open_file<P: AsRef<Path>>(filename: P) -> Result<Box<dyn BufRead>, String> {
    let file = File::open(&filename).map_err(|x| x.to_string())?;
    let inner = BufReader::new(file);
    if is_gzipped(&filename) {
        let inner = MultiGzDecoder::new(inner);
        Ok(Box::new(BufReader::new(inner)))
    } else {
        Ok(Box::new(inner))
    }
}

//-----------------------------------------------------------------------------

// Sequence encoding and decoding.

const DECODE: [u8; 6] = [0, b'A', b'C', b'G', b'T', b'N'];

/// Decodes a sequence encoded with [`encode_sequence`].
pub fn decode_sequence(encoded: &[u8]) -> Vec<u8> {
    let capacity = if encoded.is_empty() { 0 } else { 3 * encoded.len() - 1 };
    let mut result = Vec::with_capacity(capacity);

    for byte in encoded {
        let mut value = *byte as usize;
        for _ in 0..3 {
            let decoded = DECODE[value % DECODE.len()];
            if decoded == 0 {
                return result;
            }
            value /= DECODE.len();
            result.push(decoded);
        }
    }

    result
}

const fn generate_encoding() -> [u8; 256] {
    let mut result = [5; 256];
    result[b'a' as usize] = 1; result[b'A' as usize] = 1;
    result[b'c' as usize] = 2; result[b'C' as usize] = 2;
    result[b'g' as usize] = 3; result[b'G' as usize] = 3;
    result[b't' as usize] = 4; result[b'T' as usize] = 4;
    result
}

const ENCODE: [u8; 256] = generate_encoding();

// TODO: If the length is divisible by 3, we do not need the null terminator.
/// Encodes a DNA sequence into a byte array, storing three bases in a byte.
///
/// Values outside `acgtACGT` are encoded as `N`.
/// The last encoded symbol is a special 0 character in order to preserve the length.
/// See [`decode_sequence`] for decoding.
pub fn encode_sequence(sequence: &[u8]) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::with_capacity(sequence.len() / 3 + 1);

    let mut offset = 0;
    while offset + 3 <= sequence.len() {
        let byte = ENCODE[sequence[offset] as usize] +
            6 * ENCODE[sequence[offset + 1] as usize] +
            36 * ENCODE[sequence[offset + 2] as usize];
        result.push(byte);
        offset += 3;
    }
    let byte = match sequence.len() - offset {
        0 => 0,
        1 => ENCODE[sequence[offset] as usize],
        _ => ENCODE[sequence[offset] as usize] + 6 * ENCODE[sequence[offset + 1] as usize],
    };
    result.push(byte);

    result
}

//-----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequence_encoding() {
        let full_sequence = b"GATTACACACCAGATNNNNNACATTGAACCTTACACAGTCTGAC";
        for i in 0..full_sequence.len() {
            let sequence = &full_sequence[0..i];
            let encoded = encode_sequence(sequence);
            let decoded = decode_sequence(&encoded);
            assert_eq!(decoded, sequence, "Wrong sequence encoding for length {}", i);
        }
    }
}

//-----------------------------------------------------------------------------