use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

/// Escape TSV field value: replace newlines, carriage returns, and tabs
/// to prevent breaking TSV format when fields contain special characters.
pub fn escape_tsv_field(field: &str) -> String {
    field
        .replace('\\', "\\\\")  // Backslash first to avoid double-escaping
        .replace('\n', "\\n")   // Newline
        .replace('\r', "\\r")   // Carriage return
        .replace('\t', "\\t")   // Tab
}

/// Open an input reader.
/// - `"-"` → stdin
/// - path ending in `.gz` → gzip-compressed file
/// - otherwise → plain file
pub fn open_input(path: &str) -> Result<Box<dyn BufRead>> {
    match path {
        "-" => Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            io::stdin().lock(),
        ))),
        p if p.ends_with(".gz") => {
            let file = File::open(p).with_context(|| format!("cannot open '{p}'"))?;
            let decoder = GzDecoder::new(file);
            Ok(Box::new(BufReader::with_capacity(1 << 20, decoder)))
        }
        p => {
            let file = File::open(p).with_context(|| format!("cannot open '{p}'"))?;
            Ok(Box::new(BufReader::with_capacity(1 << 20, file)))
        }
    }
}

/// Open an output writer.
/// - `None` or `"-"` → stdout
/// - path ending in `.gz` → gzip-compressed file
/// - otherwise → plain file
pub fn open_output(path: Option<&str>) -> Result<Box<dyn Write>> {
    match path {
        None | Some("-") => Ok(Box::new(BufWriter::with_capacity(
            1 << 20,
            io::stdout().lock(),
        ))),
        Some(p) if p.ends_with(".gz") => {
            let file =
                File::create(p).with_context(|| format!("cannot create '{p}'"))?;
            let encoder = GzEncoder::new(file, Compression::default());
            Ok(Box::new(BufWriter::with_capacity(1 << 20, encoder)))
        }
        Some(p) => {
            let file =
                File::create(p).with_context(|| format!("cannot create '{p}'"))?;
            Ok(Box::new(BufWriter::with_capacity(1 << 20, file)))
        }
    }
}

/// Format a float to match Python's repr(float):
/// integer-valued floats get ".0" appended; NaN → "NaN".
pub fn fmt_float(v: f64) -> String {
    if v.is_nan() {
        return "NaN".to_owned();
    }
    if v.is_infinite() {
        return if v > 0.0 {
            "inf".to_owned()
        } else {
            "-inf".to_owned()
        };
    }
    let s = format!("{v}");
    // If Rust's Display produced no decimal point and no exponent marker,
    // the value is integer-valued: append ".0".
    if s.bytes().all(|b| b == b'-' || b.is_ascii_digit()) {
        format!("{s}.0")
    } else {
        s
    }
}
