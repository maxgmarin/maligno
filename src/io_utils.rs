use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
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
/// - path ending in `.gz` → gzip-compressed file (all members are read, so
///   multi-member files from bgzip, pigz or `cat a.gz b.gz` are read in full)
/// - otherwise → plain file
pub fn open_input(path: &str) -> Result<Box<dyn BufRead>> {
    match path {
        "-" => Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            io::stdin().lock(),
        ))),
        p if p.ends_with(".gz") => {
            let file = File::open(p).with_context(|| format!("cannot open '{p}'"))?;
            let decoder = MultiGzDecoder::new(file);
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn tmp_path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("maligno_io_utils_{}_{name}", std::process::id()))
    }

    fn gzip_member(data: &[u8]) -> Vec<u8> {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(data).unwrap();
        enc.finish().unwrap()
    }

    /// PAF-like text larger than one 64 KB BGZF block.
    fn paf_text() -> Vec<u8> {
        let mut s = String::new();
        for i in 0..5000 {
            s.push_str(&format!(
                "read{i:05}\t1000\t0\t1000\t+\tchr22\t50818468\t{i}\t{}\t1000\t1000\t60\n",
                i + 1000
            ));
        }
        assert!(s.len() > 3 * 65536);
        s.into_bytes()
    }

    fn read_all(path: &std::path::Path) -> Vec<u8> {
        let mut out = Vec::new();
        open_input(path.to_str().unwrap())
            .unwrap()
            .read_to_end(&mut out)
            .unwrap();
        out
    }

    #[test]
    fn reads_all_gzip_members() {
        // BGZF layout: one gzip member per <=64 KB chunk, then an empty EOF member.
        let text = paf_text();
        let mut gz = Vec::new();
        for chunk in text.chunks(65280) {
            gz.extend(gzip_member(chunk));
        }
        gz.extend(gzip_member(b""));
        let path = tmp_path("bgzf_like.paf.gz");
        std::fs::write(&path, &gz).unwrap();
        let got = read_all(&path);
        std::fs::remove_file(&path).ok();
        assert_eq!(got.len(), text.len());
        assert_eq!(got, text);
    }

    #[test]
    fn reads_concatenated_gz_files() {
        let a = b"readA\tpart one\n".to_vec();
        let b = b"readB\tpart two\n".to_vec();
        let mut gz = gzip_member(&a);
        gz.extend(gzip_member(&b));
        let path = tmp_path("cat.paf.gz");
        std::fs::write(&path, &gz).unwrap();
        let got = read_all(&path);
        std::fs::remove_file(&path).ok();
        assert_eq!(got, [a, b].concat());
    }

    #[test]
    fn single_member_gzip_unchanged() {
        let text = paf_text();
        let gz_path = tmp_path("single.paf.gz");
        let plain_path = tmp_path("single.paf");
        std::fs::write(&gz_path, gzip_member(&text)).unwrap();
        std::fs::write(&plain_path, &text).unwrap();
        let (got_gz, got_plain) = (read_all(&gz_path), read_all(&plain_path));
        std::fs::remove_file(&gz_path).ok();
        std::fs::remove_file(&plain_path).ok();
        assert_eq!(got_gz, text);
        assert_eq!(got_plain, text);
    }
}
