//! In-process external sort of a PAF by `Query_Name`, plus the read-ID set check.
//!
//! Sorting whole lines byte-lexicographically groups records by `Query_Name`
//! (the field-terminating `\t` = 0x09 sorts below any name byte) and orders the
//! groups by `Query_Name` — equivalent to `LC_ALL=C sort -t$'\t' -k1,1` for PAF,
//! which has no tabs/newlines inside names. Backed by `ext-sort`: it buffers up
//! to a byte limit in RAM, spills sorted runs to temp files, then k-way merges.

use std::cmp::Ordering;
use std::io::{BufRead, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};
use ext_sort::buffer::mem::MemoryLimitedBufferBuilder;
use ext_sort::{ExternalSorter, ExternalSorterBuilder};

use crate::io_utils::{open_input, open_output};

/// Parse a memory size like `512` (bytes), `500M`, `2G`, `1g` into bytes.
pub(crate) fn parse_mem(s: &str) -> Result<u64> {
    let s = s.trim();
    if s.is_empty() {
        bail!("empty --mem value");
    }
    let last = s.chars().last().unwrap();
    let (num, mult) = match last {
        'k' | 'K' => (&s[..s.len() - 1], 1024u64),
        'm' | 'M' => (&s[..s.len() - 1], 1024 * 1024),
        'g' | 'G' => (&s[..s.len() - 1], 1024 * 1024 * 1024),
        '0'..='9' => (s, 1u64),
        c => bail!("invalid --mem suffix '{c}' (use K/M/G, or plain bytes)"),
    };
    let n: u64 = num
        .trim()
        .parse()
        .with_context(|| format!("invalid --mem value '{s}'"))?;
    Ok(n * mult)
}

/// Sort the PAF at `input` by `Query_Name` (byte-lex, whole-line) into `output`
/// (`.gz` is gzip-compressed). Buffers up to `mem_bytes` in RAM, spilling sorted
/// runs to `tmp_dir`.
pub(crate) fn sort_paf_to_file(
    input: &str,
    output: &str,
    mem_bytes: u64,
    tmp_dir: &Path,
    threads: Option<usize>,
) -> Result<()> {
    let reader = open_input(input)?;

    let mut builder = ExternalSorterBuilder::new()
        .with_tmp_dir(tmp_dir)
        .with_buffer(MemoryLimitedBufferBuilder::new(mem_bytes));
    if let Some(t) = threads {
        builder = builder.with_threads_number(t);
    }
    let sorter: ExternalSorter<String, _, _> =
        builder.build().context("failed to build external sorter")?;

    let sorted = sorter
        .sort(reader.lines())
        .context("external sort failed")?;

    let mut out = open_output(Some(output))?;
    for item in sorted {
        let line = item.context("error while reading a line during external sort")?;
        writeln!(out, "{line}")?;
    }
    out.flush()?;
    Ok(())
}

// ── Read-ID set check ─────────────────────────────────────────────────────────

/// Counts + a few examples from comparing the `Query_Name` sets of two PAFs.
pub(crate) struct IdSetCheck {
    pub shared: u64,
    pub only_a: u64,
    pub only_b: u64,
    pub examples_a: Vec<String>,
    pub examples_b: Vec<String>,
}

/// Streams the *distinct* `Query_Name`s of a sorted PAF in ascending order.
/// Constant memory (one buffered line + the last name returned).
struct DistinctNames {
    reader: Box<dyn BufRead>,
    buf: String,
    prev: Option<String>,
}

impl DistinctNames {
    fn new(path: &str) -> Result<Self> {
        Ok(Self {
            reader: open_input(path)?,
            buf: String::new(),
            prev: None,
        })
    }

    fn next_name(&mut self) -> Result<Option<String>> {
        loop {
            self.buf.clear();
            let n = self.reader.read_line(&mut self.buf)?;
            if n == 0 {
                return Ok(None);
            }
            let trimmed = self.buf.trim_end_matches(['\n', '\r']);
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let name = trimmed.split('\t').next().unwrap_or("");
            // Skip consecutive duplicates (the file is sorted, so equal names
            // are contiguous).
            if self.prev.as_deref() == Some(name) {
                continue;
            }
            let owned = name.to_owned();
            self.prev = Some(owned.clone());
            return Ok(Some(owned));
        }
    }
}

/// Merge the distinct-name streams of two **sorted** PAFs to determine whether
/// they carry the same `Query_Name` set. O(1) memory.
pub(crate) fn read_id_set_check(a: &str, b: &str, max_examples: usize) -> Result<IdSetCheck> {
    let mut da = DistinctNames::new(a)?;
    let mut db = DistinctNames::new(b)?;
    let mut na = da.next_name()?;
    let mut nb = db.next_name()?;

    let mut c = IdSetCheck {
        shared: 0,
        only_a: 0,
        only_b: 0,
        examples_a: Vec::new(),
        examples_b: Vec::new(),
    };

    loop {
        match (na.as_deref(), nb.as_deref()) {
            (Some(x), Some(y)) => match x.as_bytes().cmp(y.as_bytes()) {
                Ordering::Equal => {
                    c.shared += 1;
                    na = da.next_name()?;
                    nb = db.next_name()?;
                }
                Ordering::Less => {
                    c.only_a += 1;
                    if c.examples_a.len() < max_examples {
                        c.examples_a.push(x.to_owned());
                    }
                    na = da.next_name()?;
                }
                Ordering::Greater => {
                    c.only_b += 1;
                    if c.examples_b.len() < max_examples {
                        c.examples_b.push(y.to_owned());
                    }
                    nb = db.next_name()?;
                }
            },
            (Some(x), None) => {
                c.only_a += 1;
                if c.examples_a.len() < max_examples {
                    c.examples_a.push(x.to_owned());
                }
                na = da.next_name()?;
            }
            (None, Some(y)) => {
                c.only_b += 1;
                if c.examples_b.len() < max_examples {
                    c.examples_b.push(y.to_owned());
                }
                nb = db.next_name()?;
            }
            (None, None) => break,
        }
    }
    Ok(c)
}
