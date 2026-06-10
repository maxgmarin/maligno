//! Shared PAF → per-read group reader.
//!
//! Turns a PAF `BufRead` into a sequence of per-read groups, where each group
//! is the contiguous run of alignment rows sharing one `Query_Name`. Each row
//! is materialized into the exact same `AlnRow` form `readinfo` works with, by
//! routing every PAF record through the unchanged conversion + serialization
//! path:
//!
//! ```text
//! PAF line ─parse_line─▶ PafRecord ─AlnInfo::from_paf─▶ AlnInfo
//!          ─write_row(&mut Vec<u8>)─▶ alninfo line bytes   (identical to paf2alninfo)
//!          ─parse_aln_row─▶ AlnRow                          (identical to readinfo input)
//! ```
//!
//! Because the bytes crossing the boundary are byte-for-byte what `paf2alninfo`
//! would have written to disk, downstream `collapse_group` produces output
//! identical to the discrete `paf2alninfo | readinfo` pipeline — without ever
//! materializing the intermediate alninfo file.
//!
//! Used by the fused `paf2readinfo` and `pafcompare` subcommands.

use std::io::BufRead;

use anyhow::Result;

use crate::paf::parse_line;
use crate::readinfo::{parse_aln_row, AlnRow};
use crate::record::AlnInfo;

/// Streaming reader that yields one per-read group (`Vec<AlnRow>`) at a time.
///
/// Requires the input PAF to be **grouped** by `Query_Name` (each read's
/// alignments contiguous). Constant memory: only the current group plus a
/// one-row lookahead are held at any time.
pub(crate) struct PafGroups<R: BufRead> {
    reader: R,
    lineno: u64,
    line: String,            // reusable line read buffer
    serialize_buf: Vec<u8>,  // reusable AlnInfo::write_row target
    pending: Option<AlnRow>, // first row of the next group (lookahead)
    eof: bool,
    warn_unsorted: bool,     // emit the readinfo-style lex-decrease warning once
    warned: bool,
    last_group_name: Option<String>, // last yielded group's Query_Name (for the warning)
}

impl<R: BufRead> PafGroups<R> {
    pub(crate) fn new(reader: R, warn_unsorted: bool) -> Self {
        PafGroups {
            reader,
            lineno: 0,
            line: String::with_capacity(4096),
            serialize_buf: Vec::with_capacity(512),
            pending: None,
            eof: false,
            warn_unsorted,
            warned: false,
            last_group_name: None,
        }
    }

    /// Read the next PAF data line and convert it to an `AlnRow`.
    /// Returns `Ok(None)` at EOF. Malformed lines emit a WARNING and are skipped
    /// (same behavior as `paf2alninfo`).
    fn next_row(&mut self) -> Result<Option<AlnRow>> {
        loop {
            if self.eof {
                return Ok(None);
            }
            self.line.clear();
            let n = self.reader.read_line(&mut self.line)?;
            if n == 0 {
                self.eof = true;
                return Ok(None);
            }
            self.lineno += 1;
            let trimmed = self
                .line
                .trim_end_matches(|c| c == '\n' || c == '\r');
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let rec = match parse_line(trimmed, self.lineno) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("WARNING: {e}");
                    continue;
                }
            };
            // PAF record → AlnInfo → exact alninfo TSV bytes → AlnRow.
            self.serialize_buf.clear();
            AlnInfo::from_paf(&rec).write_row(&mut self.serialize_buf)?;
            // write_row appends a trailing newline; trim it before parsing.
            let s = std::str::from_utf8(&self.serialize_buf)
                .expect("alninfo serialization is valid UTF-8")
                .trim_end_matches(|c| c == '\n' || c == '\r');
            match parse_aln_row(s) {
                Some(row) => return Ok(Some(row)),
                None => continue, // defensive; should not happen for our own output
            }
        }
    }

    /// Return the next contiguous-`Query_Name` group, or `None` at EOF.
    pub(crate) fn next_group(&mut self) -> Result<Option<Vec<AlnRow>>> {
        // Seed the group with the pending lookahead row, or the next read row.
        let first = match self.pending.take() {
            Some(r) => r,
            None => match self.next_row()? {
                Some(r) => r,
                None => return Ok(None),
            },
        };

        let group_name = first.query_name().to_owned();

        // One-shot byte-lex-decrease warning across group boundaries.
        if self.warn_unsorted && !self.warned {
            if let Some(prev) = &self.last_group_name {
                if group_name.as_bytes() < prev.as_bytes() {
                    eprintln!(
                        "WARNING: input PAF is not byte-lex sorted by Query_Name \
                         (saw {prev:?} then {group_name:?}). Grouping by *contiguous* \
                         Query_Name runs still works, but the downstream `compare` \
                         step requires byte-lex sort. To pre-sort:\n\
                         \n\
                         \tLC_ALL=C sort -t$'\\t' -k1,1 in.paf > sorted.paf\n"
                    );
                    self.warned = true;
                }
            }
        }
        self.last_group_name = Some(group_name.clone());

        let mut group: Vec<AlnRow> = vec![first];

        loop {
            match self.next_row()? {
                Some(row) => {
                    if row.query_name() == group_name {
                        group.push(row);
                    } else {
                        // Boundary: stash for the next call and return this group.
                        self.pending = Some(row);
                        break;
                    }
                }
                None => break, // EOF ends the final group
            }
        }

        Ok(Some(group))
    }
}
