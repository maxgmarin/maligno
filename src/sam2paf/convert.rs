use std::collections::HashMap;
use std::io::{BufRead, Write};

use anyhow::{Context, Result};

use super::cigar::{build_merged_cigar, write_cigar_no_clips, CigarStats};
use super::cs_generator::generate_cs;

/// Runtime options forwarded from the CLI.
#[derive(Debug, Clone, Default)]
pub struct Options {
    /// Skip secondary alignments (flag 0x100).
    pub pri_only: bool,
    /// Skip secondary AND supplementary alignments (flags 0x100 | 0x800).
    pub pri_pri_only: bool,
    /// Emit long-form cs (`=ACGT`) instead of short-form (`:N`).
    pub long_cs: bool,
    /// Emit placeholder PAF records for unmapped reads.
    pub convert_unaligned: bool,
}

/// Read a SAM stream line by line and write PAF records to `writer`.
pub fn convert<R: BufRead, W: Write>(
    mut reader: R,
    writer: &mut W,
    opts: &Options,
) -> Result<()> {
    let mut ctg_len: HashMap<String, u64> = HashMap::new();
    let mut line = String::with_capacity(4096);
    let mut lineno = 0u64;

    while reader.read_line(&mut line)? > 0 {
        lineno += 1;
        // Strip newline; handle both \n and \r\n.
        let record = line.trim_end_matches(|c| c == '\n' || c == '\r');

        if record.starts_with('@') {
            // SAM header line.
            if record.starts_with("@SQ\t") {
                parse_sq_line(record, &mut ctg_len);
            }
            line.clear();
            continue;
        }

        if let Err(e) = process_record(record, &ctg_len, opts, lineno, writer) {
            // Match JS behaviour: print warnings to stderr but keep going for
            // non-fatal errors.  Fatal errors (unknown contig, inconsistent MD)
            // propagate up.
            eprintln!("WARNING at line {lineno}: {e}");
        }
        line.clear();
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// @SQ header parsing
// ---------------------------------------------------------------------------

fn parse_sq_line(line: &str, ctg_len: &mut HashMap<String, u64>) {
    let bytes = line.as_bytes();
    let sn = find_tag_value(bytes, b"SN:");
    let ln = find_tag_value(bytes, b"LN:");
    if let (Some(name), Some(len_str)) = (sn, ln) {
        if let Ok(l) = len_str.parse::<u64>() {
            ctg_len.insert(name.to_owned(), l);
        }
    }
}

/// Scan `bytes` for `key` (e.g. `b"SN:"`) and return everything up to the
/// next tab or end-of-string.
fn find_tag_value<'a>(bytes: &'a [u8], key: &[u8]) -> Option<&'a str> {
    let key_len = key.len();
    bytes.windows(key_len).position(|w| w == key).map(|pos| {
        let start = pos + key_len;
        let end = bytes[start..].iter().position(|&b| b == b'\t').map_or(bytes.len(), |i| start + i);
        // SAFETY: SAM is ASCII; slicing on tab boundaries is safe.
        unsafe { std::str::from_utf8_unchecked(&bytes[start..end]) }
    })
}

// ---------------------------------------------------------------------------
// Per-record processing
// ---------------------------------------------------------------------------

/// Process one SAM alignment record and write zero or one PAF lines.
fn process_record<W: Write>(
    line: &str,
    ctg_len: &HashMap<String, u64>,
    opts: &Options,
    lineno: u64,
    out: &mut W,
) -> Result<()> {
    // Split the mandatory SAM fields.  `splitn(12, '\t')` gives:
    //   indices 0..=10 → the 11 mandatory fields
    //   index  11      → everything after QUAL (the optional tags joined by \t)
    let mut parts = line.splitn(12, '\t');
    macro_rules! next_field {
        ($name:expr) => {
            parts.next().with_context(|| {
                format!("line {lineno}: missing SAM field '{}'", $name)
            })?
        };
    }

    let qname_raw   = next_field!("QNAME");
    let flag_str    = next_field!("FLAG");
    let rname       = next_field!("RNAME");
    let pos_str     = next_field!("POS");
    let mapq_str    = next_field!("MAPQ");
    let cigar_str   = next_field!("CIGAR");
    let _rnext      = next_field!("RNEXT");
    let _pnext      = next_field!("PNEXT");
    let _tlen_sam   = next_field!("TLEN");
    let seq         = next_field!("SEQ");
    let qual        = next_field!("QUAL");
    let tags_raw    = parts.next().unwrap_or("");

    let flag: u32 = flag_str.parse().with_context(|| format!("line {lineno}: invalid FLAG"))?;

    // SEQ / QUAL length consistency (only when both are present).
    if seq != "*" && qual != "*" && seq.len() != qual.len() {
        anyhow::bail!(
            "line {lineno}: inconsistent SEQ and QUAL lengths ({} != {})",
            seq.len(), qual.len()
        );
    }

    // Unmapped record.
    if rname == "*" || (flag & 4 != 0) || cigar_str == "*" {
        if opts.convert_unaligned {
            let qlen = seq.len();
            writeln!(out, "{}\t{}\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0", qname_raw, qlen)?;
        }
        return Ok(());
    }

    // Primary / supplementary filter.
    if opts.pri_only && (flag & 0x100 != 0) { return Ok(()); }
    if opts.pri_pri_only && (flag & 0x900 != 0) { return Ok(()); }

    let tlen = ctg_len.get(rname).with_context(|| {
        format!("line {lineno}: contig '{rname}' not found in SAM header @SQ lines")
    })?;

    // Parse optional tags in one pass.
    let tags = parse_tags(tags_raw);

    // Parse CIGAR statistics in one pass.
    let cs_stats = CigarStats::parse(cigar_str);
    if cs_stats.n_ops > 65535 {
        eprintln!("WARNING at line {lineno}: {} CIGAR operations", cs_stats.n_ops);
    }

    // Reference coordinates (0-based half-open).
    let pos1: u64 = pos_str.parse().with_context(|| format!("line {lineno}: invalid POS"))?;
    let ts: u64 = pos1 - 1;
    let te: u64 = ts + cs_stats.ref_len() as u64;

    if te > *tlen {
        return Err(anyhow::anyhow!(
            "line {lineno}: alignment end ({te}) > contig length ({tlen}); skipped"
        ));
    }

    // SEQ / CIGAR length consistency.
    let ql = cs_stats.query_consumed();
    if seq != "*" && seq.len() as u32 != ql {
        return Err(anyhow::anyhow!(
            "line {lineno}: SEQ length ({}) inconsistent with CIGAR ({ql}); skipped",
            seq.len()
        ));
    }

    // Derive mismatches and calibrate NM.
    let (nm, mm) = calibrate_nm(&cs_stats, tags.nm, lineno);

    // Matching bases and block length.
    let mlen = cs_stats.m.saturating_sub(mm);
    let blen = cs_stats.block_len();

    // Query name: append /1 or /2 for paired reads (rare; avoids alloc otherwise).
    let rev = flag & 0x10 != 0;
    let qname_buf: String; // only used when paired
    let qname: &str = if flag & 1 != 0 {
        qname_buf = build_qname(qname_raw, flag);
        &qname_buf
    } else {
        qname_raw
    };

    // Query coordinates (always on the forward strand of the query).
    let qlen = cs_stats.query_len();
    let (qs, qe) = if rev {
        (cs_stats.clip_trail, qlen - cs_stats.clip_lead)
    } else {
        (cs_stats.clip_lead, qlen - cs_stats.clip_trail)
    };

    // Generate cs tag if not already present as a cs:Z: SAM tag.
    let cs_out: Option<String> = if let Some(cs_str) = tags.cs {
        Some(cs_str.to_owned())
    } else if tags.md.is_some() && seq != "*" {
        let merged = build_merged_cigar(cigar_str);
        match generate_cs(&merged, tags.md.unwrap(), seq, opts.long_cs, lineno) {
            Ok(s) => Some(s),
            Err(e) => { eprintln!("WARNING: {e}"); None }
        }
    } else {
        None
    };

    // Alignment type: 'S' = secondary, 'P' = primary/supplementary.
    let aln_type = if flag & 0x100 != 0 { b'S' } else { b'P' };
    let strand   = if rev { b'-' } else { b'+' };

    // Write PAF line: all fields tab-separated, newline at end.
    write!(out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\ttp:A:{}",
        qname, qlen, qs, qe, strand as char,
        rname, tlen, ts, te,
        mlen, blen, mapq_str, aln_type as char
    )?;
    if let Some(nm_val) = nm {
        write!(out, "\tNM:i:{nm_val}\tmm:i:{mm}")?;
    }
    if let Some(a) = tags.r#as { write!(out, "\tAS:i:{a}")?; }
    if let Some(m) = tags.ms  { write!(out, "\tms:i:{m}")?; }
    write!(out, "\tgn:i:{}\tgo:i:{}\tcg:Z:",
           cs_stats.i_bases + cs_stats.d_bases,
           cs_stats.i_count + cs_stats.d_count)?;
    // Write CIGAR without clip ops directly — avoids one String allocation.
    write_cigar_no_clips(out, cigar_str)?;
    if let Some(ref cs) = cs_out { write!(out, "\tcs:Z:{cs}")?; }
    writeln!(out)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Parsed values of the optional SAM tags we care about.
#[derive(Default)]
struct SamTags<'a> {
    nm: Option<u32>,
    _nn: u32,
    md: Option<&'a str>,
    cs: Option<&'a str>,
    r#as: Option<i64>,
    ms: Option<i64>,
}

/// Parse optional SAM tags from the raw tag substring (already past QUAL).
fn parse_tags(raw: &str) -> SamTags<'_> {
    let mut t = SamTags::default();
    // Each tag is `TT:T:VALUE` separated by tabs.
    for field in raw.split('\t') {
        let b = field.as_bytes();
        if b.len() < 6 { continue; }
        // Check the type byte (b[3]) to branch quickly.
        match &b[..5] {
            b"NM:i:" => { t.nm  = parse_int(&b[5..]).map(|v| v as u32); }
            b"nn:i:" => { t._nn = parse_int(&b[5..]).unwrap_or(0) as u32; }
            b"MD:Z:" => { t.md  = Some(&field[5..]); }
            b"cs:Z:" => { t.cs  = Some(&field[5..]); }
            b"AS:i:" => { t.r#as = parse_int(&b[5..]); }
            b"ms:i:" => { t.ms  = parse_int(&b[5..]); }
            _ => {}
        }
    }
    t
}

/// Parse an integer from a byte slice (no allocations).
#[inline]
fn parse_int(b: &[u8]) -> Option<i64> {
    // SAFETY: SAM integers are ASCII digits (possibly with leading '-').
    unsafe { std::str::from_utf8_unchecked(b) }.parse().ok()
}

/// Build the query name, appending `/1` or `/2` for paired-end reads.
fn build_qname(raw: &str, flag: u32) -> String {
    if flag & 1 != 0 {
        if flag & 0x40 != 0 { return format!("{raw}/1"); }
        if flag & 0x80 != 0 { return format!("{raw}/2"); }
    }
    raw.to_owned()
}

/// Calibrate the mismatch count and NM value from CIGAR stats.
///
/// Returns `(nm, mm)` where `nm` is `None` when there is truly no information
/// (standard CIGAR, NM absent).
fn calibrate_nm(cs: &CigarStats, nm_tag: Option<u32>, lineno: u64) -> (Option<u32>, u32) {
    let gap_total = cs.i_bases + cs.d_bases;

    if cs.have_ext && !cs.have_m {
        // Extended CIGAR (= and X only): derive everything from ops.
        let computed = gap_total + cs.mm_ext;
        if let Some(nm) = nm_tag {
            if nm != computed {
                eprintln!(
                    "WARNING at line {lineno}: NM ({nm}) != sum of gaps+mismatches ({computed})"
                );
            }
        }
        (Some(computed), cs.mm_ext)
    } else if let Some(nm) = nm_tag {
        // Standard CIGAR + NM present.
        let nm = if nm < gap_total {
            eprintln!(
                "WARNING at line {lineno}: NM ({nm}) < total gaps ({gap_total})"
            );
            gap_total
        } else {
            nm
        };
        let mm = nm - gap_total;
        (Some(nm), mm)
    } else {
        // Standard CIGAR, no NM: assume no mismatches.
        eprintln!("WARNING at line {lineno}: no NM tag; assuming zero mismatches");
        (None, 0)
    }
}
