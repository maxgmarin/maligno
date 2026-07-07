use anyhow::{bail, Context, Result};

/// The 12 mandatory PAF fields plus the three optional tags we need.
#[derive(Debug)]
pub struct PafRecord<'a> {
    pub query_name:           &'a str,
    pub query_len:            u64,
    pub query_start:          u64,
    pub query_end:            u64,
    pub strand:               char,
    pub target_name:          &'a str,
    pub target_len:           u64,
    pub target_start:         u64,
    pub target_end:           u64,
    pub num_residue_matches:  u64,
    pub aln_block_len:        u64,
    pub mapq:                 u32,
    /// SAM `ms:i` tag (mate / second-best alignment score).
    pub ms:                   i64,
    /// SAM `AS:i` tag (alignment score).
    pub aln_score:            i64,
    /// cs tag string (borrowed from the input line).
    pub cs:                   &'a str,
    /// True when Target_Name == "*" (unmapped record).
    pub is_unmapped:          bool,
}

/// Parse one PAF line into a `PafRecord`.
///
/// The returned record borrows string slices directly from `line`, so `line`
/// must outlive the record.
pub fn parse_line(line: &str, lineno: u64) -> Result<PafRecord<'_>> {
    // Split into at most 13 parts: the 12 mandatory fields + the rest (tags).
    let mut it = line.splitn(13, '\t');

    macro_rules! nxt {
        ($name:expr) => {
            it.next().ok_or_else(|| anyhow::anyhow!("line {lineno}: missing PAF field '{}'", $name))?
        };
    }
    macro_rules! parse_int {
        ($s:expr, $name:expr) => {{
            let v: &str = $s;
            v.parse::<u64>().map_err(|_| anyhow::anyhow!("line {lineno}: cannot parse field '{}' as integer (got {:?})", $name, v))?
        }};
    }

    let query_name  = nxt!("QNAME");
    let query_len   = parse_int!(nxt!("QLEN"),  "QLEN");
    let query_start = parse_int!(nxt!("QSTART"), "QSTART");
    let query_end   = parse_int!(nxt!("QEND"),   "QEND");
    let strand_str  = nxt!("STRAND");
    let target_name = nxt!("TNAME");
    let target_len  = parse_int!(nxt!("TLEN"),   "TLEN");
    let target_start= parse_int!(nxt!("TSTART"), "TSTART");
    let target_end  = parse_int!(nxt!("TEND"),   "TEND");
    let num_res     = parse_int!(nxt!("NRES"),   "NRES");
    let aln_blk     = parse_int!(nxt!("BLKLEN"), "BLKLEN");
    let mapq        = nxt!("MAPQ").parse::<u32>()
                        .with_context(|| format!("line {lineno}: cannot parse MAPQ"))?;
    let tags_raw    = it.next().unwrap_or("");

    let strand = match strand_str {
        "+" => '+',
        "-" => '-',
        "*" => '*',
        other => bail!("line {lineno}: unexpected strand '{other}'"),
    };

    let is_unmapped = target_name == "*" || strand == '*';

    // Parse optional tags: scan tab-separated fields for ms:i, AS:i, cs:Z.
    let (ms, aln_score, cs) = parse_tags(tags_raw, is_unmapped);

    Ok(PafRecord {
        query_name, query_len, query_start, query_end, strand,
        target_name, target_len, target_start, target_end,
        num_residue_matches: num_res, aln_block_len: aln_blk, mapq,
        ms, aln_score, cs,
        is_unmapped,
    })
}

/// Scan the optional-tag substring for `ms:i:`, `AS:i:`, and `cs:Z:`.
/// Returns `(ms, AS, cs_slice)`.
///
/// For unmapped records it will return (0, 0, "").
fn parse_tags(raw: &str, is_unmapped: bool) -> (i64, i64, &str) {
    if is_unmapped {
        return (0, 0, "");
    }

    let mut ms: i64  = 0;
    let mut aln: i64 = 0;
    let mut cs: &str = "";

    for field in raw.split('\t') {
        let b = field.as_bytes();
        if b.len() < 6 { continue; }
        match &b[..5] {
            b"ms:i:" => { ms  = parse_signed_int(&b[5..]); }
            b"AS:i:" => { aln = parse_signed_int(&b[5..]); }
            b"cs:Z:" => { cs  = &field[5..]; }
            _ => {}
        }
    }

    (ms, aln, cs)
}

#[inline]
fn parse_signed_int(b: &[u8]) -> i64 {
    // SAFETY: SAM/PAF integer tags are ASCII; slice is valid UTF-8.
    unsafe { std::str::from_utf8_unchecked(b) }
        .parse()
        .unwrap_or(0)
}
