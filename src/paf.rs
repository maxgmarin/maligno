use anyhow::{bail, Context, Result};

/// The 12 mandatory PAF fields plus the optional tags we need.
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
    /// SAM `ms:i` tag (mapping Score: the maximum-scoring sub-segment within a given alignment).
    pub ms:                   i64,
    /// SAM `AS:i` tag (alignment score).
    pub aln_score:            i64,
    /// cs tag string (borrowed from the input line).
    pub cs:                   &'a str,
    /// SAM/PAF `tp:A` tag (alignment type: 'P' primary, 'S' secondary,
    /// 'I' inversion of the primary). '*' when absent or unmapped — the
    /// same sentinel this codebase already uses for an unknown `strand`.
    pub tp_tag:               char,
    /// Strand of the transcript relative to the **read** (`+` sense, `-`
    /// antisense, `*` unknown), from `ts:A` or `XS:A`; see [`transcript_strand`].
    pub transcript_strand:    char,
    /// True when target_name == "*" (unmapped record).
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

    let tags = parse_tags(tags_raw, is_unmapped);
    let transcript_strand = if is_unmapped { '*' } else { transcript_strand(strand, tags.ts, tags.xs) };

    Ok(PafRecord {
        query_name, query_len, query_start, query_end, strand,
        target_name, target_len, target_start, target_end,
        num_residue_matches: num_res, aln_block_len: aln_blk, mapq,
        ms: tags.ms, aln_score: tags.aln_score, cs: tags.cs, tp_tag: tags.tp_tag,
        transcript_strand,
        is_unmapped,
    })
}

/// The transcript's strand relative to the **read**: `+` when the read is in
/// the transcript's sense orientation, `-` when antisense, `*` when unknown.
/// Being read-relative, it does not change with the reference the read is
/// aligned to (e.g. a reverse-complemented assembly flips `strand`, not this).
///
/// - `ts:A` (minimap2 splice modes) is already read-relative: used as written.
/// - `XS:A` (STAR `--outSAMstrandField intronMotif`, HISAT2; Cufflinks
///   convention) is relative to the reference: `+` when it matches the
///   alignment `strand`, `-` when it doesn't.
/// - `ts:A` wins when both are present. Any other value, or neither tag, is `*`.
pub fn transcript_strand(strand: char, ts: Option<char>, xs: Option<char>) -> char {
    match (ts, xs) {
        (Some(t @ ('+' | '-')), _) => t,
        (_, Some(x @ ('+' | '-'))) if strand == '+' || strand == '-' => {
            if x == strand { '+' } else { '-' }
        }
        _ => '*',
    }
}

/// Optional tags read from a PAF line. Unmapped records keep the defaults.
struct PafTags<'a> {
    ms: i64,
    aln_score: i64,
    cs: &'a str,
    /// '*' when absent, as for an unknown `strand`.
    tp_tag: char,
    /// Raw `ts:A` / `XS:A` values, combined by [`transcript_strand`].
    ts: Option<char>,
    xs: Option<char>,
}

/// Scan the optional-tag substring for `ms:i:`, `AS:i:`, `cs:Z:`, `tp:A:`,
/// `ts:A:` and `XS:A:`. Unmapped records skip the scan and keep the defaults
/// (`0`, `0`, `""`, `'*'`, none, none).
fn parse_tags(raw: &str, is_unmapped: bool) -> PafTags<'_> {
    let mut t = PafTags { ms: 0, aln_score: 0, cs: "", tp_tag: '*', ts: None, xs: None };
    if is_unmapped {
        return t;
    }

    for field in raw.split('\t') {
        let b = field.as_bytes();
        if b.len() < 6 { continue; }
        match &b[..5] {
            b"ms:i:" => { t.ms        = parse_signed_int(&b[5..]); }
            b"AS:i:" => { t.aln_score = parse_signed_int(&b[5..]); }
            b"cs:Z:" => { t.cs        = &field[5..]; }
            b"tp:A:" => { t.tp_tag    = b[5] as char; }
            b"ts:A:" => { t.ts        = Some(b[5] as char); }
            b"XS:A:" => { t.xs        = Some(b[5] as char); }
            _ => {}
        }
    }

    t
}

#[inline]
fn parse_signed_int(b: &[u8]) -> i64 {
    // SAFETY: SAM/PAF integer tags are ASCII; slice is valid UTF-8.
    unsafe { std::str::from_utf8_unchecked(b) }
        .parse()
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ts_is_read_relative_and_used_as_written() {
        for strand in ['+', '-'] {
            assert_eq!(transcript_strand(strand, Some('+'), None), '+');
            assert_eq!(transcript_strand(strand, Some('-'), None), '-');
        }
    }

    #[test]
    fn xs_is_reference_relative_and_combined_with_the_alignment_strand() {
        assert_eq!(transcript_strand('+', None, Some('+')), '+');
        assert_eq!(transcript_strand('-', None, Some('-')), '+');
        assert_eq!(transcript_strand('+', None, Some('-')), '-');
        assert_eq!(transcript_strand('-', None, Some('+')), '-');
    }

    #[test]
    fn ts_wins_over_xs_and_anything_else_is_unknown() {
        assert_eq!(transcript_strand('+', Some('-'), Some('+')), '-');
        assert_eq!(transcript_strand('+', None, None), '*');
        assert_eq!(transcript_strand('+', Some('?'), None), '*');
        assert_eq!(transcript_strand('+', Some('?'), Some('-')), '-', "an unusable ts falls back to XS");
        assert_eq!(transcript_strand('*', None, Some('+')), '*');
    }

    #[test]
    fn same_read_on_a_reverse_complemented_reference_keeps_its_transcript_strand() {
        // Aligned to FWD, the read maps '+'; to the RC reference it maps '-'.
        // minimap2 writes the same ts; STAR's XS flips with the reference.
        assert_eq!(transcript_strand('+', Some('-'), None), transcript_strand('-', Some('-'), None));
        assert_eq!(transcript_strand('+', None, Some('-')), transcript_strand('-', None, Some('+')));
    }

    #[test]
    fn parse_line_reads_the_tags() {
        let base = "r1\t100\t0\t100\t-\tchr1\t1000\t10\t110\t100\t100\t60\ttp:A:P";
        let ts = |tags: &str| parse_line(&format!("{base}{tags}"), 1).unwrap().transcript_strand;
        assert_eq!(ts("\tts:A:+"), '+');
        assert_eq!(ts("\tXS:A:+"), '-', "reference '+' on a '-' alignment is antisense");
        assert_eq!(ts(""), '*');
        let unmapped = parse_line("r2\t50\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\tts:A:+", 1).unwrap();
        assert_eq!(unmapped.transcript_strand, '*');
    }

    /// One real read (mate 2 of `A01045:715:HMN2KDSX2:3:1104:3224:5885`, HG002
    /// short-read RNA-seq) as aligned by minimap2 2.31 (`-ax splice:sr`) and STAR
    /// (`--outSAMstrandField intronMotif`) to HG38 and to the reverse-complemented
    /// HG38. The intron's motif is `ct..ac` on the forward reference, so the
    /// transcript is on the reference's `-` strand and this read is antisense to it.
    /// All four records must give `transcript_strand` `-`, through sam2paf. QUAL is
    /// trimmed to `*`, and SEQ for the minimap2 records (their cs:Z makes it unused).
    #[test]
    fn real_minimap2_and_star_records_agree_across_reference_orientations() {
        use crate::sam2paf::convert::{convert, Options};
        const HEADER: &str = "@SQ\tSN:chr22\tLN:50818468\n";
        let q = "A01045:715:HMN2KDSX2:3:1104:3224:5885";
        let records = [
            // minimap2, HG38: '+' alignment, ts:A:- (read-relative).
            format!("{q}\t163\tchr22\t17138248\t60\t110M192N40M\t=\t17138068\t0\t*\t*\tNM:i:0\tts:A:-\tcs:Z::110~ct192ac:40"),
            // minimap2, reverse-complemented HG38: '-' alignment, same ts:A:-.
            format!("{q}\t147\tchr22\t33679880\t60\t40M192N110M\t=\t33680060\t0\t*\t*\tNM:i:0\tts:A:-\tcs:Z::40~gt192ag:110"),
            // STAR, HG38: '+' alignment, XS:A:- (reference-relative).
            format!("{q}\t163\tchr22\t17138248\t255\t110M192N40M\t=\t17138068\t0\tGTTGCTGCCGTGTGCCCCCGGCCCCTAGTTCTGGCGCCCCATCATGCGTTGCCTTCTGCAGGTACTGGTGGAACAGGTTGGCGCCGTATACGTCAGACATAGGGTTATCACCCACAGCATAGAGCTTCCGGATGGGGGCGGCCCAGCCCC\t*\tNH:i:1\tXS:A:-\tNM:i:0\tMD:Z:150"),
            // STAR, reverse-complemented HG38: '-' alignment, XS flips to '+'.
            format!("{q}\t147\tchr22\t33679880\t255\t40M192N110M\t=\t33680060\t0\tGGGGCTGGGCCGCCCCCATCCGGAAGCTCTATGCTGTGGGTGATAACCCTATGTCTGACGTATACGGCGCCAACCTGTTCCACCAGTACCTGCAGAAGGCAACGCATGATGGGGCGCCAGAACTAGGGGCCGGGGGCACACGGCAGCAAC\t*\tNH:i:1\tXS:A:+\tNM:i:0\tMD:Z:150"),
        ];
        let expected_strand = ['+', '-', '+', '-'];
        for (rec, strand) in records.iter().zip(expected_strand) {
            let mut out = Vec::new();
            convert(format!("{HEADER}{rec}\n").as_bytes(), &mut out, &Options::default()).unwrap();
            let paf = String::from_utf8(out).unwrap();
            let r = parse_line(paf.trim_end(), 1).unwrap();
            assert_eq!(r.strand, strand, "{paf}");
            assert_eq!(r.transcript_strand, '-', "{paf}");
        }
    }
}
