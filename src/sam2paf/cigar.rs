/// Iterator over (length, op_byte) pairs in a SAM CIGAR string.
///
/// Operates directly on the raw bytes of the CIGAR string; performs no
/// heap allocation and does no regex matching.
pub struct CigarIter<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> CigarIter<'a> {
    pub fn new(s: &'a str) -> Self {
        CigarIter { bytes: s.as_bytes(), pos: 0 }
    }
}

impl<'a> Iterator for CigarIter<'a> {
    type Item = (u32, u8);

    #[inline]
    fn next(&mut self) -> Option<(u32, u8)> {
        if self.pos >= self.bytes.len() {
            return None;
        }
        // parse integer length
        let start = self.pos;
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        if self.pos >= self.bytes.len() {
            return None; // malformed: digits at end with no op
        }
        // SAFETY: all bytes in start..self.pos are ASCII digits → valid UTF-8
        let len: u32 = unsafe {
            std::str::from_utf8_unchecked(&self.bytes[start..self.pos])
                .parse()
                .unwrap_unchecked()
        };
        let op = self.bytes[self.pos];
        self.pos += 1;
        Some((len, op))
    }
}

/// Cumulative stats extracted from a single-pass CIGAR walk.
#[derive(Debug, Default)]
pub struct CigarStats {
    /// Sum of M / = / X lengths (all "aligned" reference+query ops).
    pub m: u32,
    /// Number of I operations.
    pub i_count: u32,
    /// Total inserted bases.
    pub i_bases: u32,
    /// Number of D operations.
    pub d_count: u32,
    /// Total deleted bases.
    pub d_bases: u32,
    /// Total N (intron skip) bases.
    pub n_bases: u32,
    /// Total soft-clipped bases (S ops only).
    pub soft_clip: u32,
    /// Leading clip in query coordinates (S or H at the start).
    pub clip_lead: u32,
    /// Trailing clip in query coordinates (S or H at the end).
    pub clip_trail: u32,
    /// Mismatch count from X ops (extended CIGAR only).
    pub mm_ext: u32,
    /// True if any M op was seen (standard CIGAR).
    pub have_m: bool,
    /// True if any = or X op was seen (extended CIGAR).
    pub have_ext: bool,
    /// Total number of CIGAR operations (used for overflow warning).
    pub n_ops: u32,
}

impl CigarStats {
    /// Walk `cigar_str` once, populating all fields.
    pub fn parse(cigar_str: &str) -> Self {
        let mut s = CigarStats::default();
        for (len, op) in CigarIter::new(cigar_str) {
            match op {
                b'M' => { s.m += len; s.have_m = true; }
                b'=' => { s.m += len; s.have_ext = true; }
                b'X' => { s.m += len; s.mm_ext += len; s.have_ext = true; }
                b'I' => { s.i_count += 1; s.i_bases += len; }
                b'D' => { s.d_count += 1; s.d_bases += len; }
                b'N' => { s.n_bases += len; }
                b'S' => {
                    let slot = if s.n_ops == 0 { &mut s.clip_lead } else { &mut s.clip_trail };
                    *slot = len;
                    s.soft_clip += len;
                }
                b'H' => {
                    let slot = if s.n_ops == 0 { &mut s.clip_lead } else { &mut s.clip_trail };
                    *slot = len;
                }
                _ => {}
            }
            s.n_ops += 1;
        }
        s
    }

    /// Total reference-consuming length (M + D + N).
    #[inline] pub fn ref_len(&self) -> u32 { self.m + self.d_bases + self.n_bases }
    /// Total query-consuming length (M + I + S).
    #[inline] pub fn query_consumed(&self) -> u32 { self.m + self.i_bases + self.soft_clip }
    /// Full query length including hard clips (M + I + leading + trailing).
    #[inline] pub fn query_len(&self) -> u32 { self.m + self.i_bases + self.clip_lead + self.clip_trail }
    /// Alignment block length (M + I + D, excluding introns).
    #[inline] pub fn block_len(&self) -> u32 { self.m + self.i_bases + self.d_bases }
}

/// A compacted representation of CIGAR operations used during cs-tag generation.
///
/// Hard-clip ops are dropped, and = / X are normalised to M.  Consecutive ops
/// of the same kind are merged so the MD-walk needs fewer iterations.
pub fn build_merged_cigar(cigar_str: &str) -> Vec<(u32, u8)> {
    let mut out: Vec<(u32, u8)> = Vec::new();
    for (len, op) in CigarIter::new(cigar_str) {
        if op == b'H' { continue; }
        let norm = if op == b'=' || op == b'X' { b'M' } else { op };
        if let Some(last) = out.last_mut() {
            if last.1 == norm { last.0 += len; continue; }
        }
        out.push((len, norm));
    }
    out
}

/// Write the CIGAR string with S and H ops removed directly to `w`.
/// Avoids an intermediate String allocation on the hot path.
#[inline]
pub fn write_cigar_no_clips<W: std::io::Write>(w: &mut W, cigar_str: &str) -> std::io::Result<()> {
    // Fast path: if no S or H bytes are present, write the raw string.
    if !cigar_str.bytes().any(|b| b == b'S' || b == b'H') {
        return w.write_all(cigar_str.as_bytes());
    }
    for (len, op) in CigarIter::new(cigar_str) {
        if op == b'S' || op == b'H' { continue; }
        write!(w, "{}{}", len, op as char)?;
    }
    Ok(())
}
