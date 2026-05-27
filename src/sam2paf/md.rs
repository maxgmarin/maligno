/// A single token from a SAM MD string.
#[derive(Debug)]
pub enum MdToken<'a> {
    /// `(\d+)` — consecutive matching bases.
    Match(u32),
    /// `(\^[A-Za-z]+)` — deleted reference bases (the `^` is stripped).
    Deletion(&'a str),
    /// `([A-Za-z])` — one mismatched reference base.
    Mismatch(u8),
}

/// Zero-copy iterator over the tokens in a SAM MD string.
///
/// The MD grammar (from SAM spec §1.4) is:
///   MD  := [0-9]+ ([A-Za-z]|\^[A-Za-z]+  [0-9]+)*
pub struct MdIter<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> MdIter<'a> {
    pub fn new(md: &'a str) -> Self {
        MdIter { bytes: md.as_bytes(), pos: 0 }
    }
}

impl<'a> Iterator for MdIter<'a> {
    type Item = MdToken<'a>;

    #[inline]
    fn next(&mut self) -> Option<MdToken<'a>> {
        // Skip any trailing semicolons or other stray separators that some
        // callers append (defensive).
        while self.pos < self.bytes.len()
            && !self.bytes[self.pos].is_ascii_alphanumeric()
            && self.bytes[self.pos] != b'^'
        {
            self.pos += 1;
        }

        if self.pos >= self.bytes.len() {
            return None;
        }

        let b = self.bytes[self.pos];

        if b.is_ascii_digit() {
            // Match-length token
            let start = self.pos;
            while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
                self.pos += 1;
            }
            // SAFETY: bytes in range are ASCII digits → valid UTF-8
            let n: u32 = unsafe {
                std::str::from_utf8_unchecked(&self.bytes[start..self.pos])
                    .parse()
                    .unwrap_unchecked()
            };
            Some(MdToken::Match(n))
        } else if b == b'^' {
            // Deletion token
            self.pos += 1; // consume '^'
            let start = self.pos;
            while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_alphabetic() {
                self.pos += 1;
            }
            // SAFETY: bytes in range are ASCII alpha → valid UTF-8
            let bases = unsafe { std::str::from_utf8_unchecked(&self.bytes[start..self.pos]) };
            Some(MdToken::Deletion(bases))
        } else if b.is_ascii_alphabetic() {
            // Mismatch token (single reference base)
            let base = b.to_ascii_lowercase();
            self.pos += 1;
            Some(MdToken::Mismatch(base))
        } else {
            self.pos += 1;
            self.next() // skip unexpected byte
        }
    }
}
