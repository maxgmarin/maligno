//! Alignment-input dispatch: lets `compare` and `sam2paf` take SAM/BAM as well as PAF.
//!
//! Format is sniffed per file with htslib (`hts_open`), so mixed inputs (e.g. BAM vs
//! PAF) are fine. SAM/BAM never gets its own conversion logic: a BAM is rendered
//! back to SAM text record by record with htslib's `sam_format1` (a built-in
//! `samtools view -h`) and fed to the unchanged `sam2paf` converter. The PAF lines
//! `compare` sees are therefore byte-identical to `samtools view -h | maligno sam2paf`.
//!
//! CRAM is deliberately unsupported (it would need a reference FASTA).

use std::ffi::CString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use anyhow::{bail, Context, Result};
use rust_htslib::bam::{self, Read as _, Record};
use rust_htslib::htslib;

use crate::io_utils::open_input;
use crate::sam2paf::convert::{convert, Options};

/// BGZF decompression threads for a BAM reader (htslib thread pool).
const BAM_THREADS: usize = 2;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlnFormat {
    Paf,
    /// SAM text; `compressed` = gzip/BGZF (read through htslib), else plain text.
    Sam { compressed: bool },
    Bam,
}

impl AlnFormat {
    pub fn is_paf(self) -> bool {
        self == AlnFormat::Paf
    }
    pub fn name(self) -> &'static str {
        match self {
            AlnFormat::Paf => "PAF",
            AlnFormat::Sam { .. } => "SAM",
            AlnFormat::Bam => "BAM",
        }
    }
}

/// Sniff the format of `path` with htslib. Anything htslib doesn't call SAM or BAM
/// (plain or gzipped PAF included) is treated as PAF; CRAM is an error.
pub fn detect(path: &str) -> Result<AlnFormat> {
    let c_path = CString::new(path).with_context(|| format!("invalid path '{path}'"))?;
    // SAFETY: `c_path` is a valid NUL-terminated string; the handle is closed
    // before returning and its `format` struct is copied out first.
    let fmt = unsafe {
        let fp = htslib::hts_open(c_path.as_ptr(), c"r".as_ptr());
        if fp.is_null() {
            bail!("cannot open '{path}'");
        }
        let fmt = (*fp).format;
        htslib::hts_close(fp);
        fmt
    };
    Ok(match fmt.format {
        htslib::htsExactFormat_bam => AlnFormat::Bam,
        htslib::htsExactFormat_sam => AlnFormat::Sam {
            compressed: fmt.compression != htslib::htsCompression_no_compression,
        },
        htslib::htsExactFormat_cram => {
            bail!("'{path}' is CRAM, which maligno does not support; convert it to BAM first")
        }
        _ => AlnFormat::Paf,
    })
}

/// The `SO:` value of the `@HD` header line of a SAM/BAM (e.g. "coordinate"), if any.
pub fn sort_order(path: &str) -> Result<Option<String>> {
    let reader =
        bam::Reader::from_path(path).with_context(|| format!("cannot open '{path}'"))?;
    let text = String::from_utf8_lossy(reader.header().as_bytes()).into_owned();
    Ok(text
        .lines()
        .find(|l| l.starts_with("@HD\t"))
        .and_then(|hd| hd.split('\t').find_map(|f| f.strip_prefix("SO:")))
        .map(str::to_owned))
}

/// A SAM/BAM file as a stream of SAM text: plain SAM is read as-is (exactly what
/// `sam2paf` always did); compressed SAM and BAM go through htslib (`BamAsSam`).
pub fn open_sam_text(path: &str, fmt: AlnFormat) -> Result<Box<dyn BufRead>> {
    match fmt {
        AlnFormat::Sam { compressed: false } => {
            let file = File::open(path).with_context(|| format!("cannot open '{path}'"))?;
            Ok(Box::new(BufReader::with_capacity(1 << 20, file)))
        }
        AlnFormat::Sam { compressed: true } | AlnFormat::Bam => {
            Ok(Box::new(BamAsSam::open(path)?))
        }
        AlnFormat::Paf => bail!("internal error: '{path}' is PAF, not SAM/BAM"),
    }
}

/// Hand `consume` the file at `path` as PAF lines. A PAF is opened directly; a
/// SAM/BAM is converted on a worker thread (`sam2paf` with `opts`) and streamed
/// through a pipe, so no intermediate PAF is written.
///
/// Error precedence: a converter failure (e.g. a record with no cs/MD) wins over
/// whatever it caused downstream (a truncated stream); a converter "broken pipe"
/// only means `consume` stopped early, so `consume`'s own error is reported.
pub fn with_paf_lines<T>(
    path: &str,
    fmt: AlnFormat,
    opts: &Options,
    consume: impl FnOnce(Box<dyn BufRead>) -> Result<T>,
) -> Result<T> {
    if fmt.is_paf() {
        return consume(open_input(path)?);
    }
    let (pipe_r, pipe_w) = io::pipe().context("cannot create conversion pipe")?;
    std::thread::scope(|s| {
        let worker = s.spawn(move || -> Result<()> {
            let mut w = BufWriter::with_capacity(1 << 20, pipe_w);
            convert(open_sam_text(path, fmt)?, &mut w, opts)?;
            w.flush()?;
            Ok(())
        });
        let result = consume(Box::new(BufReader::with_capacity(1 << 20, pipe_r)));
        let converted = worker.join().expect("SAM/BAM conversion thread panicked");
        match converted {
            Err(e) if !is_broken_pipe(&e) => {
                Err(e).with_context(|| format!("converting {} '{path}' to PAF", fmt.name()))
            }
            _ => result,
        }
    })
}

fn is_broken_pipe(e: &anyhow::Error) -> bool {
    e.chain()
        .filter_map(|c| c.downcast_ref::<io::Error>())
        .any(|io| io.kind() == io::ErrorKind::BrokenPipe)
}

// ── BAM → SAM text ────────────────────────────────────────────────────────────

/// Reads a BAM (or compressed SAM) through htslib and yields it as SAM text: the
/// header, then one `sam_format1` line per record — the same bytes as
/// `samtools view -h`.
pub struct BamAsSam {
    reader: bam::Reader,
    record: Record,
    line: htslib::kstring_t,
    buf: Vec<u8>,
    pos: usize,
    done: bool,
}

/// Refill `buf` with at least this many bytes of records per `fill_buf` call.
const BATCH_BYTES: usize = 1 << 20;

impl BamAsSam {
    pub fn open(path: &str) -> Result<Self> {
        let mut reader =
            bam::Reader::from_path(path).with_context(|| format!("cannot open '{path}'"))?;
        reader
            .set_threads(BAM_THREADS)
            .with_context(|| format!("cannot set read threads for '{path}'"))?;
        let mut buf = reader.header().as_bytes().to_vec();
        if !buf.is_empty() && !buf.ends_with(b"\n") {
            buf.push(b'\n');
        }
        Ok(BamAsSam {
            reader,
            record: Record::new(),
            line: htslib::kstring_t { l: 0, m: 0, s: std::ptr::null_mut() },
            buf,
            pos: 0,
            done: false,
        })
    }

    /// Append the next record's SAM line to `buf`; false at end of file.
    fn push_next_record(&mut self) -> io::Result<bool> {
        match self.reader.read(&mut self.record) {
            None => return Ok(false),
            Some(r) => r.map_err(io::Error::other)?,
        }
        self.line.l = 0;
        // SAFETY: header and record are live htslib objects owned by `self`;
        // `line` is a kstring htslib grows as needed (freed in `Drop`).
        let bytes = unsafe {
            let hdr = self.reader.header().inner_ptr() as *const htslib::sam_hdr_t;
            if htslib::sam_format1(hdr, self.record.inner(), &mut self.line) < 0 {
                return Err(io::Error::other("htslib failed to format a record as SAM"));
            }
            std::slice::from_raw_parts(self.line.s as *const u8, self.line.l as usize)
        };
        self.buf.extend_from_slice(bytes);
        self.buf.push(b'\n');
        Ok(true)
    }
}

impl io::Read for BamAsSam {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let n = {
            let avail = self.fill_buf()?;
            let n = avail.len().min(out.len());
            out[..n].copy_from_slice(&avail[..n]);
            n
        };
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for BamAsSam {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.pos == self.buf.len() && !self.done {
            self.buf.clear();
            self.pos = 0;
            while self.buf.len() < BATCH_BYTES {
                if !self.push_next_record()? {
                    self.done = true;
                    break;
                }
            }
        }
        Ok(&self.buf[self.pos..])
    }

    fn consume(&mut self, amt: usize) {
        self.pos = (self.pos + amt).min(self.buf.len());
    }
}

impl Drop for BamAsSam {
    fn drop(&mut self) {
        if !self.line.s.is_null() {
            // SAFETY: `line.s` was allocated by htslib (malloc/realloc).
            unsafe { htslib::free(self.line.s as *mut std::ffi::c_void) };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read as _;

    const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";
    const REC1: &str = "r1\t0\tchr1\t101\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tMD:Z:10";
    const REC2: &str = "r2\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII";

    fn tmp(name: &str) -> String {
        std::env::temp_dir()
            .join(format!("maligno_test_aln_input_{}_{name}", std::process::id()))
            .to_string_lossy()
            .into_owned()
    }

    /// Write `HEADER` + records as a SAM file, and as a BAM via htslib.
    fn write_fixtures(stem: &str) -> (String, String) {
        let sam = tmp(&format!("{stem}.sam"));
        std::fs::write(&sam, format!("{HEADER}{REC1}\n{REC2}\n")).unwrap();
        let bam = tmp(&format!("{stem}.bam"));
        let view = bam::HeaderView::from_bytes(HEADER.as_bytes());
        let mut w =
            bam::Writer::from_path(&bam, &bam::Header::from_template(&view), bam::Format::Bam)
                .unwrap();
        for rec in [REC1, REC2] {
            w.write(&Record::from_sam(&view, rec.as_bytes()).unwrap()).unwrap();
        }
        drop(w);
        (sam, bam)
    }

    #[test]
    fn detects_paf_sam_bam() {
        let (sam, bam) = write_fixtures("detect");
        let paf = tmp("detect.paf");
        std::fs::write(&paf, "r1\t10\t0\t10\t+\tchr1\t1000\t100\t110\t10\t10\t60\ttp:A:P\n")
            .unwrap();
        assert_eq!(detect(&paf).unwrap(), AlnFormat::Paf);
        assert_eq!(detect(&sam).unwrap(), AlnFormat::Sam { compressed: false });
        assert_eq!(detect(&bam).unwrap(), AlnFormat::Bam);
        assert_eq!(sort_order(&bam).unwrap().as_deref(), Some("coordinate"));
        for p in [sam, bam, paf] {
            let _ = std::fs::remove_file(p);
        }
    }

    #[test]
    fn bam_reads_back_as_the_original_sam_text() {
        let (sam, bam) = write_fixtures("roundtrip");
        let mut text = String::new();
        BamAsSam::open(&bam).unwrap().read_to_string(&mut text).unwrap();
        assert!(text.contains("@SQ\tSN:chr1\tLN:1000\n"));
        assert!(text.ends_with(&format!("{REC1}\n{REC2}\n")), "got:\n{text}");
        for p in [sam, bam] {
            let _ = std::fs::remove_file(p);
        }
    }
}
