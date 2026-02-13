use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// A single FASTQ record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub id: String,
    pub sequence: String,
    pub quality: Option<String>,
}

impl FastqRecord {
    /// Create a new FASTQ record
    pub fn new(id: String, sequence: String, quality: Option<String>) -> Self {
        Self { id, sequence, quality }
    }
}

/// Fast FASTQ reader with buffering and optional gzip support
pub struct FastqReader<R: BufRead> {
    reader: R,
    is_fasta: bool,
    buffer: String,
}

// Enum to hold either a plain file reader or gzipped reader
pub enum FileReader {
    Plain(BufReader<std::fs::File>),
    Gzipped(BufReader<GzDecoder<BufReader<std::fs::File>>>),
}

impl Read for FileReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            FileReader::Plain(r) => r.read(buf),
            FileReader::Gzipped(r) => r.read(buf),
        }
    }
}

impl BufRead for FileReader {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            FileReader::Plain(r) => r.fill_buf(),
            FileReader::Gzipped(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            FileReader::Plain(r) => r.consume(amt),
            FileReader::Gzipped(r) => r.consume(amt),
        }
    }
}

impl FastqReader<FileReader> {
    /// Open a FASTQ file (auto-detects gzip)
    pub fn from_path(path: impl AsRef<Path>, is_fasta: bool) -> Result<Self> {
        let file = std::fs::File::open(path.as_ref())
            .with_context(|| format!("Failed to open file: {:?}", path.as_ref()))?;

        // Check if file is gzipped by reading magic bytes
        let mut buffered = BufReader::new(file);
        let is_gzipped = {
            let peek = buffered.fill_buf()?;
            peek.len() >= 2 && peek[0] == 0x1f && peek[1] == 0x8b
        };

        let reader = if is_gzipped {
            let decoder = GzDecoder::new(buffered);
            FileReader::Gzipped(BufReader::new(decoder))
        } else {
            FileReader::Plain(buffered)
        };

        Ok(Self::new(reader, is_fasta))
    }
}

impl<R: BufRead> FastqReader<R> {
    /// Create a new FASTQ reader from any BufRead type
    pub fn new(reader: R, is_fasta: bool) -> Self {
        Self {
            reader,
            is_fasta,
            buffer: String::with_capacity(512), // Pre-allocate for typical read lengths
        }
    }

    /// Read the next FASTQ record
    pub fn next(&mut self) -> Result<Option<FastqRecord>> {
        // Read ID line
        self.buffer.clear();
        let bytes_read = self.reader.read_line(&mut self.buffer)?;
        if bytes_read == 0 {
            return Ok(None); // EOF
        }

        let id = self.buffer.trim_end().to_string();

        // Read sequence line
        self.buffer.clear();
        self.reader.read_line(&mut self.buffer)
            .context("Invalid FASTQ: missing sequence line")?;
        let sequence = self.buffer.trim_end().to_string();

        if self.is_fasta {
            return Ok(Some(FastqRecord::new(id, sequence, None)));
        }

        // Read comment line ('+')
        self.buffer.clear();
        self.reader.read_line(&mut self.buffer)
            .context("Invalid FASTQ: missing '+' line")?;

        // Read quality line
        self.buffer.clear();
        self.reader.read_line(&mut self.buffer)
            .context("Invalid FASTQ: missing quality line")?;
        let quality = self.buffer.trim_end().to_string();

        Ok(Some(FastqRecord::new(id, sequence, Some(quality))))
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_fastq_parsing() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let cursor = BufReader::new(Cursor::new(data));
        let mut reader = FastqReader::new(cursor, false);

        let record1 = reader.next().unwrap().unwrap();
        assert_eq!(record1.id, "@read1");
        assert_eq!(record1.sequence, "ACGT");
        assert_eq!(record1.quality, Some("IIII".to_string()));

        let record2 = reader.next().unwrap().unwrap();
        assert_eq!(record2.id, "@read2");
        assert!(reader.next().unwrap().is_none());
    }
}
