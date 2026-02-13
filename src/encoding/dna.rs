/// DNA base encoding (2 bits per base)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Base {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

impl Base {
    pub fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'A' => Some(Base::A),
            'C' => Some(Base::C),
            'G' => Some(Base::G),
            'T' => Some(Base::T),
            _ => None,
        }
    }

    pub fn to_char(self) -> char {
        match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
        }
    }

    #[allow(dead_code)]
    pub fn complement(self) -> Self {
        match self {
            Base::A => Base::T,
            Base::T => Base::A,
            Base::C => Base::G,
            Base::G => Base::C,
        }
    }
}

/// Encode DNA sequence to bit-packed representation
#[allow(dead_code)]
pub fn encode_sequence(seq: &str) -> Vec<u8> {
    let mut result = vec![0u8; (seq.len() + 3) / 4]; // 4 bases per byte
    
    for (i, c) in seq.chars().enumerate() {
        if let Some(base) = Base::from_char(c) {
            let byte_idx = i / 4;
            let bit_offset = (i % 4) * 2;
            result[byte_idx] |= (base as u8) << bit_offset;
        }
    }
    
    result
}

/// Decode bit-packed DNA sequence
#[allow(dead_code)]
pub fn decode_sequence(data: &[u8], length: usize) -> String {
    let mut result = String::with_capacity(length);
    
    for i in 0..length {
        let byte_idx = i / 4;
        let bit_offset = (i % 4) * 2;
        let bits = (data[byte_idx] >> bit_offset) & 0b11;
        
        let base = match bits {
            0b00 => Base::A,
            0b01 => Base::C,
            0b10 => Base::G,
            0b11 => Base::T,
            _ => unreachable!(),
        };
        
        result.push(base.to_char());
    }
    
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode() {
        let seq = "ACGTACGT";
        let encoded = encode_sequence(seq);
        let decoded = decode_sequence(&encoded, seq.len());
        assert_eq!(seq, decoded);
    }

    #[test]
    fn test_complement() {
        assert_eq!(Base::A.complement(), Base::T);
        assert_eq!(Base::T.complement(), Base::A);
        assert_eq!(Base::C.complement(), Base::G);
        assert_eq!(Base::G.complement(), Base::C);
    }
}
