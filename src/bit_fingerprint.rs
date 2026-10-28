use alloc::{vec, vec::Vec};
use bitvec::prelude::{BitSlice, BitVec, Lsb0};

/// Dense bit fingerprint backed by `bitvec`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BitFingerprint {
    bits: BitVec<usize, Lsb0>,
}

impl BitFingerprint {
    /// Creates an all-zero fingerprint with the requested bit length.
    #[inline]
    #[must_use]
    pub fn zeros(len: usize) -> Self {
        Self {
            bits: BitVec::repeat(false, len),
        }
    }

    /// Returns the underlying bit slice.
    #[inline]
    #[must_use]
    pub fn as_bitslice(&self) -> &BitSlice<usize, Lsb0> {
        self.bits.as_bitslice()
    }

    /// Sets a bit to `true` when the index is in range.
    #[inline]
    pub fn set(&mut self, index: usize) {
        if index < self.bits.len() {
            self.bits.set(index, true);
        }
    }

    /// Returns whether the bit at `index` is set.
    #[inline]
    #[must_use]
    pub fn contains(&self, index: usize) -> bool {
        self.bits.get(index).is_some_and(|bit| *bit)
    }

    /// Returns the indices of every set bit.
    #[inline]
    pub fn active_bits(&self) -> impl Iterator<Item = usize> + '_ {
        self.bits.iter_ones()
    }

    /// Returns the number of bits in the fingerprint.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.bits.len()
    }

    /// Returns whether the fingerprint stores no bits.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bits.is_empty()
    }

    /// Packs bits into bytes with bit zero in the least-significant position
    /// of the first byte. Unused high bits in the final byte are zero.
    /// This representation does not depend on host word size or endianness.
    #[must_use]
    pub fn to_le_bytes(&self) -> Vec<u8> {
        let mut bytes = vec![0; self.len().div_ceil(8)];
        for bit in self.active_bits() {
            bytes[bit / 8] |= 1 << (bit % 8);
        }
        bytes
    }
}

impl From<BitVec<usize, Lsb0>> for BitFingerprint {
    #[inline]
    fn from(bits: BitVec<usize, Lsb0>) -> Self {
        Self { bits }
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;

    use bitvec::{bitvec, prelude::BitVec};

    use super::BitFingerprint;

    #[test]
    fn bit_fingerprint_tracks_set_bits() {
        let mut fingerprint = BitFingerprint::zeros(16);
        fingerprint.set(1);
        fingerprint.set(7);
        fingerprint.set(99);

        assert_eq!(
            fingerprint.active_bits().collect::<alloc::vec::Vec<_>>(),
            vec![1, 7]
        );
        assert!(fingerprint.contains(1));
        assert!(!fingerprint.contains(2));
        assert_eq!(fingerprint.len(), 16);
    }

    #[test]
    fn bit_fingerprint_exposes_bits_and_empty_state() {
        let empty = BitFingerprint::zeros(0);
        assert!(empty.is_empty());
        assert_eq!(empty.as_bitslice().len(), 0);

        let raw: BitVec<usize, bitvec::order::Lsb0> =
            bitvec![usize, bitvec::order::Lsb0; 1, 0, 1, 0];
        let fingerprint = BitFingerprint::from(raw.clone());

        assert!(!fingerprint.is_empty());
        assert_eq!(fingerprint.as_bitslice(), raw.as_bitslice());
        assert_eq!(
            fingerprint.active_bits().collect::<alloc::vec::Vec<_>>(),
            vec![0, 2]
        );
    }

    #[test]
    fn packed_bytes_follow_lsb0_order_across_byte_boundaries() {
        let mut fingerprint = BitFingerprint::zeros(18);
        for bit in [3, 7, 8, 17] {
            fingerprint.set(bit);
        }
        assert_eq!(fingerprint.to_le_bytes(), vec![0x88, 0x01, 0x02]);
        assert!(BitFingerprint::zeros(0).to_le_bytes().is_empty());
    }
}
