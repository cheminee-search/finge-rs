//! Explicit mass metadata conventions for Morgan atom invariants.

mod rdkit2026091;

/// Mass metadata used to compute Morgan's integer mass-delta invariant.
///
/// This affects only the mass-delta field supplied by a prepared
/// [`SmilesRdkitGraph`](crate::SmilesRdkitGraph). It does not change scientific
/// mass tables, parsing, preparation, or other fingerprints. Selecting a native
/// convention does not establish complete RDKit fingerprint or sanitization
/// equivalence.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MorganMassPolicy {
    /// Existing elements-rs isotope masses and finge-rs standard-weight rules.
    #[default]
    ElementsRs,
    /// RDKit 2026.09.1 mass-invariant metadata, including its isotope fallback.
    ///
    /// This uses that version's standard weights and isotope data, truncating
    /// the mass difference toward zero. Missing isotope entries use the integer
    /// mass number. It does not extend the SMILES parser's accepted isotopes.
    Rdkit2026091,
}

/// Native metadata over the atom representation's stored isotope domain.
pub(crate) fn rdkit_delta_mass(atomic_number: u8, isotope: u16) -> i32 {
    if atomic_number == 0 || isotope == 0 {
        return 0;
    }

    let key = (u32::from(atomic_number) << 16) | u32::from(isotope);
    match rdkit2026091::ISOTOPE_KEYS.binary_search(&key) {
        Ok(index) => i32::from(rdkit2026091::MASS_DELTAS[index]),
        Err(_) => {
            (f64::from(isotope) - rdkit2026091::STANDARD_WEIGHTS[usize::from(atomic_number)]) as i32
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{MorganMassPolicy, rdkit_delta_mass};

    #[test]
    fn native_metadata_preserves_zero_dummy_fallback_and_truncation() {
        assert_eq!(MorganMassPolicy::default(), MorganMassPolicy::ElementsRs);
        for isotope in [0, 1, 13, u16::MAX] {
            assert_eq!(rdkit_delta_mass(0, isotope), 0);
        }
        for number in 1..=118 {
            assert_eq!(rdkit_delta_mass(number, 0), 0);
        }
        assert_eq!(rdkit_delta_mass(69, 162), -7);
        assert_eq!(rdkit_delta_mass(90, 239), 6);
        assert_eq!(rdkit_delta_mass(43, 119), 21);
        // This native-listed tail row is space-separated in atomic_data.cpp.
        assert_eq!(rdkit_delta_mass(111, 272), -8);
        assert_eq!(rdkit_delta_mass(6, u16::MAX), 65_522);
    }

    #[test]
    fn every_listed_isotope_matches_pinned_native_runtime_delta() {
        let mut count = 0;
        for line in include_str!("../../tests/fixtures/rdkit2026091_morgan_mass_deltas.csv").lines()
        {
            if line.starts_with('#') {
                continue;
            }
            let fields = line.split(',').collect::<std::vec::Vec<_>>();
            assert_eq!(fields.len(), 3);
            let number: u8 = fields[0].parse().unwrap();
            let isotope: u16 = fields[1].parse().unwrap();
            let expected: i32 = fields[2].parse().unwrap();
            assert_eq!(
                rdkit_delta_mass(number, isotope),
                expected,
                "{number}/{isotope}"
            );
            count += 1;
        }
        assert_eq!(count, 3111);
    }
}
