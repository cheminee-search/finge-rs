# Third-party notices

The Pattern fingerprint module (`src/fingerprints/pattern.rs`) derives its
thirteen pattern definitions and fingerprint contract from RDKit's
`Code/GraphMol/Fingerprints/PatternFingerprints.cpp` at commit
`fece8caa860bdf6c9c82bdb9253f22a0cf87150c` (RDKit 2026.09.1).

Copyright (C) 2013–2020 Greg Landrum and Rational Discovery LLC.
RDKit is distributed under the BSD 3-Clause License. The exact license text
and project copyright notice are retained in `LICENSE-RDKIT`; distributions
must retain the required notices. This attribution does not imply endorsement.

The Rust implementation reuses finge-rs's existing `morgan_hash_combine` helper.
It does not introduce a runtime dependency on RDKit, C++, or Boost. Other
finge-rs source remains covered by the MIT license in `LICENSE`; the combined
package declares `MIT AND BSD-3-Clause`.

The optional Morgan mass-invariant metadata in
`src/morgan_mass/rdkit2026091.rs` derives from RDKit's
`Code/GraphMol/atomic_data.cpp` at the same pinned commit
`fece8caa860bdf6c9c82bdb9253f22a0cf87150c` (RDKit 2026.09.1), exact Git blob
`fbd8341e740e01b18bc7e2d72f71d0bcb4b46cba`.
Copyright (C) 2001–2024 Greg Landrum and other RDKit contributors.
Its mass-delta semantics follow that version's `Atom::getMass` and
`MorganFingerprints::getConnectivityInvariants`: native weight/isotope metadata,
integer-mass fallback for missing isotopes, and truncation toward zero.
The generator and runtime-reference fixture retain source/version provenance.
This versioned compatibility metadata does not replace elements-rs scientific
tables or make RDKit a runtime dependency. The BSD 3-Clause license and notice
are retained in `LICENSE-RDKIT`.
