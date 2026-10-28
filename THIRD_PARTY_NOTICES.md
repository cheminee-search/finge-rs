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
