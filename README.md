# finge-rs

[![CI](https://github.com/earth-metabolome-initiative/finge-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/earth-metabolome-initiative/finge-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Rust 1.92+](https://img.shields.io/badge/rust-1.92%2B-orange.svg)](Cargo.toml)

finge-rs computes molecular fingerprints in `no_std` Rust (with `extern crate alloc`), built around a small set of traits for molecular atoms, bonds, and graphs so any backend can supply the chemistry. Where it claims RDKit parity it reproduces RDKit's bit output exactly, checked against a tracked corpus of RDKit-parseable SMILES from the scikit-fingerprints HIV set across ECFP radii `0` to `5` and AtomPair, RDK, and Topological Torsion sizes from `64` to `4096`.

It provides the usual connectivity fingerprints. Morgan/ECFP comes as a bit fingerprint (`EcfpFingerprint`), a folded-count fingerprint (`CountEcfpFingerprint`), and an exact-radius folded-count variant (`LayeredCountEcfpFingerprint`), and alongside them are 2D AtomPair (`AtomPairFingerprint`), RDKit topological paths (`RdkFingerprint`), Topological Torsion (`TopologicalTorsionFingerprint`), and the 167-bit MACCS keys (`MaccsFingerprint`) behind the `smarts-support` feature. Each folds into a `BitFingerprint` or `CountFingerprint`.

MAP4 (`Map4Fingerprint`, MinHashed atom pair) rounds out the set. It exposes its raw shingle set, folds into a bit or counted fingerprint like the others, and produces a MinHash signature for similarity work, validated for Tanimoto parity against the reference `map4` implementation. ECFP can be sketched to a MinHash the same way. Those signatures feed `LshIndex`, a banded locality-sensitive-hashing index for approximate nearest-neighbour search over large collections.

`SmilesRdkitScratch` turns a `smiles-rs` molecule into an RDKit-normalized graph that every fingerprint accepts. AtomPair and Topological Torsion also run on a raw `Smiles` when you do not need the normalization step.

`PatternFingerprint` implements RDKit 2026.09.1's thirteen molecule patterns
with tautomeric mode disabled. It is separate from `RdkFingerprint`, which
implements topological paths. `PatternFingerprint::default().try_compute(&graph)`
returns a complete 2048-bit result or a typed error; `with_max_matches` can bound
the non-unique embeddings per pattern without returning silently truncated bits.
`BitFingerprint::to_le_bytes()` provides portable least-significant-bit-first
packing, including 256 bytes for a 2048-bit result.

The Pattern adapter accepts prepared `SmilesRdkitGraph` values and honors their
hydrogen policy. It supports single, double, triple, quadruple and aromatic
bonds; other bond families return an error. The concrete `smiles-rs` parser
currently rejects dummy/wildcard atoms even though custom concrete graph
backends can supply atomic number zero. SMARTS queries, MolBundle intersection,
tautomeric mode and substructure-screen guarantees are outside this API.
Preparation and graph sanitization remain caller responsibilities. The fixed
pattern definitions derive from BSD-licensed RDKit; see `THIRD_PARTY_NOTICES.md`
and `LICENSE-RDKIT`.

The existing preparation methods preserve input aromatic labels when standalone kekulization fails. Use `try_prepare_checked` or `try_prepare_checked_preserving_explicit_hydrogens` to reject that fallback with `CheckedSmilesPreparationError::Kekulization`. These additive methods leave existing preparation behavior and `SmilesPreparationError` unchanged. They check localization and propagate aromaticity-assignment errors, but do not perform full molecular sanitization, certify aromaticity-perception completeness, or guarantee acceptance by RDKit. Some RDKit-accepted structures, including `N=c1[nH]c(=S)ss1` and `O=c1[nH]c(=O)ss1` from the tracked reference corpus, cannot currently be localized and are rejected by the checked path.

Fingerprinting, `LshIndex` and `TanimotoIndex` support `wasm32-unknown-unknown`, including the `smarts-support` feature.

## Usage

The example below prepares a few molecules, folds ECFP and MAP4 fingerprints, and builds a MAP4 MinHash index to retrieve a molecule's nearest neighbours.

```rust
use finge_rs::{EcfpFingerprint, Fingerprint, LshIndex, Map4Fingerprint};
use finge_rs::smiles_support::SmilesRdkitScratch;
use smiles_rs::smiles::Smiles;

let mut scratch = SmilesRdkitScratch::default();
let mut index = LshIndex::<u32, 512, 256>::new();

for smiles in ["CCO", "OCC1OC(O)C(O)C(O)C1O", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"] {
    let molecule: Smiles = smiles.parse().expect("example SMILES should parse");
    let graph = scratch.try_prepare(&molecule).expect("preparation should succeed");

    let _ecfp = EcfpFingerprint::default().compute(&graph);
    let _map4_bits = Map4Fingerprint::default().compute(&graph);
    index.insert(Map4Fingerprint::default().minhash::<_, u32, 512>(&graph));
}

let query: Smiles = "CCO".parse().expect("example SMILES should parse");
let graph = scratch.try_prepare(&query).expect("preparation should succeed");
let signature = Map4Fingerprint::default().minhash::<_, u32, 512>(&graph);
let neighbours = index.query(&signature, 3);
assert!(!neighbours.is_empty());
```
