# Pattern reference provenance

These molecule fingerprints were generated with RDKit 2026.09.1, corresponding
to RDKit source commit `fece8caa860bdf6c9c82bdb9253f22a0cf87150c`, using
`tests/generate_pattern_reference.py`. The JSONL sidecar manifests record input
and output SHA-256, width, hydrogen policy and parser options.

`pattern_curated.smi` contains 51 declared development inputs. The reference
retains all 153 outcomes at widths 2048, 1000 and 64: 138 byte-exact supported
configurations, six explicitly unsupported concrete wildcard-parser cases
(`*`, `*C` at three widths), and nine RDKit-rejected outcomes (three inputs at
three widths). Unsupported/rejected rows are not counted as fingerprint parity
successes. The checked-preparation API's narrower localization domain is
separate from the fingerprint algorithm; this fixture uses legacy preparation.

`pattern_hydrogens.smi` contains four additional declared development cases:
mapped, isotopic, charged and isolated hydrogens. Both removed and retained
hydrogen policies have four reference rows at width 2048. Their eight supported
byte comparisons bring the total to 146.

To regenerate, install RDKit 2026.09.1 in an isolated development environment:

```sh
python tests/generate_pattern_reference.py --input tests/fixtures/pattern_curated.smi --output tests/fixtures/rdkit_pattern_curated.jsonl --width 2048 --width 1000 --width 64
python tests/generate_pattern_reference.py --input tests/fixtures/pattern_hydrogens.smi --output tests/fixtures/rdkit_pattern_hydrogens_removed.jsonl
python tests/generate_pattern_reference.py --input tests/fixtures/pattern_hydrogens.smi --output tests/fixtures/rdkit_pattern_hydrogens_retained.jsonl --keep-hs
```

The generator emits an explicit oracle-rejection row instead of omitting an
input. `cargo run --example check_pattern_reference -- REFERENCE.jsonl`
checks an external development corpus and returns a failing status for any
parser/preparation/fingerprint error or byte mismatch. Its summary includes
oracle rejections separately. Optional `--checked` selects checked preparation.
The CI unit test asserts the exact known unsupported inventory above while
retaining every original oracle expectation.

Broader public-corpus development/heldout results are separate evidence;
these curated examples are not an exhaustive chemistry or filter-recall claim.
See `THIRD_PARTY_NOTICES.md` for algorithm provenance and licensing.
