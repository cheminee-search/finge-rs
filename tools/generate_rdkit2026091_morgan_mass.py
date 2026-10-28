#!/usr/bin/env python3
"""Generate versioned invariant metadata from the exact pinned RDKit source.

No RDKit installation or scientific mass dependency is used by this generator.
"""
import argparse
import hashlib
import pathlib
import re

EXPECTED_BLOB = "fbd8341e740e01b18bc7e2d72f71d0bcb4b46cba"
COMMIT = "fece8caa860bdf6c9c82bdb9253f22a0cf87150c"


def read_native(path):
    data = path.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    assert blob == EXPECTED_BLOB, "Source must be the exact pinned native Git blob"
    source = data.decode()
    statement = source.split("const std::string periodicTableAtomData =", 1)[1].split(')DAT";', 1)[0] + ')DAT"'
    weights = {}
    for line in "".join(re.findall(r'R"DAT\((.*?)\)DAT"', statement, re.S)).splitlines():
        if not line.strip():
            continue
        fields = line.split()
        assert len(fields) >= 11
        weights.setdefault(int(fields[0]), float(fields[6]))
    assert sorted(weights) == list(range(119))
    statement = source.split("const std::string isotopesAtomData[] =", 1)[1].split('"EOS"', 1)[0]
    isotopes = {}
    space_only = []
    for payload in re.findall(r'R"DAT\((.*?)\)DAT"', statement, re.S):
        for line in payload.splitlines():
            if not line.strip():
                continue
            fields = line.split()
            assert len(fields) >= 5
            key = int(fields[0]), int(fields[2])
            assert key not in isotopes
            isotopes[key] = float(fields[3])
            if "\t" not in line:
                space_only.append(key)
    # The final 54 space-separated rows are part of the compiled table.
    assert len(isotopes) == 3111
    assert len(space_only) == 54 and (111, 272) in space_only
    assert {z for z, _ in space_only} == set(range(111, 119))
    return weights, isotopes, hashlib.sha256(data).hexdigest()


def render_array(name, kind, values, width):
    lines = ["#[rustfmt::skip]", f"pub(super) static {name}: [{kind}; {len(values)}] = ["]
    for offset in range(0, len(values), width):
        lines.append("    " + ", ".join(values[offset:offset + width]) + ",")
    lines.append("];")
    return "\n".join(lines)


def generate(path, full):
    weights, masses, source_sha = read_native(path)
    deltas = {key: int(mass - weights[key[0]]) for key, mass in masses.items()}
    overrides = {key: value for key, value in deltas.items()
                 if value != int(key[1] - weights[key[0]])}
    assert len(overrides) == 626
    # Confirm compact representation reconstructs every native-listed delta.
    for (number, isotope), expected in deltas.items():
        assert overrides.get((number, isotope), int(isotope - weights[number])) == expected
    entries = sorted((deltas if full else overrides).items())
    assert all(0 < z <= 118 and 0 < isotope <= 65535 for (z, isotope), _ in entries)
    assert all(-32768 <= value <= 32767 for _, value in entries)
    keys = [(number << 16) | isotope for (number, isotope), _ in entries]
    assert keys == sorted(set(keys))
    representation = "all listed isotope deltas" if full else "delta overrides relative to integer-mass fallback"
    header = f'''// SPDX-License-Identifier: BSD-3-Clause
// Generated RDKit 2026.09.1 Morgan mass-invariant metadata, not scientific masses.
// Derived from Code/GraphMol/atomic_data.cpp at {COMMIT}.
// Native Git blob: {EXPECTED_BLOB}
// Native source SHA256: {source_sha}
// Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors.
// See LICENSE-RDKIT and THIRD_PARTY_NOTICES.md. Regenerate with the pinned-source
// generator; do not edit entries. Representation: {representation}.
// The source contains 119 weights and 3111 listed isotopes, including all 54
// space-separated tail rows. Compact representation contains 626 overrides.
'''
    return header + "\n" + "\n\n".join([
        render_array("STANDARD_WEIGHTS", "f64", [repr(weights[z]) for z in range(119)], 8),
        render_array("ISOTOPE_KEYS", "u32", [f"0x{key:08x}" for key in keys], 8),
        render_array("MASS_DELTAS", "i16", [str(value) for _, value in entries], 16),
    ]) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--full", action="store_true", help="External equal-output comparison only")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = generate(args.source, args.full)
    if args.check:
        assert args.output.read_text() == rendered, "Generated source differs"
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(hashlib.sha256(rendered.encode()).hexdigest())


if __name__ == "__main__":
    main()
