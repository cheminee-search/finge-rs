"""Generate versioned molecule-only PatternFingerprint reference rows.

Development tooling only; RDKit is not a production dependency. Every input
line produces a row for every requested width, including oracle rejections.
"""

import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, DataStructs, rdBase

SOURCE_COMMIT = "fece8caa860bdf6c9c82bdb9253f22a0cf87150c"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--width", type=int, action="append")
    parser.add_argument("--keep-hs", action="store_true")
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        raise RuntimeError(f"expected RDKit 2026.09.1, got {rdBase.rdkitVersion}")
    widths = args.width or [2048]
    if any(width <= 0 or width > 2**32 - 1 for width in widths):
        raise ValueError("fingerprint widths must fit a nonzero uint32")
    content = args.input.read_bytes()
    params = Chem.SmilesParserParams()
    params.allowCXSMILES = False
    params.parseName = False
    params.removeHs = not args.keep_hs
    rows = 0
    rejected = 0
    with args.output.open("w") as output:
        for index, smiles in enumerate(content.decode().splitlines()):
            molecule = Chem.MolFromSmiles(smiles, params)
            for width in widths:
                row = {"index": index, "smiles": smiles, "fp_size": width, "remove_hs": params.removeHs}
                if molecule is None:
                    row.update(expected_hex=None, oracle_error="RDKit sanitized parse rejected")
                    rejected += 1
                else:
                    fingerprint = Chem.PatternFingerprint(
                        molecule, fpSize=width, tautomerFingerprints=False
                    )
                    row.update(
                        expected_hex=DataStructs.BitVectToBinaryText(fingerprint).hex(),
                        oracle_error=None,
                    )
                output.write(json.dumps(row, separators=(",", ":")) + "\n")
                rows += 1
    metadata = {
        "rdkit_version": rdBase.rdkitVersion,
        "rdkit_source_commit": SOURCE_COMMIT,
        "input_sha256": hashlib.sha256(content).hexdigest(),
        "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "rows": rows,
        "oracle_rejections": rejected,
        "widths": widths,
        "tautomeric": False,
        "remove_hs": params.removeHs,
        "allow_cxsmiles": False,
        "parse_name": False,
    }
    args.output.with_suffix(args.output.suffix + ".manifest.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata))


if __name__ == "__main__":
    main()
