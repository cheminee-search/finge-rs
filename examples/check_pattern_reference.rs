//! Differential validation of a supplied, versioned molecule reference corpus.
//! Every input contributes to the summary, including unsupported inputs.

use std::{
    env,
    fs::File,
    io::{BufRead, BufReader},
};

use finge_rs::{PatternFingerprint, SmilesRdkitScratch};
use serde::Deserialize;
use smiles_rs::smiles::Smiles;

#[derive(Deserialize)]
struct Reference {
    index: usize,
    smiles: String,
    fp_size: usize,
    expected_hex: Option<String>,
    oracle_error: Option<String>,
    remove_hs: Option<bool>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let path = args.next().ok_or("expected reference JSONL path")?;
    let checked = match args.next().as_deref() {
        None => false,
        Some("--checked") => true,
        _ => return Err("only optional --checked is supported".into()),
    };
    let mut scratch = SmilesRdkitScratch::default();
    let (mut rows, mut matched, mut rejected, mut errors, mut mismatches) = (0, 0, 0, 0, 0);
    for line in BufReader::new(File::open(path)?).lines() {
        rows += 1;
        let row: Reference = serde_json::from_str(&line?)?;
        let Some(expected) = row.expected_hex else {
            if row.oracle_error.is_none() {
                return Err("missing oracle outcome".into());
            }
            rejected += 1;
            continue;
        };
        let molecule = match row.smiles.parse::<Smiles>() {
            Ok(molecule) => molecule,
            Err(error) => {
                eprintln!("index={} stage=parse error={error}", row.index);
                errors += 1;
                continue;
            }
        };
        let graph = if checked && !row.remove_hs.unwrap_or(true) {
            scratch
                .try_prepare_checked_preserving_explicit_hydrogens(&molecule)
                .map_err(|error| error.to_string())
        } else if checked {
            scratch
                .try_prepare_checked(&molecule)
                .map_err(|error| error.to_string())
        } else if !row.remove_hs.unwrap_or(true) {
            scratch
                .try_prepare_preserving_explicit_hydrogens(&molecule)
                .map_err(|error| error.to_string())
        } else {
            scratch
                .try_prepare(&molecule)
                .map_err(|error| error.to_string())
        };
        let graph = match graph {
            Ok(graph) => graph,
            Err(error) => {
                eprintln!("index={} stage=prepare error={error}", row.index);
                errors += 1;
                continue;
            }
        };
        match PatternFingerprint::new(row.fp_size).try_compute(&graph) {
            Ok(fingerprint) => {
                let observed = fingerprint
                    .to_le_bytes()
                    .iter()
                    .map(|byte| format!("{byte:02x}"))
                    .collect::<String>();
                if observed == expected {
                    matched += 1;
                } else {
                    mismatches += 1;
                    eprintln!(
                        "index={} width={} fingerprint_mismatch",
                        row.index, row.fp_size
                    );
                }
            }
            Err(error) => {
                eprintln!("index={} stage=fingerprint error={error}", row.index);
                errors += 1;
            }
        }
    }
    println!(
        "{}",
        serde_json::json!({
            "rows": rows, "matched": matched, "oracle_rejected": rejected,
            "errors": errors, "mismatches": mismatches, "checked_preparation": checked,
        })
    );
    if errors != 0 || mismatches != 0 {
        std::process::exit(1);
    }
    Ok(())
}
