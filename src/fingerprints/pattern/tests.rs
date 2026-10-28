use alloc::{string::String, vec, vec::Vec};

use geometric_traits::traits::{Graph, MonopartiteGraph, MonoplexGraph};
use serde::Deserialize;
use smiles_rs::{
    bond::{Bond, bond_edge::BondEdge},
    smiles::Smiles,
};

use super::{PatternFingerprint, PatternFingerprintError, Target};
use crate::{
    MolecularGraph, PatternBondProperties, PatternFingerprintGraph, SmilesRdkitGraph,
    SmilesRdkitScratch,
};

fn fingerprint(input: &str) -> crate::BitFingerprint {
    let molecule: Smiles = input.parse().unwrap();
    let mut scratch = SmilesRdkitScratch::default();
    PatternFingerprint::default()
        .try_compute(&scratch.try_prepare(&molecule).unwrap())
        .unwrap()
}

#[test]
fn curated_rdkit_2026_09_1_bytes_match() {
    #[derive(Deserialize)]
    struct Case {
        smiles: String,
        fp_size: usize,
        expected_hex: Option<String>,
        oracle_error: Option<String>,
        remove_hs: Option<bool>,
    }
    let mut matched = 0;
    let mut oracle_rejected = 0;
    let mut parser_unsupported = 0;
    let mut scratch = SmilesRdkitScratch::default();
    for line in include_str!("../../../tests/fixtures/rdkit_pattern_curated.jsonl")
        .lines()
        .chain(
            include_str!("../../../tests/fixtures/rdkit_pattern_hydrogens_removed.jsonl").lines(),
        )
        .chain(
            include_str!("../../../tests/fixtures/rdkit_pattern_hydrogens_retained.jsonl").lines(),
        )
    {
        let case: Case = serde_json::from_str(line).unwrap();
        let Some(expected) = case.expected_hex else {
            assert!(case.oracle_error.is_some());
            oracle_rejected += 1;
            continue;
        };
        let molecule = match case.smiles.parse::<Smiles>() {
            Ok(molecule) => molecule,
            Err(error) => {
                // The current concrete Smiles adapter does not admit dummy
                // atoms. Retain both oracle rows as explicit domain gaps.
                assert!(matches!(case.smiles.as_str(), "*" | "*C"));
                assert!(alloc::format!("{error}").contains("Wildcard atom not allowed"));
                parser_unsupported += 1;
                continue;
            }
        };
        let graph = if case.remove_hs.unwrap_or(true) {
            scratch.try_prepare(&molecule)
        } else {
            scratch.try_prepare_preserving_explicit_hydrogens(&molecule)
        }
        .unwrap();
        let observed = PatternFingerprint::new(case.fp_size)
            .try_compute(&graph)
            .unwrap();
        let observed_hex = observed
            .to_le_bytes()
            .iter()
            .map(|byte| alloc::format!("{byte:02x}"))
            .collect::<String>();
        assert_eq!(
            observed_hex, expected,
            "input={} size={}",
            case.smiles, case.fp_size
        );
        matched += 1;
    }
    // All declared inputs remain in the fixture; RDKit-rejected inputs are
    // recorded outcomes, not counted as molecule-fingerprint parity successes.
    assert_eq!((matched, parser_unsupported, oracle_rejected), (146, 6, 9));
}

#[test]
fn count_bits_include_nonunique_embeddings_and_limits_fail_explicitly() {
    let molecule: Smiles = "CC".parse().unwrap();
    let mut scratch = SmilesRdkitScratch::default();
    let graph = scratch.try_prepare_checked(&molecule).unwrap();
    assert_eq!(
        PatternFingerprint::default()
            .with_max_matches(1)
            .try_compute(&graph),
        Err(PatternFingerprintError::MatchLimitExceeded {
            pattern: 1,
            limit: 1
        })
    );
    assert_eq!(
        PatternFingerprint::default()
            .with_max_matches(2)
            .try_compute(&graph)
            .unwrap(),
        fingerprint("CC")
    );
    assert_eq!(
        PatternFingerprint::new(0).try_compute(&graph),
        Err(PatternFingerprintError::InvalidSize(0))
    );
    assert_eq!(
        PatternFingerprint::default()
            .with_max_matches(0)
            .try_compute(&graph),
        Err(PatternFingerprintError::InvalidMatchLimit(0))
    );
    assert_eq!(
        PatternFingerprint::default()
            .with_max_matches(100_000_001)
            .try_compute(&graph),
        Err(PatternFingerprintError::InvalidMatchLimit(100_000_001))
    );
}

#[test]
fn equivalent_atom_orders_preserve_all_pattern_bits() {
    for (left, right) in [
        ("CC(O)C", "OC(C)C"),
        ("CC.CN", "NC.CC"),
        ("c1ccncc1", "n1ccccc1"),
    ] {
        assert_eq!(fingerprint(left), fingerprint(right));
    }
    assert_eq!(fingerprint("C"), fingerprint("[13CH4]"));
    assert_eq!(fingerprint("N"), fingerprint("[NH4+]"));
}

#[test]
fn prepared_hydrogen_policy_is_respected() {
    let molecule: Smiles = "[H]CO".parse().unwrap();
    let mut scratch = SmilesRdkitScratch::default();
    let collapsed = PatternFingerprint::default()
        .try_compute(&scratch.try_prepare(&molecule).unwrap())
        .unwrap();
    assert_eq!(collapsed, fingerprint("CO"));
    let retained = PatternFingerprint::default()
        .try_compute(
            &scratch
                .try_prepare_preserving_explicit_hydrogens(&molecule)
                .unwrap(),
        )
        .unwrap();
    assert_ne!(retained, collapsed);
}

#[test]
fn bridges_do_not_become_ring_bonds_between_ring_atoms() {
    let molecule: Smiles = "C1CC1-C2CC2".parse().unwrap();
    let mut scratch = SmilesRdkitScratch::default();
    let target = Target::new(&scratch.try_prepare(&molecule).unwrap()).unwrap();
    assert!(target.in_ring.iter().all(|&value| value));
    assert_eq!(
        target
            .adjacency
            .iter()
            .flatten()
            .filter(|edge| edge.in_ring)
            .count(),
        12
    );
    assert_eq!(
        target
            .adjacency
            .iter()
            .flatten()
            .filter(|edge| !edge.in_ring)
            .count(),
        2
    );
}

#[derive(Clone, Copy)]
enum Fault {
    MissingAtom,
    Outside,
    SelfLoop,
    Asymmetric,
    Duplicate,
    Inconsistent,
    AtomicNumber,
    BondType,
}

struct Malformed<'a> {
    graph: SmilesRdkitGraph<'a>,
    fault: Fault,
}

impl Graph for Malformed<'_> {
    fn has_nodes(&self) -> bool {
        self.graph.has_nodes()
    }
    fn has_edges(&self) -> bool {
        self.graph.has_edges()
    }
}
impl MonoplexGraph for Malformed<'_> {
    type Edge = <Smiles as MonoplexGraph>::Edge;
    type Edges = <Smiles as MonoplexGraph>::Edges;
    fn edges(&self) -> &Self::Edges {
        self.graph.edges()
    }
}
impl MonopartiteGraph for Malformed<'_> {
    type NodeId = usize;
    type NodeSymbol = <Smiles as MonopartiteGraph>::NodeSymbol;
    type Nodes = <Smiles as MonopartiteGraph>::Nodes;
    fn nodes_vocabulary(&self) -> &Self::Nodes {
        self.graph.nodes_vocabulary()
    }
}
impl MolecularGraph for Malformed<'_> {
    type Bond = BondEdge;
    fn atom(&self, atom: usize) -> Option<&Self::NodeSymbol> {
        if matches!(self.fault, Fault::MissingAtom) && atom == 1 {
            None
        } else {
            self.graph.atom(atom)
        }
    }
    fn bonds(&self, atom: usize) -> impl Iterator<Item = Self::Bond> + '_ {
        let mut bonds = self.graph.bonds(atom).collect::<Vec<_>>();
        match self.fault {
            Fault::Outside => bonds = vec![BondEdge::new(atom, 2, Bond::Single, None)],
            Fault::SelfLoop => bonds = vec![BondEdge::new(atom, atom, Bond::Single, None)],
            Fault::Asymmetric if atom == 1 => bonds.clear(),
            Fault::Duplicate if atom == 0 => bonds.extend(bonds.clone()),
            Fault::Inconsistent if atom == 1 => {
                bonds = vec![BondEdge::new(0, 1, Bond::Double, None)]
            }
            _ => {}
        }
        bonds.into_iter()
    }
}
impl PatternFingerprintGraph for Malformed<'_> {
    fn pattern_atomic_number(&self, atom: usize) -> u32 {
        if matches!(self.fault, Fault::AtomicNumber) {
            200
        } else {
            self.graph.pattern_atomic_number(atom)
        }
    }
    fn pattern_bond_properties(&self, bond: &Self::Bond) -> PatternBondProperties {
        let mut properties = self.graph.pattern_bond_properties(bond);
        if matches!(self.fault, Fault::BondType) {
            properties.bond_type = 17;
        }
        properties
    }
}

#[test]
fn malformed_or_unsupported_backend_data_never_returns_a_partial_fingerprint() {
    let molecule: Smiles = "CC".parse().unwrap();
    let mut scratch = SmilesRdkitScratch::default();
    let graph = scratch.try_prepare(&molecule).unwrap();
    for fault in [
        Fault::MissingAtom,
        Fault::Outside,
        Fault::SelfLoop,
        Fault::Asymmetric,
        Fault::Duplicate,
        Fault::Inconsistent,
        Fault::AtomicNumber,
        Fault::BondType,
    ] {
        assert!(
            PatternFingerprint::default()
                .try_compute(&Malformed { graph, fault })
                .is_err()
        );
    }
}
