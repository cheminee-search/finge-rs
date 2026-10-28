// SPDX-License-Identifier: MIT AND BSD-3-Clause
// Pattern definitions and fingerprint contract derived from RDKit:
// Copyright (C) 2013-2020 Greg Landrum and Rational Discovery LLC.
// See LICENSE-RDKIT and THIRD_PARTY_NOTICES.md.

use alloc::{vec, vec::Vec};

use crate::{
    BitFingerprint, Fingerprint, MolecularAtom, MolecularBond, PatternBondProperties,
    PatternFingerprintGraph, morgan_hash_combine,
};

const RDKIT_MAX_MATCHES: usize = 100_000_000;

#[cfg(test)]
mod tests;

/// RDKit-style Pattern fingerprint for concrete molecular graphs.
///
/// Implements the thirteen patterns of RDKit 2026.09.1, with tautomeric mode
/// disabled. This is distinct from the RDK path fingerprint. SMARTS query
/// graphs, query bundles and tautomeric fingerprints are outside this API.
/// Graph preparation is the caller's responsibility; this does not sanitize.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PatternFingerprint {
    fp_size: usize,
    max_matches: usize,
}

/// Failure to produce a complete Pattern fingerprint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PatternFingerprintError {
    /// Fingerprint width must be nonzero and fit RDKit's unsigned 32-bit size.
    InvalidSize(usize),
    /// The per-pattern limit must be between one and 100,000,000.
    InvalidMatchLimit(usize),
    /// The graph did not provide an atom at a dense node identifier.
    MissingAtom(usize),
    /// Atomic numbers outside zero (dummy) through 118 are unsupported.
    UnsupportedAtomicNumber { atom: usize, atomic_number: u32 },
    /// A bond has invalid endpoints, is duplicated, or has inconsistent views.
    InvalidBond { source: usize, target: usize },
    /// Only single, double, triple, quadruple and aromatic bonds are supported.
    UnsupportedBondType {
        source: usize,
        target: usize,
        bond_type: u32,
    },
    /// Enumeration exceeded the configured limit; no partial result is returned.
    MatchLimitExceeded { pattern: usize, limit: usize },
}

impl core::fmt::Display for PatternFingerprintError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidSize(size) => write!(f, "invalid Pattern fingerprint size: {size}"),
            Self::InvalidMatchLimit(limit) => write!(f, "invalid Pattern match limit: {limit}"),
            Self::MissingAtom(atom) => write!(f, "missing Pattern graph atom {atom}"),
            Self::UnsupportedAtomicNumber {
                atom,
                atomic_number,
            } => write!(
                f,
                "unsupported atomic number {atomic_number} at atom {atom}"
            ),
            Self::InvalidBond { source, target } => {
                write!(f, "invalid Pattern graph bond {source}-{target}")
            }
            Self::UnsupportedBondType {
                source,
                target,
                bond_type,
            } => write!(
                f,
                "unsupported bond type {bond_type} at bond {source}-{target}"
            ),
            Self::MatchLimitExceeded { pattern, limit } => {
                write!(f, "Pattern {pattern} exceeded its match limit of {limit}")
            }
        }
    }
}

impl core::error::Error for PatternFingerprintError {}

impl PatternFingerprint {
    /// Creates a non-tautomeric fingerprint with the given width.
    /// Invalid parameters are reported by [`try_compute`](Self::try_compute).
    #[must_use]
    pub const fn new(fp_size: usize) -> Self {
        Self {
            fp_size,
            max_matches: RDKIT_MAX_MATCHES,
        }
    }

    /// Limits embeddings per pattern, including automorphic mappings.
    /// Exceeding the limit returns an error rather than a partial fingerprint.
    #[must_use]
    pub const fn with_max_matches(mut self, max_matches: usize) -> Self {
        self.max_matches = max_matches;
        self
    }

    /// Returns the fingerprint width.
    #[must_use]
    pub const fn fp_size(self) -> usize {
        self.fp_size
    }

    /// Returns the per-pattern embedding limit.
    #[must_use]
    pub const fn max_matches(self) -> usize {
        self.max_matches
    }

    /// Computes a complete molecule fingerprint from all thirteen patterns.
    ///
    /// The backend must supply a concrete molecular graph under its intended
    /// aromaticity and hydrogen-removal policy. Query predicates must not be
    /// represented as ordinary atoms or bonds.
    ///
    /// # Errors
    /// Returns an error for invalid parameters, malformed or unsupported graph
    /// properties, or an exceeded embedding limit. No partial result is returned.
    pub fn try_compute<G>(&self, graph: &G) -> Result<BitFingerprint, PatternFingerprintError>
    where
        G: PatternFingerprintGraph,
        G::NodeSymbol: MolecularAtom,
    {
        if self.fp_size == 0 || self.fp_size > u32::MAX as usize {
            return Err(PatternFingerprintError::InvalidSize(self.fp_size));
        }
        if self.max_matches == 0 || self.max_matches > RDKIT_MAX_MATCHES {
            return Err(PatternFingerprintError::InvalidMatchLimit(self.max_matches));
        }
        let target = Target::new(graph)?;
        let mut fingerprint = BitFingerprint::zeros(self.fp_size);
        for (index, pattern) in PATTERNS.iter().enumerate() {
            Matcher {
                target: &target,
                pattern,
                pattern_id: index + 1,
                count_seed: (index + 1 + pattern.atoms + pattern.bonds.len()) as u32,
                matches: 0,
                max_matches: self.max_matches,
                mapping: [0; 8],
                used: vec![false; target.atoms.len()],
                fingerprint: &mut fingerprint,
            }
            .visit(0)?;
        }
        Ok(fingerprint)
    }
}

impl Default for PatternFingerprint {
    fn default() -> Self {
        Self::new(2048)
    }
}

impl<G> Fingerprint<G> for PatternFingerprint
where
    G: PatternFingerprintGraph,
    G::NodeSymbol: MolecularAtom,
{
    type Output = Result<BitFingerprint, PatternFingerprintError>;
    fn compute(&self, graph: &G) -> Self::Output {
        self.try_compute(graph)
    }
}

#[derive(Debug, Clone, Copy)]
struct Edge {
    target: usize,
    properties: PatternBondProperties,
    in_ring: bool,
}

struct Target {
    atoms: Vec<Option<u32>>,
    adjacency: Vec<Vec<Edge>>,
    in_ring: Vec<bool>,
}

impl Target {
    fn new<G>(graph: &G) -> Result<Self, PatternFingerprintError>
    where
        G: PatternFingerprintGraph,
        G::NodeSymbol: MolecularAtom,
    {
        let count = graph.atom_count();
        let mut atoms = Vec::with_capacity(count);
        for atom in 0..count {
            if graph.atom(atom).is_none() {
                return Err(PatternFingerprintError::MissingAtom(atom));
            }
            let number = graph.pattern_atomic_number(atom);
            if number > 118 {
                return Err(PatternFingerprintError::UnsupportedAtomicNumber {
                    atom,
                    atomic_number: number,
                });
            }
            atoms.push((!graph.pattern_atom_is_ignored(atom)).then_some(number));
        }
        let mut adjacency = vec![Vec::<Edge>::new(); count];
        for (atom, neighbors) in adjacency.iter_mut().enumerate() {
            for bond in graph.bonds(atom) {
                let source = bond.source();
                let target = bond.target();
                if source >= count
                    || target >= count
                    || source == target
                    || (source != atom && target != atom)
                {
                    return Err(PatternFingerprintError::InvalidBond { source, target });
                }
                if atoms[source].is_none() || atoms[target].is_none() {
                    continue;
                }
                let other = if source == atom { target } else { source };
                let properties = graph.pattern_bond_properties(&bond);
                if !matches!(properties.bond_type, 1..=4 | 12) {
                    return Err(PatternFingerprintError::UnsupportedBondType {
                        source,
                        target,
                        bond_type: properties.bond_type,
                    });
                }
                neighbors.push(Edge {
                    target: other,
                    properties,
                    in_ring: true,
                });
            }
            neighbors.sort_unstable_by_key(|edge| edge.target);
            if let Some(duplicate) = neighbors
                .windows(2)
                .find(|pair| pair[0].target == pair[1].target)
            {
                return Err(PatternFingerprintError::InvalidBond {
                    source: atom,
                    target: duplicate[0].target,
                });
            }
        }
        for (source, edges) in adjacency.iter().enumerate() {
            for edge in edges {
                let reverse = adjacency[edge.target]
                    .binary_search_by_key(&source, |edge| edge.target)
                    .ok()
                    .map(|index| adjacency[edge.target][index]);
                if reverse.is_none_or(|reverse| reverse.properties != edge.properties) {
                    return Err(PatternFingerprintError::InvalidBond {
                        source,
                        target: edge.target,
                    });
                }
            }
        }
        mark_bridges(&mut adjacency);
        let in_ring = adjacency
            .iter()
            .map(|edges| edges.iter().any(|edge| edge.in_ring))
            .collect();
        Ok(Self {
            atoms,
            adjacency,
            in_ring,
        })
    }

    fn edge(&self, source: usize, target: usize) -> Option<Edge> {
        self.adjacency[source]
            .binary_search_by_key(&target, |edge| edge.target)
            .ok()
            .map(|index| self.adjacency[source][index])
    }
}

// Iterative low-link traversal: an undirected edge belongs to a cycle exactly
// when it is not a bridge. The traversal stack does not grow the call stack.
fn mark_bridges(adjacency: &mut [Vec<Edge>]) {
    let count = adjacency.len();
    let mut discovered = vec![usize::MAX; count];
    let mut low = vec![0; count];
    let mut parent = vec![None; count];
    let mut stack = Vec::<(usize, usize)>::new();
    let mut time = 0;
    for root in 0..count {
        if discovered[root] != usize::MAX {
            continue;
        }
        discovered[root] = time;
        low[root] = time;
        time += 1;
        stack.push((root, 0));
        while let Some((atom, next)) = stack.last_mut() {
            let atom = *atom;
            if *next < adjacency[atom].len() {
                let other = adjacency[atom][*next].target;
                *next += 1;
                if discovered[other] == usize::MAX {
                    parent[other] = Some(atom);
                    discovered[other] = time;
                    low[other] = time;
                    time += 1;
                    stack.push((other, 0));
                } else if parent[atom] != Some(other) {
                    low[atom] = low[atom].min(discovered[other]);
                }
            } else {
                stack.pop();
                if let Some(parent) = parent[atom] {
                    low[parent] = low[parent].min(low[atom]);
                    if low[atom] > discovered[parent] {
                        for edge in &mut adjacency[parent] {
                            if edge.target == atom {
                                edge.in_ring = false;
                            }
                        }
                        for edge in &mut adjacency[atom] {
                            if edge.target == parent {
                                edge.in_ring = false;
                            }
                        }
                    }
                }
            }
        }
    }
}

#[derive(Clone, Copy)]
enum BondConstraint {
    Any,
    Ring,
    SingleOrAromatic,
}

struct PatternBond {
    source: usize,
    target: usize,
    constraint: BondConstraint,
}
struct Pattern {
    atoms: usize,
    ring_atoms: u8,
    bonds: &'static [PatternBond],
}

struct Matcher<'a> {
    target: &'a Target,
    pattern: &'a Pattern,
    pattern_id: usize,
    count_seed: u32,
    matches: usize,
    max_matches: usize,
    mapping: [usize; 8],
    used: Vec<bool>,
    fingerprint: &'a mut BitFingerprint,
}

impl Matcher<'_> {
    fn visit(&mut self, depth: usize) -> Result<(), PatternFingerprintError> {
        if depth == self.pattern.atoms {
            return self.emit();
        }
        if let Some(parent) = self
            .pattern
            .bonds
            .iter()
            .find(|bond| bond.target == depth)
            .map(|bond| bond.source)
        {
            let anchor = self.mapping[parent];
            for index in 0..self.target.adjacency[anchor].len() {
                self.try_atom(depth, self.target.adjacency[anchor][index].target)?;
            }
        } else {
            for atom in 0..self.target.atoms.len() {
                self.try_atom(depth, atom)?;
            }
        }
        Ok(())
    }

    fn try_atom(&mut self, depth: usize, atom: usize) -> Result<(), PatternFingerprintError> {
        if self.used[atom]
            || self.target.atoms[atom].is_none()
            || (self.pattern.ring_atoms & (1 << depth) != 0 && !self.target.in_ring[atom])
        {
            return Ok(());
        }
        for bond in self
            .pattern
            .bonds
            .iter()
            .filter(|bond| bond.target == depth)
        {
            let Some(edge) = self.target.edge(self.mapping[bond.source], atom) else {
                return Ok(());
            };
            let matches = match bond.constraint {
                BondConstraint::Any => true,
                BondConstraint::Ring => edge.in_ring,
                BondConstraint::SingleOrAromatic => matches!(edge.properties.bond_type, 1 | 12),
            };
            if !matches {
                return Ok(());
            }
        }
        self.mapping[depth] = atom;
        self.used[atom] = true;
        let result = self.visit(depth + 1);
        self.used[atom] = false;
        result
    }

    fn emit(&mut self) -> Result<(), PatternFingerprintError> {
        if self.matches == self.max_matches {
            return Err(PatternFingerprintError::MatchLimitExceeded {
                pattern: self.pattern_id,
                limit: self.max_matches,
            });
        }
        self.matches += 1;
        morgan_hash_combine(&mut self.count_seed, 0xbeef);
        self.fingerprint
            .set(self.count_seed as usize % self.fingerprint.len());
        let mut typed = self.pattern_id as u32;
        for &atom in &self.mapping[..self.pattern.atoms] {
            morgan_hash_combine(
                &mut typed,
                self.target.atoms[atom].expect("matched atom is present"),
            );
        }
        for bond in self.pattern.bonds {
            let edge = self
                .target
                .edge(self.mapping[bond.source], self.mapping[bond.target])
                .expect("matched bond is present");
            morgan_hash_combine(
                &mut typed,
                if edge.properties.aromatic {
                    12
                } else {
                    edge.properties.bond_type
                },
            );
        }
        self.fingerprint
            .set(typed as usize % self.fingerprint.len());
        Ok(())
    }
}

// Bond order is RDKit SMARTS construction order, including ring closures.
const fn bond(source: usize, target: usize, constraint: BondConstraint) -> PatternBond {
    PatternBond {
        source,
        target,
        constraint,
    }
}

const PATTERNS: &[Pattern] = &[
    // 1: [*]~[*]
    Pattern {
        atoms: 2,
        ring_atoms: 0,
        bonds: &[bond(0, 1, BondConstraint::Any)],
    },
    // 2: [*]~[*]~[*]
    Pattern {
        atoms: 3,
        ring_atoms: 0,
        bonds: &[
            bond(0, 1, BondConstraint::Any),
            bond(1, 2, BondConstraint::Any),
        ],
    },
    // 3: [R]~1~[R]~[R]~1
    Pattern {
        atoms: 3,
        ring_atoms: 7,
        bonds: &[
            bond(0, 1, BondConstraint::Any),
            bond(1, 2, BondConstraint::Any),
            bond(0, 2, BondConstraint::Any),
        ],
    },
    // 4: [*]~[*](~[*])~[*]
    Pattern {
        atoms: 4,
        ring_atoms: 0,
        bonds: &[
            bond(0, 1, BondConstraint::Any),
            bond(1, 2, BondConstraint::Any),
            bond(1, 3, BondConstraint::Any),
        ],
    },
    // 5: [R]~1[R]~[R]~[R]~1
    Pattern {
        atoms: 4,
        ring_atoms: 15,
        bonds: &[
            bond(0, 1, BondConstraint::SingleOrAromatic),
            bond(1, 2, BondConstraint::Any),
            bond(2, 3, BondConstraint::Any),
            bond(0, 3, BondConstraint::Any),
        ],
    },
    // 6: [*]~[*]~[*](~[*])~[*]
    Pattern {
        atoms: 5,
        ring_atoms: 0,
        bonds: &[
            bond(0, 1, BondConstraint::Any),
            bond(1, 2, BondConstraint::Any),
            bond(2, 3, BondConstraint::Any),
            bond(2, 4, BondConstraint::Any),
        ],
    },
    // 7: [R]~1~[R]~[R]~[R]~[R]~1
    Pattern {
        atoms: 5,
        ring_atoms: 31,
        bonds: &[
            bond(0, 1, BondConstraint::Any),
            bond(1, 2, BondConstraint::Any),
            bond(2, 3, BondConstraint::Any),
            bond(3, 4, BondConstraint::Any),
            bond(0, 4, BondConstraint::Any),
        ],
    },
    // 8: [R]~1~[R]~[R]~[R]~[R]~[R]~1
    Pattern {
        atoms: 6,
        ring_atoms: 63,
        bonds: &[
            bond(0, 1, BondConstraint::Any),
            bond(1, 2, BondConstraint::Any),
            bond(2, 3, BondConstraint::Any),
            bond(3, 4, BondConstraint::Any),
            bond(4, 5, BondConstraint::Any),
            bond(0, 5, BondConstraint::Any),
        ],
    },
    // 9: [R](@[R])(@[R])~[R]~[R](@[R])(@[R])
    Pattern {
        atoms: 7,
        ring_atoms: 127,
        bonds: &[
            bond(0, 1, BondConstraint::Ring),
            bond(0, 2, BondConstraint::Ring),
            bond(0, 3, BondConstraint::Any),
            bond(3, 4, BondConstraint::Any),
            bond(4, 5, BondConstraint::Ring),
            bond(4, 6, BondConstraint::Ring),
        ],
    },
    // 10: [R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])
    Pattern {
        atoms: 8,
        ring_atoms: 255,
        bonds: &[
            bond(0, 1, BondConstraint::Ring),
            bond(0, 2, BondConstraint::Ring),
            bond(0, 3, BondConstraint::Any),
            bond(3, 4, BondConstraint::Ring),
            bond(4, 5, BondConstraint::Any),
            bond(5, 6, BondConstraint::Ring),
            bond(5, 7, BondConstraint::Ring),
        ],
    },
    // 11: [*]~[R](@[R])@[R](@[R])~[*]
    Pattern {
        atoms: 6,
        ring_atoms: 30,
        bonds: &[
            bond(0, 1, BondConstraint::Any),
            bond(1, 2, BondConstraint::Ring),
            bond(1, 3, BondConstraint::Ring),
            bond(3, 4, BondConstraint::Ring),
            bond(3, 5, BondConstraint::Any),
        ],
    },
    // 12: [*]~[R](@[R])@[R]@[R](@[R])~[*]
    Pattern {
        atoms: 7,
        ring_atoms: 62,
        bonds: &[
            bond(0, 1, BondConstraint::Any),
            bond(1, 2, BondConstraint::Ring),
            bond(1, 3, BondConstraint::Ring),
            bond(3, 4, BondConstraint::Ring),
            bond(4, 5, BondConstraint::Ring),
            bond(4, 6, BondConstraint::Any),
        ],
    },
    // 13: [*]
    Pattern {
        atoms: 1,
        ring_atoms: 0,
        bonds: &[],
    },
];
