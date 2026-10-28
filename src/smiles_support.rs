//! Public `smiles-rs` adapter surface for RDKit-parity fingerprinting.

pub use crate::MorganMassPolicy;

pub use crate::smiles_support_impl::{
    CheckedSmilesPreparationError, SmilesPreparationError, SmilesRdkitGraph, SmilesRdkitScratch,
};

/// Backward-compatible alias for the RDKit-normalized `smiles-rs` graph view.
pub type SmilesEcfpGraph<'a> = SmilesRdkitGraph<'a>;

/// Backward-compatible alias for the reusable RDKit-normalization scratch state.
pub type SmilesEcfpScratch = SmilesRdkitScratch;
