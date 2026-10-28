use finge_rs::{
    CountEcfpFingerprint, EcfpFingerprint, EcfpGraph, Fingerprint, InvalidatedGraph,
    LayeredCountEcfpFingerprint, MorganMassPolicy, PatternFingerprint, RdkFingerprint,
    SmilesRdkitScratch,
};
use smiles_rs::Smiles;

#[test]
fn opt_in_native_mass_metadata_matches_fixed_runtime_regressions() {
    // RDKit 2026.09.1 wheel c700cb37..., source fece8caa...; sanitized parsing,
    // chirality=false, ring membership/bond types=true. Baseline evidence keeps
    // all native bytes and unfolded hashes for radii 0..=3 at 2048 bits.
    for (text, scientific_delta, native_delta, seed, bit) in [
        ("[162Tm]", -6, -7, 3_275_783_702, 1558),
        ("[239Th]", 7, 6, 496_497_325, 685),
        ("[119Tc]", 20, 21, 625_568_618, 874),
        ("[255Db]", -12, -13, 1_386_265_684, 1108),
    ] {
        let molecule: Smiles = text.parse().unwrap();
        let mut scratch = SmilesRdkitScratch::default();
        for retain_h in [false, true] {
            let graph = if retain_h {
                scratch.try_prepare_checked_preserving_explicit_hydrogens(&molecule)
            } else {
                scratch.try_prepare_checked(&molecule)
            }
            .unwrap();
            let native = graph.with_morgan_mass_policy(MorganMassPolicy::Rdkit2026091);
            assert_eq!(
                graph.ecfp_atom_invariant_fields(0).delta_mass,
                scientific_delta
            );
            assert_eq!(
                native.ecfp_atom_invariant_fields(0).delta_mass,
                native_delta
            );
            assert_eq!(native.ecfp_atom_invariant(0, true), seed);
            let mut native_fields = native.ecfp_atom_invariant_fields(0);
            native_fields.delta_mass = scientific_delta;
            assert_eq!(native_fields, graph.ecfp_atom_invariant_fields(0));
            assert_eq!(
                PatternFingerprint::default().compute(&native),
                PatternFingerprint::default().compute(&graph)
            );
            assert_eq!(
                RdkFingerprint::default().compute(&native),
                RdkFingerprint::default().compute(&graph)
            );
            for radius in 0..=3 {
                let fingerprint = EcfpFingerprint::new(radius, 2048);
                assert_eq!(fingerprint.unfolded(&native).as_slice(), &[seed]);
                assert_eq!(
                    fingerprint
                        .compute(&native)
                        .active_bits()
                        .collect::<Vec<_>>(),
                    vec![bit]
                );
                assert_eq!(
                    CountEcfpFingerprint::new(radius, 2048)
                        .compute(&native)
                        .active_counts()
                        .collect::<Vec<_>>(),
                    vec![(bit, 1)]
                );
                let layered = LayeredCountEcfpFingerprint::new(radius, 2048).compute(&native);
                assert_eq!(
                    layered.formula().active_counts().collect::<Vec<_>>(),
                    vec![(bit, 1)]
                );
                assert!(
                    layered.layers()[1..]
                        .iter()
                        .all(|layer| layer.active_counts().next().is_none())
                );
                let restored = native.with_morgan_mass_policy(MorganMassPolicy::ElementsRs);
                assert_eq!(fingerprint.compute(&restored), fingerprint.compute(&graph));
                let wrapped = InvalidatedGraph::new(native);
                assert_eq!(fingerprint.compute(&wrapped), fingerprint.compute(&native));
            }
        }
    }
}

#[test]
fn mass_policy_does_not_extend_parser_acceptance() {
    for text in ["[0C]", "[21C]", "[999C]", "[1000C]", "[65535C]"] {
        assert!(text.parse::<Smiles>().is_err(), "{text}");
    }
}
