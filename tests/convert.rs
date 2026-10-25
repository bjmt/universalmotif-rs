//! `convert_type()` against universalmotif for every pair of types.

mod common;

use common::*;
use universalmotif::MotifType;
use universalmotif::convert::convert_type;

#[test]
fn convert_type_matches_r() {
    let cases = load("convert.json");
    for (i, case) in cases.iter().enumerate() {
        let input = motif_from_record(&case["input"]);
        let a = &case["args"];
        let ctx = format!("case {i} ({} -> args {a})", input.name);
        let got = convert_type(
            &input,
            motif_type(a["type"].as_str().unwrap()),
            opt_num(&a["pseudocount"]),
            a["relative_entropy"].as_bool().unwrap_or(false),
        );
        match (&case["motif"], got) {
            (serde_json::Value::Null, Err(_)) => {}
            (serde_json::Value::Null, Ok(m)) => {
                panic!("{ctx}: R failed with {}, but got\n{m}", case["error"])
            }
            (want, Ok(m)) => {
                // A new pseudocount changes the total IC, which the R package
                // leaves as it was.
                let new_pseudo = opt_num(&a["pseudocount"])
                    .filter(|&p| p != input.pseudocount && m.motif_type != input.motif_type);
                // A CWM has no total IC, so one is computed for the new type.
                let from_cwm = input.motif_type == MotifType::Cwm;
                let skip = Skip {
                    derived: new_pseudo.is_some() || from_cwm,
                    ..Default::default()
                };
                assert_motif_eq_skip(&m, want, &ctx, skip);
                if from_cwm && m.motif_type != MotifType::Cwm {
                    assert!(close(m.icscore, own_icscore(&m)), "{ctx}: icscore");
                }
                if let Some(p) = new_pseudo.filter(|_| !from_cwm) {
                    let mut parts = input.clone().into_parts();
                    parts.pseudocount = p;
                    let want_ic = own_icscore(&universalmotif::Motif::from_parts(parts).unwrap());
                    // The consensus follows the new pseudocount, which
                    // random_reference.rs checks against the R package.
                    assert!(close(m.icscore, want_ic), "{ctx}: icscore");
                }
            }
            (_, Err(e)) => panic!("{ctx}: {e}"),
        }
    }
}
