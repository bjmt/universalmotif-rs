//! Motif operations against universalmotif.

mod common;

use common::*;
use universalmotif::ops::{
    TrimFrom, normalize, reverse_complement, round_motif, trim_cwm, trim_motif,
};

#[test]
fn ops_match_r() {
    for (i, case) in load("ops.json").iter().enumerate() {
        let op = case["op"].as_str().unwrap();
        let input = motif_from_record(&case["input"]);
        let a = &case["args"];
        let ctx = format!("case {i} ({op} on {}, args {a})", input.name);
        let trim_from = || {
            a["trim_from"]
                .as_str()
                .unwrap()
                .parse::<TrimFrom>()
                .unwrap()
        };
        let mut skip = Skip::default();
        let got = match op {
            "motif_rc" => reverse_complement(&input, false).map(Some),
            "round_motif" => {
                skip.derived = true;
                round_motif(&input, num(&a["pct_tolerance"])).map(Some)
            }
            "normalize" => {
                skip.derived = true;
                normalize(&input).map(Some)
            }
            "trim_motifs" => {
                // universalmotif stores nsites = 100 in trimmed motifs that had
                // none; this crate leaves it empty.
                skip.nsites = input.nsites.is_none();
                trim_motif(&input, num(&a["min_ic"]), trim_from())
            }
            "trim_cwm" => trim_cwm(
                &input,
                num(&a["trim_threshold"]),
                opt_num(&a["abs_threshold"]),
                trim_from(),
            ),
            _ => panic!("unknown op {op}"),
        };
        let got = got.unwrap_or_else(|e| panic!("{ctx}: {e}"));
        match (&case["motif"], got) {
            (serde_json::Value::Null, None) => {
                assert!(
                    case["error"]
                        .as_str()
                        .unwrap()
                        .contains("completely trimmed"),
                    "{ctx}"
                );
            }
            (serde_json::Value::Null, Some(m)) => {
                panic!("{ctx}: R failed with {}, but got\n{m}", case["error"])
            }
            (want, Some(m)) => {
                assert_motif_eq_skip(&m, want, &ctx, skip);
                if skip.derived {
                    // The port recomputes the total IC from the new matrix.
                    assert!(close(m.icscore, own_icscore(&m)), "{ctx}: icscore");
                }
                if skip.nsites {
                    assert_eq!(m.nsites, None, "{ctx}: nsites");
                }
            }
            (_, None) => panic!("{ctx}: unexpectedly trimmed away"),
        }
    }
}
