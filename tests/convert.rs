//! `convert_type()` against universalmotif for every pair of types.

mod common;

use common::*;
use umotif::convert::convert_type;

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
            (want, Ok(m)) => assert_motif_eq(&m, want, &ctx),
            (_, Err(e)) => panic!("{ctx}: {e}"),
        }
    }
}
