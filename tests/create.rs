//! `create_motif()` against universalmotif.

mod common;

use common::*;
use universalmotif::create::{CreateOptions, from_matrix, from_strings};
use universalmotif::{Matrix, MotifType};

#[test]
fn create_matches_r() {
    for (i, case) in load("create.json").iter().enumerate() {
        let inp = &case["input"];
        let a = &case["args"];
        let ctx = format!("case {i} (input {inp}, args {a})");
        let mut opts = CreateOptions::default();
        if let Some(t) = a["type"].as_str() {
            opts.motif_type = motif_type(t);
        }
        if let Some(n) = a["name"].as_str() {
            opts.name = n.into();
        }
        if let Some(p) = opt_num(&a["pseudocount"]) {
            opts.pseudocount = p;
        }
        opts.nsites = opt_num(&a["nsites"]);
        if !a["bkg"].is_null() {
            opts.bkg = Some(nums(&a["bkg"]));
        }
        let alph = a["alphabet"].as_str().map(alphabet);

        let got = match inp["kind"].as_str().unwrap() {
            "matrix" => {
                let cols = columns(&inp["matrix"]);
                let rows: Vec<char> = inp["rows"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|r| r.as_str().unwrap().chars().next().unwrap())
                    .collect();
                let m = Matrix::from_columns(&cols).unwrap();
                let got = from_matrix(m, Some(&rows), alph, &opts);
                // universalmotif drops the row names of single-column
                // matrices while sorting them and then fails (see
                // DIVERGENCES.md); this crate handles them.
                if cols.len() == 1 && case["motif"].is_null() {
                    let m = got.unwrap_or_else(|e| panic!("{ctx}: {e}"));
                    assert_eq!(m.motif_type, MotifType::Pcm);
                    continue;
                }
                got
            }
            "strings" => {
                let s: Vec<&str> = inp["strings"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|s| s.as_str().unwrap())
                    .collect();
                from_strings(&s, alph, &opts)
            }
            k => panic!("unknown input kind {k}"),
        };
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
