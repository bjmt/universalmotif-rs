//! Motif construction defaults against universalmotif's `universalmotif_cpp()`.

mod common;

use common::*;
use universalmotif::{GapInfo, Matrix, Motif, Strand};

#[test]
fn builder_matches_r() {
    for (i, case) in load("build.json").iter().enumerate() {
        let inp = &case["input"];
        let ctx = format!("case {i} (input {inp})");
        let alph = alphabet(inp["alphabet"].as_str().unwrap_or("DNA"));
        let matrix = Matrix::from_columns(&columns(&inp["matrix"])).unwrap();
        let mut b = Motif::builder(matrix, alph)
            .motif_type(inp["type"].as_str().map(motif_type))
            .nsites(opt_num(&inp["nsites"]))
            .icscore(opt_num(&inp["icscore"]))
            .bkg(inp.get("bkg").filter(|v| !v.is_null()).map(nums))
            .altname(opt_str(&inp["altname"]))
            .family(opt_str(&inp["family"]))
            .organism(opt_str(&inp["organism"]))
            .bkgsites(opt_num(&inp["bkgsites"]))
            .pval(opt_num(&inp["pval"]))
            .qval(opt_num(&inp["qval"]))
            .eval(opt_num(&inp["eval"]))
            .extrainfo(extrainfo(&inp["extrainfo"]));
        if let Some(p) = opt_num(&inp["pseudocount"]) {
            b = b.pseudocount(p);
        }
        if let Some(n) = inp["name"].as_str() {
            b = b.name(n);
        }
        if let Some(s) = inp["strand"].as_str() {
            b = b.strand(s.parse::<Strand>().unwrap());
        }
        if inp["isgapped"].as_bool() == Some(true) {
            let u = |k: &str| nums(&inp[k]).iter().map(|&x| x as usize).collect();
            b = b.gap(GapInfo {
                isgapped: true,
                gaploc: u("gaploc"),
                mingap: u("mingap"),
                maxgap: u("maxgap"),
            });
        }
        let got = b.build();
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
