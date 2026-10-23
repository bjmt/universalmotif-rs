//! Helpers for reading the reference values exported from universalmotif by
//! `scripts/make_reference.R`.

#![allow(dead_code)]

use serde_json::Value;
use universalmotif::{Alphabet, GapInfo, Matrix, Motif, MotifParts, MotifType, Strand};

pub fn load(file: &str) -> Vec<Value> {
    load_key(file, "cases")
}

pub fn load_key(file: &str, key: &str) -> Vec<Value> {
    let path = format!("{}/tests/reference/{file}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let v: Value = serde_json::from_str(&text).unwrap();
    v[key].as_array().unwrap().clone()
}

/// The text of one of the input files in tests/reference/files.
pub fn reference_file(name: &str) -> String {
    let path = format!(
        "{}/tests/reference/files/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// Function arguments recorded by R; an empty R list() is written as [].
pub fn args(v: &Value) -> serde_json::Map<String, Value> {
    v.as_object().cloned().unwrap_or_default()
}

/// A number, allowing for the strings jsonlite writes for non-finite values.
pub fn num(v: &Value) -> f64 {
    match v {
        Value::Number(n) => n.as_f64().unwrap(),
        Value::String(s) => match s.as_str() {
            "NaN" | "NA" => f64::NAN,
            "Inf" => f64::INFINITY,
            "-Inf" => f64::NEG_INFINITY,
            _ => panic!("not a number: {s}"),
        },
        Value::Bool(b) => *b as u8 as f64,
        _ => panic!("not a number: {v}"),
    }
}

pub fn opt_num(v: &Value) -> Option<f64> {
    if v.is_null() { None } else { Some(num(v)) }
}

pub fn nums(v: &Value) -> Vec<f64> {
    match v {
        Value::Array(a) => a.iter().map(num).collect(),
        Value::Null => Vec::new(),
        _ => vec![num(v)],
    }
}

pub fn opt_str(v: &Value) -> Option<String> {
    v.as_str().map(String::from)
}

pub fn columns(v: &Value) -> Vec<Vec<f64>> {
    v.as_array().unwrap().iter().map(nums).collect()
}

pub fn close(a: f64, b: f64) -> bool {
    if a.is_nan() || b.is_nan() {
        return a.is_nan() && b.is_nan();
    }
    if a.is_infinite() || b.is_infinite() {
        return a == b;
    }
    (a - b).abs() <= 1e-12 * f64::max(1.0, f64::max(a.abs(), b.abs()))
}

pub fn assert_close_vec(got: &[f64], want: &[f64], ctx: &str) {
    assert_eq!(got.len(), want.len(), "{ctx}: length differs");
    for (i, (&g, &w)) in got.iter().zip(want).enumerate() {
        assert!(close(g, w), "{ctx}: element {i}: got {g:?}, want {w:?}");
    }
}

pub fn alphabet(name: &str) -> Alphabet {
    Alphabet::parse(name).unwrap()
}

pub fn motif_type(name: &str) -> MotifType {
    MotifType::parse(name).unwrap()
}

/// Rebuild a motif exactly as recorded, without recomputing anything.
pub fn motif_from_record(r: &Value) -> Motif {
    let cols = columns(&r["matrix"]);
    let usizes = |v: &Value| nums(v).iter().map(|&x| x as usize).collect::<Vec<_>>();
    let parts = MotifParts {
        name: r["name"].as_str().unwrap().into(),
        altname: opt_str(&r["altname"]),
        family: opt_str(&r["family"]),
        organism: opt_str(&r["organism"]),
        matrix: Matrix::from_columns(&cols).unwrap(),
        alphabet: alphabet(r["alphabet"].as_str().unwrap()),
        motif_type: motif_type(r["type"].as_str().unwrap()),
        // The port gives CWMs no total IC.
        icscore: if r["type"] == "CWM" {
            f64::NAN
        } else {
            num(&r["icscore"])
        },
        nsites: opt_num(&r["nsites"]),
        pseudocount: num(&r["pseudocount"]),
        bkg: nums(&r["bkg"]),
        bkgsites: opt_num(&r["bkgsites"]),
        consensus: opt_str(&r["consensus"]),
        strand: Strand::parse(r["strand"].as_str().unwrap()).unwrap(),
        pval: opt_num(&r["pval"]),
        qval: opt_num(&r["qval"]),
        eval: opt_num(&r["eval"]),
        extrainfo: extrainfo(&r["extrainfo"]),
        gap: GapInfo {
            isgapped: r["gap"]["isgapped"].as_bool().unwrap(),
            gaploc: usizes(&r["gap"]["gaploc"]),
            mingap: usizes(&r["gap"]["mingap"]),
            maxgap: usizes(&r["gap"]["maxgap"]),
        },
    };
    Motif::from_parts(parts).unwrap()
}

pub fn extrainfo(v: &Value) -> Vec<(String, String)> {
    match v {
        Value::Array(a) => a
            .iter()
            .map(|p| {
                let p = p.as_array().unwrap();
                (p[0].as_str().unwrap().into(), p[1].as_str().unwrap().into())
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// Numbers read from text can differ from R in the last bit, because R's
/// parser is not always correctly rounded.
pub fn assert_opt_close(got: Option<f64>, want: Option<f64>, ctx: &str) {
    match (got, want) {
        (Some(g), Some(w)) => assert!(close(g, w), "{ctx}: got {g}, want {w}"),
        _ => assert_eq!(got, want, "{ctx}"),
    }
}

/// Fields that are allowed to differ from R for a given comparison.
#[derive(Default, Clone, Copy)]
pub struct Skip {
    pub nsites: bool,
    /// The R package leaves the total IC and consensus as they were after
    /// some operations where the port recomputes them.
    pub derived: bool,
    /// The R package rounds some columns to counts differently, so the
    /// matrix is checked by the caller.
    pub matrix: bool,
}

/// The total IC computed afresh from a motif's own matrix and settings.
pub fn own_icscore(m: &Motif) -> f64 {
    let nsites = m.nsites.unwrap_or(0.0);
    m.matrix
        .columns()
        .map(|c| {
            universalmotif::convert::position_icscore(
                c,
                &m.bkg,
                m.motif_type,
                m.pseudocount,
                nsites,
                false,
            )
        })
        .sum()
}

/// Whether R's rounding of a column of probabilities to `nsites` counts
/// overshoots by more than one count, in which case the port rounds
/// differently.
pub fn r_rounding_differs(pos: &[f64], nsites: f64) -> bool {
    let nsites = if nsites.is_nan() || nsites <= 1.0 {
        100.0
    } else {
        nsites
    };
    let s: f64 = pos.iter().map(|&p| (p * nsites).round()).sum();
    s - nsites > 1.0
}

/// Whether the total IC of a PPM motif goes through a column that R rounds
/// differently.
fn ic_rounds_differently(m: &Motif) -> bool {
    m.motif_type == MotifType::Ppm
        && m.matrix
            .columns()
            .any(|c| r_rounding_differs(c, m.nsites.unwrap_or(0.0)))
}

/// Compare a motif against a record from R, field by field.
pub fn assert_motif_eq(got: &Motif, want: &Value, ctx: &str) {
    assert_motif_eq_skip(got, want, ctx, Skip::default());
}

pub fn assert_motif_eq_skip(got: &Motif, want: &Value, ctx: &str, skip: Skip) {
    let want_motif = motif_from_record(want);
    let w = &*want_motif;
    assert_eq!(got.name, w.name, "{ctx}: name");
    assert_eq!(got.altname, w.altname, "{ctx}: altname");
    assert_eq!(got.family, w.family, "{ctx}: family");
    assert_eq!(got.organism, w.organism, "{ctx}: organism");
    assert_eq!(got.alphabet, w.alphabet, "{ctx}: alphabet");
    assert_eq!(got.motif_type, w.motif_type, "{ctx}: type");
    if got.motif_type == MotifType::Cwm {
        // The R package gives CWMs a total IC; the port has none.
        assert!(got.icscore.is_nan(), "{ctx}: CWM icscore {}", got.icscore);
    } else if skip.derived || skip.matrix || ic_rounds_differently(&want_motif) {
        // Checked by the caller where it matters.
    } else {
        assert!(
            close(got.icscore, w.icscore),
            "{ctx}: icscore: got {}, want {}",
            got.icscore,
            w.icscore
        );
    }
    if !skip.nsites {
        assert_opt_close(got.nsites, w.nsites, &format!("{ctx}: nsites"));
    }
    assert!(close(got.pseudocount, w.pseudocount), "{ctx}: pseudocount");
    assert_close_vec(&got.bkg, &w.bkg, &format!("{ctx}: bkg"));
    assert_opt_close(got.bkgsites, w.bkgsites, &format!("{ctx}: bkgsites"));
    if !skip.derived && !skip.matrix {
        assert_eq!(got.consensus, w.consensus, "{ctx}: consensus");
    }
    assert_eq!(got.strand, w.strand, "{ctx}: strand");
    assert_opt_close(got.pval, w.pval, &format!("{ctx}: pval"));
    assert_opt_close(got.qval, w.qval, &format!("{ctx}: qval"));
    assert_opt_close(got.eval, w.eval, &format!("{ctx}: eval"));
    assert_eq!(got.extrainfo, w.extrainfo, "{ctx}: extrainfo");
    assert_eq!(got.gap, w.gap, "{ctx}: gap");
    assert_eq!(got.matrix.nrow(), w.matrix.nrow(), "{ctx}: nrow");
    assert_eq!(got.matrix.ncol(), w.matrix.ncol(), "{ctx}: ncol");
    if !skip.matrix {
        assert_close_vec(
            got.matrix.as_slice(),
            w.matrix.as_slice(),
            &format!("{ctx}: matrix"),
        );
    }
}

/// The PWM the R package makes from a PPM motif. R rounds a column to counts
/// and adds any rounding error to the largest count, where the port uses
/// largest remainders, so the two differ for columns that do not round to
/// `nsites` exactly, which is common with twenty letters. Building the PWM
/// R's way lets scans be compared with R's without that difference.
pub fn r_style_pwm(m: &Motif) -> Motif {
    use universalmotif::convert::{pcm_to_ppm, ppm_to_pwm};
    assert_eq!(m.motif_type, MotifType::Ppm);
    let nsites = m.nsites.filter(|&n| n > 1.0).unwrap_or(100.0);
    let pseudo = m.pseudocount;
    let mut p = m.clone().into_parts();
    p.matrix = m.matrix.map_columns(|c| {
        let mut counts: Vec<f64> = c.iter().map(|&x| (x * nsites).round()).collect();
        let s: f64 = counts.iter().sum();
        if s != nsites {
            let mut i = 0;
            for j in 1..counts.len() {
                if counts[j] > counts[i] {
                    i = j;
                }
            }
            counts[i] += nsites - s;
        }
        ppm_to_pwm(&pcm_to_ppm(&counts, pseudo), &m.bkg, 0.0, nsites)
    });
    p.motif_type = MotifType::Pwm;
    Motif::from_parts(p).unwrap()
}
