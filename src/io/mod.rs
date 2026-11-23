//! Reading and writing motif and FASTA files.
//!
//! Each format has its own module with `read` and `write` functions that work
//! on text, so callers decide how files are opened, skipped or appended to;
//! FASTA files can also be read from a path, gzipped or not, and MEME output
//! has readers for the sites and statistics it holds.

pub mod cisbp;
pub mod fasta;
pub mod homer;
pub mod jaspar;
pub mod matrix;
pub mod meme;
pub mod transfac;
pub mod uniprobe;
pub mod universalmotif;

use crate::error::{Error, Result};
use crate::motif::{GapInfo, Motif};

/// A gapped motif's gaps written into its consensus, as in
/// `AGGTCA[1-5]AGGTCA`: each gap's range of sizes, in brackets, after the
/// letters of the columns before it. `None` for motifs without gaps, or
/// without a consensus (those of custom alphabets).
///
/// The MEME, JASPAR and TRANSFAC writers give gapped motifs this as their
/// alternative name, and the readers restore the gaps from it, as no motif
/// format has a place for them.
pub fn gap_notation(m: &Motif) -> Option<String> {
    if !m.gap.isgapped {
        return None;
    }
    let c: Vec<char> = m.consensus.as_ref()?.chars().collect();
    let mut out = String::new();
    let mut at = 0;
    for ((&loc, &lo), &hi) in m.gap.gaploc.iter().zip(&m.gap.mingap).zip(&m.gap.maxgap) {
        let loc = loc.min(c.len());
        out.extend(&c[at..loc]);
        out.push_str(&format!("[{lo}-{hi}]"));
        at = loc;
    }
    out.extend(&c[at..]);
    Some(out)
}

/// The gaps that `text` gives in the notation of [`gap_notation`], if it is
/// exactly that notation and its letters number `width`.
pub fn parse_gap_notation(text: &str, width: usize) -> Option<GapInfo> {
    let mut gap = GapInfo {
        isgapped: true,
        ..GapInfo::default()
    };
    let mut letters = 0;
    let mut rest = text;
    while !rest.is_empty() {
        let run = rest.find('[').unwrap_or(rest.len());
        if !rest[..run].chars().all(|c| c.is_ascii_alphabetic()) {
            return None;
        }
        letters += run;
        rest = &rest[run..];
        if rest.is_empty() {
            break;
        }
        let close = rest.find(']')?;
        let (lo, hi) = rest[1..close].split_once('-')?;
        let digits = |s: &str| !s.is_empty() && s.chars().all(|c| c.is_ascii_digit());
        if run == 0 || !digits(lo) || !digits(hi) {
            return None;
        }
        gap.gaploc.push(letters);
        gap.mingap.push(lo.parse().ok()?);
        gap.maxgap.push(hi.parse().ok()?);
        rest = &rest[close + 1..];
    }
    (!gap.gaploc.is_empty() && letters == width).then_some(gap)
}

/// The alternative name a motif file gets: the gap notation for a gapped
/// motif, and the motif's own alternative name otherwise, unless that name
/// is itself a gap notation that fits the motif, as when a motif read with
/// gaps has had them removed, since reading it back would restore them.
pub(crate) fn file_altname(m: &Motif) -> Option<String> {
    gap_notation(m).or_else(|| {
        m.altname
            .clone()
            .filter(|a| parse_gap_notation(a, m.width()).is_none())
    })
}

/// A motif read from a file, with the gaps its alternative name gives if
/// that name is in the notation of [`gap_notation`]. Gaps that do not fit
/// the motif leave it as it is.
pub(crate) fn with_noted_gaps(m: Motif) -> Motif {
    let Some(gap) = m
        .altname
        .as_deref()
        .filter(|_| !m.gap.isgapped)
        .and_then(|a| parse_gap_notation(a, m.width()))
    else {
        return m;
    };
    let mut parts = m.clone().into_parts();
    parts.gap = gap;
    Motif::from_parts(parts).unwrap_or(m)
}

/// The lines with any content; lines holding only whitespace count as
/// empty.
pub(crate) fn content_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.trim().is_empty()).collect()
}

/// Parse a number as found in motif files: surrounding whitespace is ignored,
/// `Inf` and `-Inf` are accepted, and `NA` or anything unparsable gives `None`.
pub(crate) fn parse_num(s: &str) -> Option<f64> {
    let s = s.trim();
    if s == "NA" {
        return None;
    }
    s.parse::<f64>().ok().filter(|v| !v.is_nan())
}

pub(crate) fn num(s: &str, what: &str) -> Result<f64> {
    parse_num(s).ok_or_else(|| Error::Argument(format!("could not read {what} '{s}' as a number")))
}

pub(crate) fn fmt_err<T>(msg: impl Into<String>) -> Result<T> {
    Err(Error::Argument(msg.into()))
}

/// Format a number with at most 15 significant digits, with scientific notation used when it is
/// shorter than fixed notation.
pub(crate) fn format_r_number(x: f64) -> String {
    let target = format!("{x:.14e}").parse().unwrap_or(x);
    format_significant(x, 15, target)
}

/// Format a number with the fewest digits that read back as exactly the
/// same value, in scientific notation only when it is shorter than fixed
/// notation.
pub(crate) fn format_shortest(x: f64) -> String {
    format_significant(x, 17, x)
}

/// Use the fewest significant digits, at most `max_digits`, that give `target`
/// when parsed, in scientific notation only when it is shorter than fixed
/// notation.
fn format_significant(x: f64, max_digits: usize, target: f64) -> String {
    if x.is_nan() {
        return "NaN".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    if x == 0.0 {
        return "0".into();
    }
    let mut sci = String::new();
    for d in 1..=max_digits {
        sci = format!("{:.*e}", d - 1, x);
        if sci.parse::<f64>().unwrap() == target {
            break;
        }
    }
    let (mant, exp) = sci.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let neg = mant.starts_with('-');
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    let sign = if neg { "-" } else { "" };

    let sci_str = format!(
        "{sign}{}{}e{}{:02}",
        &digits[..1],
        if digits.len() > 1 {
            format!(".{}", &digits[1..])
        } else {
            String::new()
        },
        if exp < 0 { '-' } else { '+' },
        exp.abs()
    );
    // Fixed notation prints the integer part exactly, as sprintf("%.*f") does.
    let decimals = if exp >= 0 {
        digits.len().saturating_sub(exp as usize + 1)
    } else {
        (-exp - 1) as usize + digits.len()
    };
    let fixed_str = format!("{x:.decimals$}");
    if fixed_str.len() <= sci_str.len() {
        fixed_str
    } else {
        sci_str
    }
}

/// Fixed notation with `digits` decimals, writing non-finite values as
/// `NaN`, `Inf` and `-Inf`.
pub(crate) fn format_fixed(x: f64, digits: usize) -> String {
    if x.is_nan() {
        "NaN".into()
    } else if x.is_infinite() {
        if x > 0.0 { "Inf".into() } else { "-Inf".into() }
    } else {
        format!("{x:.digits$}")
    }
}

/// Join lines, each one followed by a newline.
pub(crate) fn join_lines(lines: &[String]) -> String {
    let mut out = String::with_capacity(lines.iter().map(|l| l.len() + 1).sum());
    for l in lines {
        out.push_str(l);
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::format_r_number as f;

    #[test]
    fn r_number_formatting() {
        assert_eq!(f(30.0), "30");
        assert_eq!(f(0.5), "0.5");
        assert_eq!(f(1.0 / 3.0), "0.333333333333333");
        assert_eq!(f(1e5), "1e+05");
        assert_eq!(f(123456.0), "123456");
        assert_eq!(f(1e-4), "1e-04");
        assert_eq!(f(0.001), "0.001");
        assert_eq!(f(0.00012), "0.00012");
        assert_eq!(f(-2.5e-12), "-2.5e-12");
        assert_eq!(f(4.70742403722102), "4.70742403722102");
        assert_eq!(f(1e15), "1e+15");
        assert_eq!(f(1.5e300), "1.5e+300");
        assert_eq!(f(f64::NEG_INFINITY), "-Inf");
        assert_eq!(f(100.0), "100");
        assert_eq!(f(1234567890123456.0), "1234567890123456");
        assert_eq!(f(0.1 + 0.2), "0.3");
        assert_eq!(f(2.0 / 3.0 * 1e-7), "6.66666666666667e-08");
    }

    #[test]
    fn shortest_formatting_round_trips() {
        use super::format_shortest as s;
        assert_eq!(s(0.1), "0.1");
        assert_eq!(s(1.0 / 3.0), "0.3333333333333333");
        assert_eq!(s(1e-20), "1e-20");
        assert_eq!(s(25.0), "25");
        for x in [
            1.0 / 3.0,
            0.17647058823529413,
            2.5e-300,
            123456.789,
            -7.0 / 9.0,
        ] {
            assert_eq!(s(x).parse::<f64>().unwrap(), x);
        }
    }
}
