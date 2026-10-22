//! Reading and writing motif files.
//!
//! Each format has its own module with `read` and `write` functions that work
//! on text, so callers decide how files are opened, skipped or appended to.

pub mod cisbp;
pub mod fasta;
pub mod homer;
pub mod jaspar;
pub mod matrix;
pub mod meme;
pub mod transfac;
pub mod uniprobe;

use crate::error::{Error, Result};

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
pub fn format_r_number(x: f64) -> String {
    let target = format!("{x:.14e}").parse().unwrap_or(x);
    format_significant(x, 15, target)
}

/// Format a number with the fewest digits that read back as exactly the
/// same value, in scientific notation only when it is shorter than fixed
/// notation.
pub fn format_shortest(x: f64) -> String {
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
