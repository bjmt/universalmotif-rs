//! The text format of the R package universalmotif, written by its
//! `write_motifs()` and read by its `read_motifs()`.
//!
//! Each motif is a YAML document, after a `# universalmotif version` header
//! line, holding every field of the motif: its names, type, background,
//! pseudocount, `nsites`, P-values, extra information and gaps, with the
//! matrix stored one column per line. Only the plain subset of YAML that the
//! R package writes is read: block mappings and sequences, plain and quoted
//! scalars, and long scalars wrapped onto indented lines.

use super::{fmt_err, format_shortest, join_lines};
use crate::alphabet::Alphabet;
use crate::error::Result;
use crate::matrix::Matrix;
use crate::motif::{GapInfo, Motif, MotifType, Strand};

/// The version of the R package written in the header line. The R package
/// reads files from before 1.1.67 with an older parser, and wrote gaps from
/// 1.31.51, so this is the version whose format is written here.
pub const FORMAT_VERSION: &str = "1.31.51";

/// A parsed YAML node.
#[derive(Debug, Clone, PartialEq)]
enum Node {
    Scalar(String),
    Seq(Vec<Node>),
    Map(Vec<(String, Node)>),
}

impl Node {
    fn scalar(&self) -> Option<&str> {
        match self {
            Node::Scalar(s) => Some(s),
            _ => None,
        }
    }
}

struct Line<'a> {
    indent: usize,
    text: &'a str,
}

/// Read motifs in the R package's universalmotif format.
///
/// Returns the motifs and notes on what could not be kept: higher-order
/// matrices (`multifreq`) and the extra P-values of gapped motifs, which
/// this package does not hold, are left out with a note, as are background
/// entries for words longer than one letter. The total IC and the consensus
/// are recomputed from the matrix.
pub fn read(text: &str) -> Result<(Vec<Motif>, Vec<String>)> {
    let mut lines: Vec<&str> = text.lines().collect();
    let first = lines.iter().find(|l| !l.trim().is_empty());
    if let Some(v) = first.and_then(|l| l.trim().strip_prefix("# universalmotif version"))
        && older_than(v.trim(), [1, 1, 67])
    {
        return fmt_err(format!(
            "universalmotif files from versions before 1.1.67 cannot be read (this one is from \
             {}); read and write it again with a newer version of the R package",
            v.trim()
        ));
    }
    lines.retain(|l| !l.trim_start().starts_with('#'));
    let mut docs: Vec<Vec<&str>> = Vec::new();
    for l in lines {
        if l == "---" || (l.starts_with("---") && !l.starts_with("----")) {
            docs.push(Vec::new());
        } else if let Some(d) = docs.last_mut() {
            d.push(l);
        } else if !l.trim().is_empty() {
            return fmt_err("expected '---' before the first motif");
        }
    }
    if docs.is_empty() {
        return fmt_err("no motifs found");
    }
    let mut motifs = Vec::with_capacity(docs.len());
    let mut notes = Vec::new();
    for (i, d) in docs.iter().enumerate() {
        let lines: Vec<Line> = d
            .iter()
            .filter(|l| !l.trim().is_empty())
            .map(|l| Line {
                indent: l.len() - l.trim_start_matches(' ').len(),
                text: l.trim(),
            })
            .collect();
        if lines.is_empty() {
            continue;
        }
        let mut pos = 0;
        let node = parse_block(&lines, &mut pos, lines[0].indent)?;
        if pos < lines.len() {
            return fmt_err(format!(
                "could not read motif {}: unexpected line '{}'",
                i + 1,
                lines[pos].text
            ));
        }
        let Node::Map(fields) = node else {
            return fmt_err(format!("motif {} is not a set of fields", i + 1));
        };
        motifs.push(motif(&fields, i + 1, &mut notes)?);
    }
    if motifs.is_empty() {
        return fmt_err("no motifs found");
    }
    Ok((motifs, notes))
}

fn older_than(v: &str, than: [u64; 3]) -> bool {
    let parts: Vec<u64> = v
        .split(['.', '-'])
        .map(|p| p.parse().unwrap_or(0))
        .collect();
    for (k, &t) in than.iter().enumerate() {
        let p = parts.get(k).copied().unwrap_or(0);
        if p != t {
            return p < t;
        }
    }
    false
}

fn is_item(text: &str) -> bool {
    text == "-" || text.starts_with("- ")
}

/// Parse the block starting at `lines[*pos]`, whose lines sit at `indent`.
fn parse_block(lines: &[Line], pos: &mut usize, indent: usize) -> Result<Node> {
    if is_item(lines[*pos].text) {
        let mut items = Vec::new();
        while *pos < lines.len() && lines[*pos].indent == indent && is_item(lines[*pos].text) {
            let rest = lines[*pos].text[1..].trim_start();
            *pos += 1;
            if rest.is_empty() {
                if *pos < lines.len() && lines[*pos].indent > indent {
                    let inner = lines[*pos].indent;
                    items.push(parse_block(lines, pos, inner)?);
                } else {
                    items.push(Node::Scalar(String::new()));
                }
            } else {
                items.push(Node::Scalar(scalar(rest, lines, pos, indent)?));
            }
        }
        return Ok(Node::Seq(items));
    }
    let mut fields = Vec::new();
    while *pos < lines.len() && lines[*pos].indent == indent && !is_item(lines[*pos].text) {
        let (key, rest) = split_key(lines[*pos].text)?;
        *pos += 1;
        if !rest.is_empty() {
            fields.push((key, Node::Scalar(scalar(rest, lines, pos, indent)?)));
        } else if *pos < lines.len()
            && (lines[*pos].indent > indent
                || (lines[*pos].indent == indent && is_item(lines[*pos].text)))
        {
            let inner = lines[*pos].indent;
            fields.push((key, parse_block(lines, pos, inner)?));
        } else {
            fields.push((key, Node::Scalar(String::new())));
        }
    }
    Ok(Node::Map(fields))
}

/// Split `key: value` into the unquoted key and the rest.
fn split_key(text: &str) -> Result<(String, &str)> {
    let bad = || fmt_err(format!("expected 'name: value', found '{text}'"));
    let (key, rest) = if text.starts_with('\'') || text.starts_with('"') {
        let q = text.as_bytes()[0];
        let mut end = None;
        let b = text.as_bytes();
        let mut i = 1;
        while i < b.len() {
            if b[i] == q {
                if q == b'\'' && b.get(i + 1) == Some(&b'\'') {
                    i += 2;
                    continue;
                }
                end = Some(i);
                break;
            }
            if q == b'"' && b[i] == b'\\' {
                i += 1;
            }
            i += 1;
        }
        let Some(end) = end else { return bad() };
        let after = &text[end + 1..];
        let Some(rest) = after.strip_prefix(':') else {
            return bad();
        };
        (unquote(&text[..=end])?, rest)
    } else {
        let Some(i) = text
            .find(": ")
            .or_else(|| text.ends_with(':').then(|| text.len() - 1))
        else {
            return bad();
        };
        (text[..i].trim().to_string(), &text[i + 1..])
    };
    Ok((key, rest.trim()))
}

/// A scalar starting with `first`, joined with any lines wrapped below it.
fn scalar(first: &str, lines: &[Line], pos: &mut usize, indent: usize) -> Result<String> {
    let mut text = first.to_string();
    while *pos < lines.len() && lines[*pos].indent > indent {
        text.push(' ');
        text.push_str(lines[*pos].text);
        *pos += 1;
    }
    unquote(&text)
}

fn unquote(s: &str) -> Result<String> {
    let s = s.trim();
    if s.len() >= 2 && s.starts_with('\'') && s.ends_with('\'') {
        return Ok(s[1..s.len() - 1].replace("''", "'"));
    }
    if s.len() >= 2 && s.starts_with('"') && s.ends_with('"') {
        let mut out = String::new();
        let mut chars = s[1..s.len() - 1].chars();
        while let Some(c) = chars.next() {
            if c != '\\' {
                out.push(c);
                continue;
            }
            match chars.next() {
                Some('n') => out.push('\n'),
                Some('t') => out.push('\t'),
                Some('r') => out.push('\r'),
                Some('0') => out.push('\0'),
                Some(' ') => out.push(' '),
                Some('/') => out.push('/'),
                Some(c @ ('x' | 'u' | 'U')) => {
                    let n = match c {
                        'x' => 2,
                        'u' => 4,
                        _ => 8,
                    };
                    let hex: String = chars.by_ref().take(n).collect();
                    match u32::from_str_radix(&hex, 16).ok().and_then(char::from_u32) {
                        Some(ch) => out.push(ch),
                        None => return fmt_err(format!("bad escape in '{s}'")),
                    }
                }
                Some(other) => out.push(other),
                None => return fmt_err(format!("bad escape in '{s}'")),
            }
        }
        return Ok(out);
    }
    Ok(if s == "~" || s == "null" {
        String::new()
    } else {
        s.to_string()
    })
}

fn number(s: &str, what: &str) -> Result<f64> {
    let t = s.trim();
    match t {
        ".inf" | ".Inf" | ".INF" | "+.inf" => return Ok(f64::INFINITY),
        "-.inf" | "-.Inf" | "-.INF" => return Ok(f64::NEG_INFINITY),
        _ => {}
    }
    match t.parse::<f64>() {
        Ok(v) if !v.is_nan() => Ok(v),
        _ => fmt_err(format!("could not read {what} '{s}' as a number")),
    }
}

fn whole(s: &str, what: &str) -> Result<usize> {
    let v = number(s, what)?;
    if v < 0.0 || v.fract() != 0.0 {
        return fmt_err(format!("{what} must be a whole number, not '{s}'"));
    }
    Ok(v as usize)
}

/// A scalar or a sequence of scalars, as a list.
fn list<'a>(n: &'a Node, what: &str) -> Result<Vec<&'a str>> {
    match n {
        Node::Scalar(s) => Ok(vec![s.as_str()]),
        Node::Seq(items) => items
            .iter()
            .map(|i| i.scalar().ok_or(()))
            .collect::<std::result::Result<_, _>>()
            .or_else(|_| fmt_err(format!("{what} must be a list of values"))),
        Node::Map(_) => fmt_err(format!("{what} must be a list of values")),
    }
}

fn motif(fields: &[(String, Node)], i: usize, notes: &mut Vec<String>) -> Result<Motif> {
    let get = |k: &str| fields.iter().find(|(f, _)| f == k).map(|(_, n)| n);
    let text = |k: &str| -> Result<Option<String>> {
        match get(k) {
            None => Ok(None),
            Some(Node::Scalar(s)) if s.is_empty() => Ok(None),
            Some(Node::Scalar(s)) => Ok(Some(s.clone())),
            Some(_) => fmt_err(format!("field '{k}' of motif {i} must be a single value")),
        }
    };
    let num = |k: &str| -> Result<Option<f64>> { text(k)?.map(|s| number(&s, k)).transpose() };
    let Some(name) = text("name")? else {
        return fmt_err(format!("motif {i} has no name"));
    };
    let alphabet: Alphabet = match text("alphabet")? {
        Some(a) => a.parse()?,
        None => return fmt_err(format!("motif '{name}' has no alphabet")),
    };
    let motif_type: MotifType = match text("type")? {
        Some(t) => t.parse()?,
        None => return fmt_err(format!("motif '{name}' has no type")),
    };
    let Some(rows) = get("motif") else {
        return fmt_err(format!("motif '{name}' has no matrix"));
    };
    let letters = alphabet.letters();
    let mut data = Vec::new();
    let columns = list(rows, "motif")?;
    for c in &columns {
        let vals: Vec<f64> = c
            .split_whitespace()
            .map(|v| number(v, "a matrix value"))
            .collect::<Result<_>>()?;
        if vals.len() != letters.len() {
            return fmt_err(format!(
                "motif '{name}' has a column of {} values for {} letters",
                vals.len(),
                letters.len()
            ));
        }
        data.extend(vals);
    }
    let matrix = Matrix::from_vec(letters.len(), columns.len(), data)?;

    let bkg = match get("bkg") {
        None => None,
        Some(Node::Map(entries)) => {
            let mut b = vec![f64::NAN; letters.len()];
            let mut longer = false;
            for (k, v) in entries {
                let mut cs = k.chars();
                match (cs.next(), cs.next()) {
                    (Some(c), None) => {
                        let Some(idx) = alphabet.index_of(c) else {
                            return fmt_err(format!(
                                "motif '{name}' has a background for '{k}', not in its alphabet"
                            ));
                        };
                        b[idx] = number(v.scalar().unwrap_or(""), "a background value")?;
                    }
                    _ => longer = true,
                }
            }
            if b.iter().any(|v| v.is_nan()) {
                return fmt_err(format!(
                    "motif '{name}' lacks a background for every letter"
                ));
            }
            if longer {
                notes.push(format!(
                    "motif '{name}': background frequencies of longer words were left out"
                ));
            }
            Some(b)
        }
        Some(_) => {
            return fmt_err(format!(
                "the background of motif '{name}' must be a mapping"
            ));
        }
    };
    let strand: Strand = match text("strand")? {
        Some(s) => s.parse()?,
        None => Strand::Both,
    };
    let extrainfo = match get("extrainfo") {
        None => Vec::new(),
        Some(Node::Map(entries)) => entries
            .iter()
            .map(|(k, v)| (k.clone(), v.scalar().unwrap_or("").to_string()))
            .collect(),
        Some(Node::Seq(items)) => items
            .iter()
            .map(|v| (String::new(), v.scalar().unwrap_or("").to_string()))
            .collect(),
        Some(Node::Scalar(s)) => vec![(String::new(), s.clone())],
    };
    let isgapped = matches!(
        text("isgapped")?
            .as_deref()
            .map(str::to_ascii_lowercase)
            .as_deref(),
        Some("yes" | "true" | "y" | "on")
    );
    let gap = if isgapped {
        let field = |k: &str| -> Result<Vec<usize>> {
            match get(k) {
                Some(n) => list(n, k)?.iter().map(|s| whole(s, k)).collect(),
                None => fmt_err(format!("gapped motif '{name}' has no {k}")),
            }
        };
        GapInfo {
            isgapped: true,
            gaploc: field("gaploc")?,
            mingap: field("mingap")?,
            maxgap: field("maxgap")?,
        }
    } else {
        GapInfo::default()
    };
    if get("multifreq").is_some() {
        notes.push(format!(
            "motif '{name}': its higher-order matrices (multifreq) were left out"
        ));
    }
    if get("extrapvals").is_some() {
        notes.push(format!(
            "motif '{name}': the extra P-values of its gaps were left out"
        ));
    }
    let nsites = num("nsites")?;
    let mut b = Motif::builder(matrix, alphabet)
        .name(name)
        .altname(text("altname")?)
        .family(text("family")?)
        .organism(text("organism")?)
        .motif_type(Some(motif_type))
        .nsites(nsites)
        .bkg(bkg)
        .bkgsites(num("bkgsites")?)
        .strand(strand)
        .pval(num("pval")?)
        .qval(num("qval")?)
        .eval(num("eval")?)
        .extrainfo(extrainfo)
        .gap(gap);
    if let Some(p) = num("pseudocount")? {
        b = b.pseudocount(p);
    }
    b.build()
}

/// Write motifs in the R package's universalmotif format, with a header line
/// unless `header` is false (as when adding to a file). Every field is
/// written, with numbers given in full so that they read back exactly.
pub fn write(motifs: &[Motif], header: bool) -> Result<String> {
    let mut lines = Vec::new();
    if header {
        lines.push(format!("# universalmotif version {FORMAT_VERSION}"));
    }
    for m in motifs {
        lines.push(String::new());
        lines.push("---".into());
        lines.push(format!("name: {}", quote(&m.name)));
        for (k, v) in [
            ("altname", &m.altname),
            ("family", &m.family),
            ("organism", &m.organism),
        ] {
            if let Some(v) = v {
                lines.push(format!("{k}: {}", quote(v)));
            }
        }
        lines.push(format!("alphabet: {}", quote(&m.alphabet.name())));
        lines.push(format!("type: {}", m.motif_type.as_str()));
        if m.icscore.is_finite() {
            lines.push(format!("icscore: {}", float(m.icscore)));
        }
        if let Some(n) = m.nsites {
            lines.push(format!("nsites: {}", float(n)));
        }
        lines.push(format!("pseudocount: {}", float(m.pseudocount)));
        lines.push("bkg:".into());
        for (l, v) in m.alphabet.letters().iter().zip(&m.bkg) {
            lines.push(format!("  {}: {}", quote(&l.to_string()), float(*v)));
        }
        if let Some(n) = m.bkgsites {
            lines.push(format!("bkgsites: {}", float(n)));
        }
        if let Some(c) = &m.consensus {
            lines.push(format!("consensus: {}", quote(c)));
        }
        lines.push(format!("strand: {}", quote(m.strand.as_str())));
        for (k, v) in [("pval", m.pval), ("qval", m.qval), ("eval", m.eval)] {
            if let Some(v) = v {
                lines.push(format!("{k}: {}", float(v)));
            }
        }
        if !m.extrainfo.is_empty() {
            lines.push("extrainfo:".into());
            for (k, v) in &m.extrainfo {
                lines.push(format!("  {}: {}", quote(k), quote(v)));
            }
        }
        if m.gap.isgapped {
            lines.push("isgapped: yes".into());
            for (k, v) in [
                ("gaploc", &m.gap.gaploc),
                ("mingap", &m.gap.mingap),
                ("maxgap", &m.gap.maxgap),
            ] {
                if v.len() == 1 {
                    lines.push(format!("{k}: {}", v[0]));
                } else {
                    lines.push(format!("{k}:"));
                    lines.extend(v.iter().map(|x| format!("- {x}")));
                }
            }
        }
        lines.push("motif:".into());
        for c in m.matrix.columns() {
            let vals: Vec<String> = c.iter().map(|v| format_shortest(*v)).collect();
            lines.push(format!("- {}", vals.join(" ")));
        }
    }
    Ok(join_lines(&lines))
}

/// A number as YAML reads it as a number: with a decimal point, as a plain
/// `1e-05` would be read as text.
fn float(x: f64) -> String {
    if x.is_nan() {
        return ".nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 {
            ".inf".into()
        } else {
            "-.inf".into()
        };
    }
    let s = format_shortest(x);
    if s.contains('.') {
        return s;
    }
    match s.split_once('e') {
        Some((m, e)) => format!("{m}.0e{e}"),
        None => format!("{s}.0"),
    }
}

/// A string as a YAML scalar, quoted whenever it would otherwise be read as
/// something else.
fn quote(s: &str) -> String {
    if s.chars().any(|c| c.is_control()) {
        let mut out = String::from('"');
        for c in s.chars() {
            match c {
                '"' => out.push_str("\\\""),
                '\\' => out.push_str("\\\\"),
                '\n' => out.push_str("\\n"),
                '\t' => out.push_str("\\t"),
                '\r' => out.push_str("\\r"),
                c if c.is_control() => out.push_str(&format!("\\u{:04x}", c as u32)),
                c => out.push(c),
            }
        }
        out.push('"');
        return out;
    }
    let lower = s.to_ascii_lowercase();
    let special = s.is_empty()
        || s != s.trim()
        || s.contains(": ")
        || s.contains(" #")
        || s.ends_with(':')
        || s.starts_with(|c: char| "-?:,[]{}#&*!|>'\"%@`".contains(c))
        || matches!(
            lower.as_str(),
            "y" | "n" | "yes" | "no" | "true" | "false" | "on" | "off" | "null" | "~"
        )
        || lower.trim_start_matches(['+', '-']).starts_with('.')
            && matches!(lower.trim_start_matches(['+', '-']), ".inf" | ".nan")
        || s.parse::<f64>().is_ok();
    if special {
        format!("'{}'", s.replace('\'', "''"))
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quoting() {
        assert_eq!(quote("DR"), "DR");
        assert_eq!(quote("MA0108.2 TBP"), "MA0108.2 TBP");
        assert_eq!(quote("-"), "'-'");
        assert_eq!(quote("+-"), "+-");
        assert_eq!(quote("N"), "'N'");
        assert_eq!(quote("yes"), "'yes'");
        assert_eq!(quote("it's"), "it's");
        assert_eq!(quote("'x'"), "'''x'''");
        assert_eq!(quote("a: b"), "'a: b'");
        assert_eq!(quote("1.5"), "'1.5'");
        assert_eq!(quote("a\tb"), "\"a\\tb\"");
        assert_eq!(float(50.0), "50.0");
        assert_eq!(float(1e-5), "1.0e-05");
        assert_eq!(float(0.25), "0.25");
        assert_eq!(float(f64::INFINITY), ".inf");
    }

    #[test]
    fn versions() {
        assert!(older_than("1.1.66", [1, 1, 67]));
        assert!(!older_than("1.1.67", [1, 1, 67]));
        assert!(!older_than("1.31.66", [1, 1, 67]));
        assert!(older_than("0.6.0", [1, 1, 67]));
    }
}
