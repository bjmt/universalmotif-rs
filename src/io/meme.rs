//! MEME motif files.
//!
//! Both full MEME/DREME/STREME output and the minimal MEME motif format are
//! read. Besides the motifs, full MEME output can be mined for the sites each
//! motif was built from with [`read_sites`], the per-site P-values with
//! [`read_site_pvalues`] and the combined block diagrams with
//! [`read_combined`].

use super::{content_lines, fmt_err, format_fixed, format_r_number, join_lines, num, parse_num};
use crate::alphabet::Alphabet;
use crate::convert::convert_type;
use crate::error::{Error, Result};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType, Strand};

/// The order MEME uses for the letters of custom alphabets.
const MEME_ORDER: &str = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789*-.";

fn meme_rank(c: char) -> usize {
    MEME_ORDER.find(c).unwrap_or(usize::MAX)
}

/// Content lines without comments and the `****` / `----` rulers of MEME
/// output.
fn meme_lines(text: &str) -> Vec<&str> {
    content_lines(text)
        .into_iter()
        .filter(|l| !l.starts_with('#') && !l.contains("**********") && !l.contains("------------"))
        .collect()
}

fn alphabet_from_letters(letters: &str) -> Result<Alphabet> {
    match letters {
        "ACGT" => Ok(Alphabet::Dna),
        "ACGU" => Ok(Alphabet::Rna),
        "ACDEFGHIKLMNPQRSTVWY" => Ok(Alphabet::Aa),
        _ => Alphabet::custom(letters.chars()),
    }
}

struct Header {
    alphabet: Alphabet,
    /// Letters in the order the file lists them in matrices.
    file_letters: Vec<char>,
    strand: Strand,
    bkg: Vec<f64>,
}

fn parse_header(lines: &[&str]) -> Result<Header> {
    let (alphabet, file_letters) =
        if let Some(l) = lines.iter().find(|l| l.starts_with("ALPHABET=")) {
            let letters = l["ALPHABET=".len()..]
                .split_whitespace()
                .next()
                .unwrap_or("");
            (
                alphabet_from_letters(letters)?,
                letters.chars().collect::<Vec<_>>(),
            )
        } else if let Some(i) = lines.iter().position(|l| l.starts_with("ALPHABET")) {
            let end = lines
                .iter()
                .position(|l| l.starts_with("END ALPHABET"))
                .or_else(|| lines.iter().position(|l| l.starts_with("TRAINING SET")))
                .ok_or_else(|| {
                    Error::Argument("could not find the end of the custom alphabet".into())
                })?;
            let mut letters: Vec<char> = lines[i + 1..end.max(i + 1)]
                .iter()
                .filter(|l| !l.contains('='))
                .flat_map(|l| l.split('~'))
                .filter_map(|p| p.split_whitespace().next())
                .filter_map(|t| t.chars().next())
                .collect();
            letters.sort_by_key(|&c| meme_rank(c));
            let name: String = letters.iter().collect();
            (alphabet_from_letters(&name)?, letters)
        } else {
            (Alphabet::Dna, crate::alphabet::DNA.to_vec())
        };
    if file_letters.len() != alphabet.len() {
        return fmt_err("could not parse the MEME alphabet");
    }

    let strand = match lines.iter().find(|l| l.starts_with("strands:")) {
        Some(l) => {
            let s: Vec<&str> = l.split_whitespace().skip(1).collect();
            if s.contains(&"+") && s.contains(&"-") {
                Strand::Both
            } else {
                Strand::parse(s.first().copied().unwrap_or("+"))?
            }
        }
        None => Strand::Plus,
    };

    let n = alphabet.len();
    let bkg = match lines
        .iter()
        .position(|l| l.starts_with("Background letter frequencies"))
    {
        None => vec![1.0 / n as f64; n],
        Some(i) => {
            let mut pairs: Vec<(char, f64)> = Vec::new();
            for l in lines.iter().skip(i + 1) {
                let t: Vec<&str> = l.split_whitespace().collect();
                for p in t.chunks(2) {
                    let (Some(letter), Some(v)) = (p[0].chars().next(), p.get(1)) else {
                        return fmt_err("could not parse background frequencies");
                    };
                    pairs.push((letter, num(v, "a background frequency")?));
                }
                if pairs.len() >= n {
                    break;
                }
            }
            if pairs.len() < n {
                return fmt_err(
                    "could not parse background frequencies, check that they match the alphabet",
                );
            }
            let mut out = vec![f64::NAN; n];
            if pairs.iter().all(|(l, _)| alphabet.index_of(*l).is_some()) {
                for (l, v) in &pairs {
                    out[alphabet.index_of(*l).unwrap()] = *v;
                }
            }
            if out.iter().any(|v| v.is_nan()) {
                // Unlabelled or unexpected letters: take the values in order.
                out = pairs.iter().take(n).map(|p| p.1).collect();
            }
            out
        }
    };
    Ok(Header {
        alphabet,
        file_letters,
        strand,
        bkg,
    })
}

/// Read `key= value` or `key=value` from a `letter-probability matrix:` line.
fn matrix_field<'a>(tokens: &[&'a str], key: &str) -> Option<&'a str> {
    let with_eq = format!("{key}=");
    for (i, t) in tokens.iter().enumerate() {
        if *t == with_eq {
            return tokens.get(i + 1).copied();
        }
        if let Some(v) = t.strip_prefix(&with_eq)
            && !v.is_empty()
        {
            return Some(v);
        }
    }
    None
}

struct RawMotif<'a> {
    name: &'a str,
    altname: Option<&'a str>,
    width: usize,
    nsites: Option<f64>,
    eval: Option<&'a str>,
    rows: &'a [&'a str],
}

fn raw_motifs<'a>(lines: &'a [&'a str]) -> Result<Vec<RawMotif<'a>>> {
    let motif_lines: Vec<usize> = (0..lines.len())
        .filter(|&i| lines[i].starts_with("MOTIF "))
        .collect();
    if motif_lines.is_empty() {
        return fmt_err("no motifs found");
    }
    let mut out = Vec::new();
    for p in (0..lines.len()).filter(|&i| lines[i].starts_with("letter-probability matrix:")) {
        let Some(&m) = motif_lines.iter().rev().find(|&&m| m < p) else {
            return fmt_err("found a letter-probability matrix without a MOTIF line");
        };
        let mt: Vec<&str> = lines[m].split_whitespace().collect();
        let name = mt.get(1).copied().unwrap_or("");
        let altname = mt.get(2).copied().filter(|a| *a != "width");
        let pt: Vec<&str> = lines[p].split_whitespace().collect();
        let width = matrix_field(&pt, "w")
            .and_then(parse_num)
            .filter(|w| *w >= 1.0 && w.fract() == 0.0)
            .ok_or_else(|| Error::Argument(format!("motif '{name}' has no valid width (w=)")))?
            as usize;
        if p + 1 + width > lines.len() {
            return fmt_err(format!("motif '{name}' is truncated"));
        }
        out.push(RawMotif {
            name,
            altname,
            width,
            nsites: matrix_field(&pt, "nsites").and_then(parse_num),
            eval: matrix_field(&pt, "E"),
            rows: &lines[p + 1..p + 1 + width],
        });
    }
    if out.is_empty() {
        return fmt_err("no motifs found");
    }
    Ok(out)
}

/// Options for [`read()`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ReadOptions {
    /// Read the matrices as contribution weight matrices, as written by
    /// TF-MoDISco, rather than probabilities.
    pub cwm: bool,
}

/// Read MEME-formatted motifs.
///
/// # Example
///
/// ```
/// use universalmotif::io::meme;
///
/// let text = "MEME version 4
///
/// ALPHABET= ACGT
///
/// strands: + -
///
/// Background letter frequencies
/// A 0.25 C 0.25 G 0.25 T 0.25
///
/// MOTIF m1
/// letter-probability matrix: alength= 4 w= 2 nsites= 10
///  0.9 0.0 0.1 0.0
///  0.0 0.0 0.0 1.0
/// ";
/// let motifs = meme::read(text, &meme::ReadOptions::default())?;
/// assert_eq!((motifs[0].name.as_str(), motifs[0].width()), ("m1", 2));
/// let again = meme::write(&motifs, &meme::WriteOptions::default())?;
/// assert!(again.contains("MOTIF m1"));
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn read(text: &str, opts: &ReadOptions) -> Result<Vec<Motif>> {
    let lines = meme_lines(text);
    let h = parse_header(&lines)?;
    let n = h.alphabet.len();
    let order: Vec<usize> = h
        .alphabet
        .letters()
        .iter()
        .map(|l| h.file_letters.iter().position(|f| f == l).unwrap())
        .collect();
    let mut motifs = Vec::new();
    for r in raw_motifs(&lines)? {
        let mut data = Vec::with_capacity(n * r.width);
        for row in r.rows {
            let vals: Vec<f64> = row
                .split_whitespace()
                .map(|v| num(v, "a probability"))
                .collect::<Result<_>>()?;
            if vals.len() != n {
                return fmt_err(format!(
                    "motif '{}': alphabet length does not match motif length",
                    r.name
                ));
            }
            data.extend(order.iter().map(|&i| vals[i]));
        }
        let matrix = Matrix::from_vec(n, r.width, data)?;
        let extrainfo = r
            .eval
            .map(|e| vec![("eval.string".to_string(), e.to_string())])
            .unwrap_or_default();
        motifs.push(
            Motif::builder(matrix, h.alphabet.clone())
                .name(r.name)
                .altname(r.altname.map(String::from))
                .motif_type(Some(if opts.cwm {
                    MotifType::Cwm
                } else {
                    MotifType::Ppm
                }))
                .nsites(r.nsites)
                .eval(r.eval.and_then(parse_num))
                .bkg(Some(h.bkg.clone()))
                .strand(h.strand)
                .extrainfo(extrainfo)
                .build()?,
        );
    }
    Ok(motifs)
}

/// A site from the BLOCKS section of MEME output.
#[derive(Debug, Clone, PartialEq)]
pub struct Site {
    pub sequence: String,
    pub site: String,
}

/// The sites one motif was built from.
#[derive(Debug, Clone, PartialEq)]
pub struct MotifSites {
    pub motif: String,
    pub sites: Vec<Site>,
}

/// Read the sites listed in the BLOCKS sections of MEME output.
pub fn read_sites(text: &str) -> Result<Vec<MotifSites>> {
    let lines = meme_lines(text);
    let mut out = Vec::new();
    for i in (0..lines.len()).filter(|&i| lines[i].contains("in BLOCKS format")) {
        let Some(bl) = lines.get(i + 1) else {
            return fmt_err("BLOCKS section is truncated");
        };
        let t: Vec<&str> = bl.split_whitespace().collect();
        let motif = t.get(2).copied().unwrap_or("").to_string();
        let nseq = t
            .iter()
            .find_map(|x| x.strip_prefix("seqs="))
            .and_then(parse_num)
            .ok_or_else(|| Error::Argument(format!("BLOCKS section for '{motif}' has no seqs=")))?
            as usize;
        let rows = lines
            .get(i + 2..i + 2 + nseq)
            .ok_or_else(|| Error::Argument(format!("BLOCKS section for '{motif}' is truncated")))?;
        let sites = rows
            .iter()
            .map(|r| {
                let t: Vec<&str> = r.split_whitespace().collect();
                let k = t.iter().position(|x| x.ends_with(')'));
                match (t.first(), k.and_then(|k| t.get(k + 1))) {
                    (Some(s), Some(site)) => Ok(Site {
                        sequence: s.to_string(),
                        site: site.to_string(),
                    }),
                    _ => fmt_err(format!("could not read BLOCKS line '{r}'")),
                }
            })
            .collect::<Result<_>>()?;
        out.push(MotifSites { motif, sites });
    }
    Ok(out)
}

/// A site with its P-value, from the "sites sorted by position p-value"
/// tables of MEME output.
#[derive(Debug, Clone, PartialEq)]
pub struct SitePvalue {
    pub sequence: String,
    /// `'+'` or `'-'`, when the table has a strand column.
    pub strand: Option<char>,
    pub position: f64,
    pub pvalue: f64,
    pub site: String,
}

/// The site P-value table of one motif.
#[derive(Debug, Clone, PartialEq)]
pub struct MotifSitePvalues {
    pub motif: String,
    pub sites: Vec<SitePvalue>,
}

/// Read the per-site P-value tables of MEME output.
pub fn read_site_pvalues(text: &str) -> Result<Vec<MotifSitePvalues>> {
    let lines = meme_lines(text);
    let widths: Vec<(&str, usize)> = raw_motifs(&lines)
        .map(|r| r.iter().map(|m| (m.name, m.width)).collect())
        .unwrap_or_default();
    let mut out = Vec::new();
    for i in (0..lines.len()).filter(|&i| lines[i].contains("sites sorted by position p-value")) {
        let motif = lines[i]
            .split_whitespace()
            .skip_while(|t| *t != "Motif")
            .nth(1)
            .unwrap_or("")
            .to_string();
        let width = widths.iter().find(|(n, _)| *n == motif).map(|w| w.1);
        let has_strand = lines.get(i + 1).is_some_and(|h| h.contains("Strand"));
        let mut sites = Vec::new();
        for r in lines.iter().skip(i + 2) {
            if r.ends_with("block diagrams") || r.contains("in BLOCKS format") {
                break;
            }
            let t: Vec<&str> = r.split_whitespace().collect();
            let k = if has_strand { 2 } else { 1 };
            if t.len() < k + 3 {
                return fmt_err(format!("could not read site line '{r}'"));
            }
            let rest = &t[k + 2..];
            let site = if rest.len() == 3 {
                rest[1]
            } else {
                rest.iter()
                    .find(|s| Some(s.len()) == width)
                    .copied()
                    .unwrap_or(rest[rest.len() / 2])
            };
            sites.push(SitePvalue {
                sequence: t[0].to_string(),
                strand: if has_strand {
                    t[1].chars().next()
                } else {
                    None
                },
                position: num(t[k], "a site position")?,
                pvalue: num(t[k + 1], "a site P-value")?,
                site: site.to_string(),
            });
        }
        out.push(MotifSitePvalues { motif, sites });
    }
    Ok(out)
}

/// A row of the combined block diagrams of MEME output.
#[derive(Debug, Clone, PartialEq)]
pub struct CombinedSite {
    pub sequence: String,
    pub pvalue: f64,
    pub diagram: String,
}

/// Read the combined block diagrams of MEME output, which give the
/// non-overlapping sites in each sequence. Diagrams continued over several
/// lines are joined.
pub fn read_combined(text: &str) -> Result<Vec<CombinedSite>> {
    let lines = meme_lines(text);
    let Some(start) = lines
        .iter()
        .position(|l| l.contains("Combined block diagrams: non-overlapping sites"))
    else {
        return Ok(Vec::new());
    };
    let mut out: Vec<CombinedSite> = Vec::new();
    let mut continuing = false;
    for l in lines.iter().skip(start + 2) {
        if continuing {
            let last = out.last_mut().unwrap();
            last.diagram.push_str(l.trim());
        } else {
            let t: Vec<&str> = l.split_whitespace().collect();
            let pvalue = t.get(1).and_then(|p| parse_num(p));
            match (t.len(), pvalue) {
                (3, Some(p)) => out.push(CombinedSite {
                    sequence: t[0].to_string(),
                    pvalue: p,
                    diagram: t[2].to_string(),
                }),
                _ => break,
            }
        }
        let last = out.last_mut().unwrap();
        continuing = last.diagram.ends_with('\\');
        if continuing {
            last.diagram.pop();
        }
    }
    Ok(out)
}

/// Options for [`write()`].
#[derive(Debug, Clone, PartialEq)]
pub struct WriteOptions {
    pub version: u32,
    /// Background frequencies in the order letters are written, which for DNA
    /// is A, C, G, T. By default the motifs' shared background is used, or a
    /// uniform one if they differ.
    pub bkg: Option<Vec<f64>>,
    /// The `strands:` value. By default it comes from the motifs.
    pub strand: Option<String>,
    /// Write contribution weight matrices as they are instead of converting
    /// to probabilities. Every motif must then be a CWM.
    pub cwm: bool,
    /// Write the file header. Turn this off to append motifs to an existing
    /// file.
    pub header: bool,
}

impl Default for WriteOptions {
    fn default() -> Self {
        WriteOptions {
            version: 5,
            bkg: None,
            strand: None,
            cwm: false,
            header: true,
        }
    }
}

/// Replace characters MEME does not allow in motif names.
fn meme_name(name: &str) -> String {
    let name = name.replace('=', ".");
    let mut parts: Vec<&str> = name.split(' ').collect();
    if parts.len() > 1 && parts.last() == Some(&"") {
        parts.pop();
    }
    parts.join("-")
}

/// Write motifs in MEME format. All motifs must share an alphabet.
pub fn write(motifs: &[Motif], opts: &WriteOptions) -> Result<String> {
    let Some(first) = motifs.first() else {
        return fmt_err("no motifs to write");
    };
    let alphabet = first.alphabet.clone();
    if motifs.iter().any(|m| m.alphabet != alphabet) {
        return fmt_err("all motifs must use the same alphabet");
    }
    let motifs: Vec<Motif> = if opts.cwm {
        if motifs.iter().any(|m| m.motif_type != MotifType::Cwm) {
            return fmt_err("cwm = true requires every motif to be of type CWM");
        }
        motifs.to_vec()
    } else {
        motifs
            .iter()
            .map(|m| convert_type(m, MotifType::Ppm, None, false))
            .collect::<Result<_>>()?
    };

    let mut letters: Vec<char> = alphabet.letters().to_vec();
    let alph_lines: Vec<String> = match &alphabet {
        Alphabet::Dna => vec!["ALPHABET= ACGT".into()],
        Alphabet::Rna => vec!["ALPHABET= ACGU".into()],
        Alphabet::Aa => vec!["ALPHABET= ACDEFGHIKLMNPQRSTVWY".into()],
        Alphabet::Custom(_) => {
            if letters.iter().any(|&c| meme_rank(c) == usize::MAX) {
                return fmt_err(format!(
                    "MEME format motifs only allow the following alphabet symbols: {MEME_ORDER}"
                ));
            }
            letters.sort_by_key(|&c| meme_rank(c));
            let mut v = vec![format!(
                "ALPHABET \"{}\"",
                letters.iter().collect::<String>()
            )];
            v.extend(letters.iter().map(|c| c.to_string()));
            v.push("END ALPHABET".into());
            v
        }
    };
    let order: Vec<usize> = letters
        .iter()
        .map(|l| alphabet.index_of(*l).unwrap())
        .collect();

    let mut lines: Vec<String> = Vec::new();
    if opts.header {
        let strand = match &opts.strand {
            Some(s) => s.clone(),
            None => {
                let s = motifs[0].strand;
                match s {
                    Strand::Plus | Strand::Minus if motifs.iter().all(|m| m.strand == s) => {
                        s.as_str().into()
                    }
                    _ => "+ -".into(),
                }
            }
        };
        let bkg = match &opts.bkg {
            Some(b) => {
                if b.len() != letters.len() {
                    return fmt_err(format!("bkg must have {} values", letters.len()));
                }
                b.clone()
            }
            None => {
                let per: Vec<Vec<f64>> = motifs
                    .iter()
                    .map(|m| order.iter().map(|&i| m.bkg[i]).collect())
                    .collect();
                if per.iter().all(|b| *b == per[0]) {
                    per[0].clone()
                } else {
                    vec![1.0 / letters.len() as f64; letters.len()]
                }
            }
        };
        lines.push(format!("MEME version {}", opts.version));
        lines.push(String::new());
        lines.extend(alph_lines);
        lines.push(String::new());
        lines.push(format!("strands: {strand}"));
        lines.push(String::new());
        lines.push("Background letter frequencies".into());
        let pairs: Vec<String> = letters
            .iter()
            .zip(&bkg)
            .map(|(l, b)| format!("{l} {}", format_fixed(*b, 6)))
            .collect();
        lines.push(pairs.join(" "));
        lines.push(String::new());
    }

    for m in &motifs {
        let mut head = format!("MOTIF {}", meme_name(&m.name));
        if let Some(alt) = &m.altname {
            head.push(' ');
            head.push_str(&meme_name(alt));
        }
        lines.push(head);
        // The E-value is optional in MEME files, so it is left out rather
        // than written as 0 when the motif has none.
        let mut info = format!(
            "letter-probability matrix: alength= {} w= {} nsites= {}",
            letters.len(),
            m.width(),
            format_r_number(m.nsites.unwrap_or(100.0)),
        );
        if let Some(e) = m.eval {
            info.push_str(&format!(" E= {}", format_r_number(e)));
        }
        lines.push(info);
        for c in m.matrix.columns() {
            let row: Vec<String> = order.iter().map(|&i| format_fixed(c[i], 6)).collect();
            lines.push(format!(" {}", row.join("  ")));
        }
        lines.push(String::new());
    }
    Ok(join_lines(&lines))
}

/// A MEME alphabet definition, written by [`write_alphabet`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AlphabetDefinition {
    /// The core letters, in order.
    pub core: String,
    /// The complement of each core letter, for alphabets with one.
    pub complements: Option<String>,
    /// Ambiguous letters and the core letters each stands for.
    pub ambiguity: Vec<(String, String)>,
    /// `"DNA"`, `"RNA"` or `"PROTEIN"`, for alphabets that extend one.
    pub like: Option<String>,
    pub name: Option<String>,
    /// A name for each letter. If any are given, all core and complement
    /// letters need one.
    pub letter_names: Vec<(char, String)>,
    /// A colour for each letter, as hexadecimal RGB. If any are given, all
    /// letters need one.
    pub colours: Vec<(char, String)>,
}

/// Write a MEME custom alphabet definition.
pub fn write_alphabet(def: &AlphabetDefinition) -> Result<String> {
    let mut head = String::from("ALPHABET");
    if let Some(name) = &def.name {
        head.push_str(&format!(" \"{name}\""));
    }
    if let Some(like) = &def.like {
        if !matches!(like.as_str(), "DNA" | "RNA" | "PROTEIN") {
            return Err(Error::Argument(
                "like must be one of DNA, RNA and PROTEIN".into(),
            ));
        }
        head.push_str(&format!(" {like}-LIKE"));
    }
    let core: Vec<char> = def.core.chars().collect();
    let compl: Option<Vec<char>> = def.complements.as_ref().map(|c| c.chars().collect());
    if compl.as_ref().is_some_and(|c| c.len() != core.len()) {
        return Err(Error::Argument(
            "core and complements do not have the same number of letters".into(),
        ));
    }
    let describe = |l: char| -> Result<String> {
        let mut s = l.to_string();
        if !def.letter_names.is_empty() {
            match def.letter_names.iter().find(|(k, _)| *k == l) {
                Some((_, n)) => s.push_str(&format!(" \"{n}\"")),
                None => {
                    return Err(Error::Argument(format!(
                        "the letter '{l}' has no name in letter_names"
                    )));
                }
            }
        }
        if !def.colours.is_empty() {
            match def.colours.iter().find(|(k, _)| *k == l) {
                Some((_, c)) => s.push_str(&format!(" {c}")),
                None => {
                    return Err(Error::Argument(format!("the letter '{l}' has no colour")));
                }
            }
        }
        Ok(s)
    };
    let mut lines = vec![head, String::new(), "# Core symbols".into()];
    for (i, &l) in core.iter().enumerate() {
        // Each complementary pair is declared once, on the line of its first
        // letter.
        if let Some(c) = &compl
            && core[..i].contains(&c[i])
            && c[core.iter().position(|&x| x == c[i]).unwrap()] == l
        {
            continue;
        }
        let mut line = describe(l)?;
        if let Some(c) = &compl {
            line.push_str(&format!(" ~ {}", describe(c[i])?));
        }
        lines.push(line);
    }
    if !def.ambiguity.is_empty() {
        lines.push(String::new());
        lines.push("# Ambiguous symbols".into());
        lines.extend(def.ambiguity.iter().map(|(k, v)| format!("{k} = {v}")));
    }
    Ok(lines.join("\n") + "\n")
}
