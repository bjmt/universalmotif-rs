//! Tests at scale: long sequences and many of them, wide motifs, big
//! alphabets and many motifs, where chunking, batching, block boundaries,
//! integer widths and the merging of parallel work could go wrong. Results
//! must not depend on the number of threads either.

use std::io::Write;
use std::path::PathBuf;

use flate2::Compression;
use flate2::write::{DeflateEncoder, GzEncoder};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use universalmotif::cluster::{Linkage, TreeOptions, hclust, motif_tree};
use universalmotif::compare::{CompareOptions, compare_matrix, compare_motifs};
use universalmotif::complexity::{ComplexityMethod, ComplexityOptions, sequence_complexity};
use universalmotif::coocc::{CooccHit, CooccOptions, motif_coocc};
use universalmotif::create::{CreateOptions, random};
use universalmotif::enrich::{EnrichOptions, EnrichTest, enrich_motifs};
use universalmotif::finder::{FinderOptions, motif_finder};
use universalmotif::genome::Genome;
use universalmotif::io::fasta;
use universalmotif::kmer::{BkgOptions, WindowOptions, count_kmers, get_bkg};
use universalmotif::merge::{MergeSimilarOptions, merge_similar};
use universalmotif::ops::{add_gap, ungap};
use universalmotif::packed::EncodedSequences;
use universalmotif::peaks::{PeakHit, PeaksOptions, motif_peaks};
use universalmotif::pvalue::{PvalueOptions, min_pvalues, motif_pvalue, pvalue_to_score};
use universalmotif::ranges::GRange;
use universalmotif::scan::{
    HitTable, MatchValue, NoOverlaps, OverlapPriority, ScanOptions, motif_matches, scan_encoded,
    scan_fasta, scan_table,
};
use universalmotif::scoring::score_match;
use universalmotif::shuffle::{ShuffleMethod, ShuffleOptions, create_sequences, shuffle_sequences};
use universalmotif::{Alphabet, Matrix, Motif, MotifType};

/// The scanner's unit of work, in window starts.
const CHUNK: usize = 1 << 16;

fn rng(seed: u64) -> StdRng {
    StdRng::seed_from_u64(seed)
}

/// Random DNA. A messy sequence also has lower case stretches, runs of N and
/// the odd letter outside the alphabet.
fn dna(rng: &mut StdRng, len: usize, messy: bool) -> String {
    let mut out: Vec<u8> = Vec::with_capacity(len + 300);
    while out.len() < len {
        let r = if messy {
            rng.random_range(0..1000)
        } else {
            999
        };
        if r == 0 {
            let n = rng.random_range(1..300);
            out.extend(std::iter::repeat_n(b'N', n));
        } else if r == 1 {
            out.push(b"RYnX-"[rng.random_range(0..5)]);
        } else if r < 10 {
            for _ in 0..rng.random_range(1..60) {
                out.push(b"acgt"[rng.random_range(0..4)]);
            }
        } else {
            out.push(b"ACGT"[rng.random_range(0..4)]);
        }
    }
    out.truncate(len);
    String::from_utf8(out).unwrap()
}

/// Random letters of an alphabet, with the odd letter outside it.
fn letters_of(rng: &mut StdRng, letters: &[char], len: usize, odd: Option<char>) -> String {
    (0..len)
        .map(|_| match odd {
            Some(c) if rng.random_range(0..500) == 0 => c,
            _ => letters[rng.random_range(0..letters.len())],
        })
        .collect()
}

fn random_motif(width: usize, alphabet: Alphabet, seed: u64) -> Motif {
    let n = alphabet.len();
    let opts = CreateOptions {
        pseudocount: 1.0,
        nsites: Some(100.0),
        bkg: Some(vec![1.0 / n as f64; n]),
        ..Default::default()
    };
    random(width, alphabet, &opts, seed).unwrap()
}

fn consensus(m: &Motif) -> String {
    m.matrix
        .columns()
        .map(|c| {
            let best = (0..c.len()).max_by(|&a, &b| c[a].total_cmp(&c[b])).unwrap();
            m.alphabet.letters()[best]
        })
        .collect()
}

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            o => o as char,
        })
        .collect()
}

/// A temporary directory of its own, removed when dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> TempDir {
        let dir =
            std::env::temp_dir().join(format!("universalmotif-scale-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        TempDir(dir)
    }

    fn join(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Fail on the first row where two hit tables differ, rather than printing
/// both whole.
fn assert_same_table(a: &HitTable, b: &HitTable, what: &str) {
    assert_eq!(a.len(), b.len(), "{what}: different numbers of hits");
    for i in 0..a.len() {
        let row = |t: &HitTable| {
            (
                t.motif[i],
                t.sequence[i],
                t.start[i],
                t.end[i],
                t.minus[i],
                t.score[i].to_bits(),
                t.score_pct[i].to_bits(),
                t.pvalue[i].to_bits(),
            )
        };
        assert_eq!(row(a), row(b), "{what}: row {i}");
    }
    assert!(a == b, "{what}: tables differ outside the columns compared");
}

fn fasta_text(records: &[(String, String)], width: usize, crlf: bool) -> String {
    let nl = if crlf { "\r\n" } else { "\n" };
    let mut out = String::new();
    for (name, seq) in records {
        out.push('>');
        out.push_str(name);
        out.push_str(nl);
        if width == 0 {
            out.push_str(seq);
            out.push_str(nl);
        } else {
            for chunk in seq.as_bytes().chunks(width) {
                out.push_str(std::str::from_utf8(chunk).unwrap());
                out.push_str(nl);
            }
        }
    }
    out
}

/// Compress text in BGZF blocks of at most `block` bytes, as bgzip does,
/// returning the file and the (compressed, uncompressed) offset of each block.
fn bgzip(text: &[u8], block: usize) -> (Vec<u8>, Vec<(u64, u64)>) {
    let mut out = Vec::new();
    let mut blocks = Vec::new();
    let put = |data: &[u8], out: &mut Vec<u8>| {
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::fast());
        enc.write_all(data).unwrap();
        let deflated = enc.finish().unwrap();
        let bsize = u16::try_from(deflated.len() + 25).unwrap();
        out.extend_from_slice(&[
            0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
        ]);
        out.extend_from_slice(&bsize.to_le_bytes());
        out.extend_from_slice(&deflated);
        let mut crc = flate2::Crc::new();
        crc.update(data);
        out.extend_from_slice(&crc.sum().to_le_bytes());
        out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    };
    let mut uoff = 0;
    for chunk in text.chunks(block) {
        blocks.push((out.len() as u64, uoff));
        put(chunk, &mut out);
        uoff += chunk.len() as u64;
    }
    put(&[], &mut out);
    (out, blocks)
}

fn gzi(blocks: &[(u64, u64)]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&((blocks.len() - 1) as u64).to_le_bytes());
    for &(c, u) in &blocks[1..] {
        out.extend_from_slice(&c.to_le_bytes());
        out.extend_from_slice(&u.to_le_bytes());
    }
    out
}

/// A version 0 2bit file, little-endian, with each sequence's runs of N and
/// of lower case letters as its N and mask blocks.
fn twobit(chroms: &[(String, String)]) -> Vec<u8> {
    let runs = |seq: &str, f: &dyn Fn(u8) -> bool| {
        let mut out: Vec<(u32, u32)> = Vec::new();
        for (i, b) in seq.bytes().enumerate() {
            if f(b) {
                match out.last_mut() {
                    Some(r) if (r.0 + r.1) as usize == i => r.1 += 1,
                    _ => out.push((i as u32, 1)),
                }
            }
        }
        out
    };
    let mut records = Vec::new();
    for (_, seq) in chroms {
        let n = runs(seq, &|b| b == b'N' || b == b'n');
        let m = runs(seq, &|b| b.is_ascii_lowercase());
        let mut r = Vec::new();
        r.extend_from_slice(&(seq.len() as u32).to_le_bytes());
        for blocks in [&n, &m] {
            r.extend_from_slice(&(blocks.len() as u32).to_le_bytes());
            for b in blocks.iter() {
                r.extend_from_slice(&b.0.to_le_bytes());
            }
            for b in blocks.iter() {
                r.extend_from_slice(&b.1.to_le_bytes());
            }
        }
        r.extend_from_slice(&0u32.to_le_bytes());
        let mut packed = vec![0u8; seq.len().div_ceil(4)];
        for (i, b) in seq.bytes().enumerate() {
            let code = match b.to_ascii_uppercase() {
                b'C' => 1,
                b'A' => 2,
                b'G' => 3,
                _ => 0,
            };
            packed[i / 4] |= code << (6 - 2 * (i % 4));
        }
        r.extend_from_slice(&packed);
        records.push(r);
    }
    let mut out = Vec::new();
    out.extend_from_slice(&0x1A412743u32.to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes());
    out.extend_from_slice(&(chroms.len() as u32).to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes());
    let index_len: usize = chroms.iter().map(|c| 1 + c.0.len() + 4).sum();
    let mut offset = 16 + index_len;
    for ((name, _), r) in chroms.iter().zip(&records) {
        out.push(name.len() as u8);
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(&(offset as u32).to_le_bytes());
        offset += r.len();
    }
    for r in records {
        out.extend_from_slice(&r);
    }
    out
}

// ---------------------------------------------------------------------------
// Packed sequences

/// What decoding gives back: the alphabet's letters in upper case and
/// anything else as it was.
fn decoded_dna(s: &str) -> Vec<u8> {
    s.bytes()
        .map(|b| match b {
            b'a' | b'c' | b'g' | b't' => b.to_ascii_uppercase(),
            _ => b,
        })
        .collect()
}

#[test]
fn packed_pieces_cross_word_boundaries_at_every_offset() {
    let mut r = rng(11);
    let seq = dna(&mut r, 3_001, true);
    let mut seq = seq.into_bytes();
    // Unknown letters at every offset within a byte of the packed form.
    for (k, at) in (40..60).enumerate() {
        seq[at * 7 + k % 4] = b'N';
    }
    seq[1_000..1_013].fill(b'n');
    let seq = String::from_utf8(seq).unwrap();
    let whole = EncodedSequences::new(&[&seq], &Alphabet::Dna, 1).unwrap();
    assert_eq!(whole.decode(0).as_bytes(), decoded_dna(&seq));

    // The same sequence fed in pieces of every size from 1 to 9, after a first
    // piece that puts the rest at every offset within a byte.
    for first in 0..8 {
        for size in 1..10 {
            let mut set = EncodedSequences::empty(&Alphabet::Dna).unwrap();
            let mut rec = set.record();
            rec.push(&seq.as_bytes()[..first]);
            for piece in seq.as_bytes()[first..].chunks(size) {
                rec.push(piece);
            }
            assert_eq!(rec.len(), seq.len());
            set.push(rec);
            assert!(set == whole, "first piece {first}, then pieces of {size}");
        }
    }
    // Random piece sizes, some empty.
    let mut set = EncodedSequences::empty(&Alphabet::Dna).unwrap();
    let mut rec = set.record();
    let mut at = 0;
    while at < seq.len() {
        let n = r.random_range(0..40).min(seq.len() - at);
        rec.push(&seq.as_bytes()[at..at + n]);
        at += n;
    }
    set.push(rec);
    assert!(set == whole);

    // Every range starting and ending near byte boundaries decodes to the
    // same letters.
    let want = decoded_dna(&seq);
    let mut out = Vec::new();
    for a in (0..64).chain(996..1_020).chain(2_990..3_001) {
        for b in a..(a + 70).min(seq.len() + 1) {
            out.clear();
            out.extend_from_slice(b"xy");
            whole.decode_into(0, a, b, &mut out);
            assert_eq!(&out[2..], &want[a..b], "{a}..{b}");
            assert_eq!(&out[..2], b"xy", "decoding appends");
        }
    }

    // Other alphabets take a byte per letter; check them the same way, with
    // amino acids matched in either case and custom letters exactly.
    let aa: Vec<char> = Alphabet::Aa.letters().to_vec();
    let protein: String = letters_of(&mut r, &aa, 2_000, Some('X'))
        .chars()
        .enumerate()
        .map(|(i, c)| {
            if i % 7 == 0 {
                c.to_ascii_lowercase()
            } else {
                c
            }
        })
        .collect();
    let custom = Alphabet::custom(('!'..='~').filter(|c| *c != '>')).unwrap();
    let cletters = custom.letters().to_vec();
    let text = letters_of(&mut r, &cletters, 2_000, Some('\u{7f}'));
    for (alphabet, s) in [(Alphabet::Aa, &protein), (custom, &text)] {
        let whole = EncodedSequences::new(&[s], &alphabet, 2).unwrap();
        let want: String = if alphabet == Alphabet::Aa {
            s.chars()
                .map(|c| {
                    if aa.contains(&c.to_ascii_uppercase()) {
                        c.to_ascii_uppercase()
                    } else {
                        c
                    }
                })
                .collect()
        } else {
            s.clone()
        };
        assert_eq!(whole.decode(0), want, "{alphabet}");
        for size in [1, 2, 3, 5, 64, 1_999] {
            let mut set = EncodedSequences::empty(&alphabet).unwrap();
            let mut rec = set.record();
            s.as_bytes().chunks(size).for_each(|p| rec.push(p));
            set.push(rec);
            assert!(set == whole, "{alphabet}, pieces of {size}");
        }
    }
}

#[test]
fn long_packed_sequence_round_trips() {
    // A few megabases, longer than any internal buffer, with all sorts of
    // letters.
    let mut r = rng(12);
    let seqs = [
        dna(&mut r, 3_000_017, true),
        String::new(),
        "N".repeat(9),
        dna(&mut r, 5, false),
    ];
    let one = EncodedSequences::new(&seqs, &Alphabet::Dna, 1).unwrap();
    let four = EncodedSequences::new(&seqs, &Alphabet::Dna, 4).unwrap();
    assert!(one == four);
    assert_eq!(one.total_len(), seqs.iter().map(String::len).sum::<usize>());
    for (i, s) in seqs.iter().enumerate() {
        assert_eq!(one.seq_len(i), s.len());
        assert!(one.decode(i).as_bytes() == decoded_dna(s), "sequence {i}");
    }
    // About a quarter of the text.
    assert!(one.heap_bytes() < seqs[0].len() / 3, "{}", one.heap_bytes());
}

#[test]
fn fasta_batches_match_whole_reads() {
    let dir = TempDir::new("batches");
    let mut r = rng(13);
    let records: Vec<(String, String)> = (0..60)
        .map(|i| {
            let len = match i % 6 {
                0 => 0,
                1 => r.random_range(1..10),
                2 => r.random_range(1_000..50_000),
                _ => r.random_range(10..2_000),
            };
            (format!("seq{i} desc"), dna(&mut r, len, true))
        })
        .collect();
    let path = dir.join("b.fa");
    std::fs::write(&path, fasta_text(&records, 61, true)).unwrap();
    let (names, whole) = fasta::read_encoded(&path, &Alphabet::Dna, 0).unwrap();
    assert_eq!(names.len(), records.len());
    for (i, (n, s)) in records.iter().enumerate() {
        assert_eq!(&names[i], n);
        assert!(whole.decode(i).as_bytes() == decoded_dna(s), "record {i}");
    }
    let total: usize = records.iter().map(|r| r.1.len()).sum();
    for batch in [1, 999, 50_000, total / 3, total, usize::MAX] {
        let mut got: Vec<(String, String)> = Vec::new();
        let mut sizes = Vec::new();
        let f = std::io::BufReader::new(std::fs::File::open(&path).unwrap());
        fasta::read_encoded_batches(f, &Alphabet::Dna, 0, batch, |n, b| {
            sizes.push(b.total_len());
            assert_eq!(n.len(), b.len());
            for (k, name) in n.into_iter().enumerate() {
                got.push((name, b.decode(k)));
            }
            Ok(())
        })
        .unwrap();
        assert_eq!(got.len(), records.len(), "batch size {batch}");
        for (i, (n, s)) in got.iter().enumerate() {
            assert_eq!(n, &names[i]);
            assert_eq!(*s, whole.decode(i), "batch size {batch}, record {i}");
        }
        // Every batch but the last holds at least the batch size.
        let last = sizes.len() - 1;
        assert!(
            sizes[..last].iter().all(|&s| s >= batch),
            "{batch}: {sizes:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// Scanning

fn dna_motifs(widths: &[usize], seed: u64) -> Vec<Motif> {
    widths
        .iter()
        .enumerate()
        .map(|(i, &w)| random_motif(w, Alphabet::Dna, seed + i as u64))
        .collect()
}

/// Check every hit against `score_match` and `motif_pvalue`.
fn check_hits(motifs: &[Motif], seqs: &[&str], t: &HitTable, cutoff: f64) {
    for (mi, m) in motifs.iter().enumerate() {
        let rows: Vec<usize> = (0..t.len()).filter(|&i| t.motif[i] == mi as u32).collect();
        let words: Vec<String> = rows.iter().map(|&i| t.matched(i, seqs)).collect();
        let scores = score_match(m, &words, false).unwrap();
        let queries: Vec<(usize, f64)> = scores.iter().map(|&s| (0, s)).collect();
        let p = motif_pvalue(std::slice::from_ref(m), &queries, &PvalueOptions::default()).unwrap();
        for (k, &i) in rows.iter().enumerate() {
            assert_eq!(t.score[i], scores[k], "motif {mi}, hit {i}: {}", words[k]);
            assert_eq!(t.pvalue[i], p[k], "motif {mi}, hit {i}");
            assert!(t.pvalue[i] <= cutoff);
            assert_eq!((t.end[i] - t.start[i]) as usize, m.width());
        }
    }
}

/// Whether the window of `m` at `start` on each strand should be a hit,
/// worked out from `score_match` and `motif_pvalue` alone.
fn naive_windows(m: &Motif, seq: &str, starts: &[usize], cutoff: f64) -> Vec<(usize, bool)> {
    let w = m.width();
    let mut words = Vec::new();
    let mut keys = Vec::new();
    for &s in starts {
        let word = &seq[s..s + w];
        if !word.bytes().all(|b| b"ACGTacgt".contains(&b)) {
            continue;
        }
        words.push(word.to_ascii_uppercase());
        keys.push((s, false));
        words.push(revcomp(word));
        keys.push((s, true));
    }
    let scores = score_match(m, &words, false).unwrap();
    let queries: Vec<(usize, f64)> = scores.iter().map(|&s| (0, s)).collect();
    let p = motif_pvalue(std::slice::from_ref(m), &queries, &PvalueOptions::default()).unwrap();
    keys.into_iter()
        .zip(p)
        .filter(|(_, p)| *p <= cutoff)
        .map(|(k, _)| k)
        .collect()
}

#[test]
fn scans_do_not_depend_on_threads_or_chunks() {
    let mut r = rng(21);
    let motifs = dna_motifs(&[5, 8, 12, 17, 24, 31], 100);
    // One long clean sequence with sites planted across the scanner's chunk
    // boundaries, one long messy one, and many short ones, some shorter than
    // the motifs.
    let mut long = dna(&mut r, 7 * CHUNK + 1_234, false).into_bytes();
    let mut planted = Vec::new();
    for (k, m) in motifs.iter().enumerate() {
        let b = k + 1;
        let at = b * CHUNK - 1 - (k * 5 + b) % m.width();
        let minus = b % 2 == 0;
        let site = if minus {
            revcomp(&consensus(m))
        } else {
            consensus(m)
        };
        long[at..at + m.width()].copy_from_slice(site.as_bytes());
        planted.push((k, at, minus));
    }
    let long = String::from_utf8(long).unwrap();
    let mut seqs = vec![long, dna(&mut r, 250_000, true)];
    for _ in 0..2_000 {
        let len = r.random_range(0..250);
        let messy = r.random_range(0..4) == 0;
        seqs.push(dna(&mut r, len, messy));
    }
    let refs: Vec<&str> = seqs.iter().map(String::as_str).collect();
    let cutoff = 1e-4;
    let opts = |nthreads| ScanOptions {
        pvalue: cutoff,
        nthreads,
        ..Default::default()
    };
    let one = scan_table(&motifs, &refs, &opts(1)).unwrap();
    assert!(one.len() > 200, "{} hits", one.len());
    for n in [2, 4] {
        assert_same_table(
            &one,
            &scan_table(&motifs, &refs, &opts(n)).unwrap(),
            "threads",
        );
    }
    let enc = EncodedSequences::new(&refs, &Alphabet::Dna, 3).unwrap();
    assert_same_table(
        &one,
        &scan_encoded(&motifs, &enc, &opts(4)).unwrap(),
        "encoded",
    );
    for i in (0..one.len()).step_by(7) {
        assert_eq!(one.matched_encoded(i, &enc), one.matched(i, &refs));
    }
    check_hits(&motifs, &refs, &one, cutoff);

    // Every window starting near a chunk boundary of the long sequence is
    // reported exactly when it should be, planted sites included.
    for (mi, m) in motifs.iter().enumerate() {
        let w = m.width();
        let starts: Vec<usize> = (1..=7)
            .flat_map(|b| b * CHUNK - w - 1..=b * CHUNK + 1)
            .collect();
        let want = naive_windows(m, refs[0], &starts, cutoff);
        let got: Vec<(usize, bool)> = (0..one.len())
            .filter(|&i| one.motif[i] == mi as u32 && one.sequence[i] == 0)
            .map(|i| (one.start[i] as usize, one.minus[i]))
            .filter(|(s, _)| starts.contains(s))
            .collect();
        let mut want_sorted = want.clone();
        want_sorted.sort();
        assert_eq!(got, want_sorted, "motif {mi}");
        // The planted sites score the best possible, so they are hits
        // whenever the motif can reach the cutoff at all.
        let reachable =
            min_pvalues(std::slice::from_ref(m), &PvalueOptions::default()).unwrap()[0] <= cutoff;
        for &(k, at, minus) in &planted {
            if k == mi {
                assert_eq!(got.contains(&(at, minus)), reachable, "motif {mi} at {at}");
            }
        }
    }

    // Overlap removal at scale, by P-value and by score, is the same on
    // any number of threads.
    for by in [OverlapPriority::Pvalue, OverlapPriority::Score] {
        let no = ScanOptions {
            no_overlaps: Some(NoOverlaps {
                by,
                ignore_strand: true,
                ignore_motif: false,
            }),
            ..opts(1)
        };
        let a = scan_table(&motifs, &refs, &ScanOptions { nthreads: 4, ..no }).unwrap();
        assert!(a.len() < one.len());
        // Overlap removal is the same on one thread, and leaves no two hits
        // of a motif overlapping on a sequence.
        if by == OverlapPriority::Pvalue {
            assert_same_table(&a, &scan_table(&motifs, &refs, &no).unwrap(), "no overlaps");
        }
        for i in 1..a.len() {
            if a.motif[i] == a.motif[i - 1] && a.sequence[i] == a.sequence[i - 1] {
                assert!(a.start[i] >= a.end[i - 1], "{by:?}: overlap at row {i}");
            }
        }
    }

    // Match matrices summarize the same hits on any number of threads.
    for value in [MatchValue::Matches, MatchValue::Counts, MatchValue::Scores] {
        let m1 = motif_matches(&motifs, &enc, &opts(1), value).unwrap();
        if value == MatchValue::Counts {
            assert!(m1 == motif_matches(&motifs, &enc, &opts(4), value).unwrap());
        }
        let mut want: Vec<(u32, u32, f64)> = Vec::new();
        let mut order: Vec<usize> = (0..one.len()).collect();
        order.sort_by_key(|&i| (one.sequence[i], one.motif[i]));
        for i in order {
            let key = (one.sequence[i] as u32, one.motif[i] as u32);
            match want.last_mut() {
                Some(l) if (l.0, l.1) == key => match value {
                    MatchValue::Matches => {}
                    MatchValue::Counts => l.2 += 1.0,
                    MatchValue::Scores => l.2 = l.2.max(one.score[i]),
                },
                _ => want.push((
                    key.0,
                    key.1,
                    if value == MatchValue::Scores {
                        one.score[i]
                    } else {
                        1.0
                    },
                )),
            }
        }
        let got: Vec<(u32, u32, f64)> = (0..m1.value.len())
            .map(|k| (m1.sequence[k], m1.motif[k], m1.value[k]))
            .collect();
        assert_eq!(got, want, "{value}");
    }
}

#[test]
fn match_matrices_span_blocks() {
    // More than the 4M positions motif_matches scans at a time, as many
    // sequences, so blocks end between sequences, and one of a few
    // megabases split into many chunks.
    let mut r = rng(22);
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 300, 15_000, None, 5, 4).unwrap();
    seqs.push(dna(&mut r, 100_000, true));
    seqs.push(dna(&mut r, 2_500_000, false));
    let motifs = dna_motifs(&[10], 300);
    let enc = EncodedSequences::new(&seqs, &Alphabet::Dna, 4).unwrap();
    assert!(enc.total_len() > 1 << 22);
    let opts = ScanOptions {
        pvalue: 1e-3,
        nthreads: 4,
        ..Default::default()
    };
    let t = scan_encoded(&motifs, &enc, &opts).unwrap();
    let t1 = scan_encoded(
        &motifs,
        &enc,
        &ScanOptions {
            nthreads: 1,
            ..opts.clone()
        },
    )
    .unwrap();
    assert_same_table(&t, &t1, "threads");
    let m = motif_matches(&motifs, &enc, &opts, MatchValue::Counts).unwrap();
    let mut counts = vec![0.0; seqs.len()];
    for &s in &t.sequence {
        counts[s as usize] += 1.0;
    }
    let want: Vec<(u32, f64)> = counts
        .iter()
        .enumerate()
        .filter(|(_, c)| **c > 0.0)
        .map(|(i, &c)| (i as u32, c))
        .collect();
    let got: Vec<(u32, f64)> = m.sequence.iter().copied().zip(m.value.clone()).collect();
    assert_eq!(got, want);
}

#[test]
fn fasta_scans_match_table_scans() {
    let dir = TempDir::new("scanfasta");
    let mut r = rng(31);
    let motifs = dna_motifs(&[6, 9, 14, 21], 400);
    let mut records: Vec<(String, String)> = Vec::new();
    for i in 0..40 {
        let len = match i % 8 {
            0 => 0,
            1 => r.random_range(1..8),
            2 => r.random_range(100_000..300_000),
            _ => r.random_range(10..5_000),
        };
        records.push((format!("rec{i} some words"), dna(&mut r, len, true)));
    }
    let seqs: Vec<&str> = records.iter().map(|r| r.1.as_str()).collect();
    let opts = ScanOptions {
        pvalue: 1e-4,
        nthreads: 1,
        ..Default::default()
    };
    let want = scan_table(&motifs, &seqs, &opts).unwrap();
    assert!(want.len() > 100);

    let plain = fasta_text(&records, 61, true);
    let mut gz = GzEncoder::new(Vec::new(), Compression::fast());
    gz.write_all(fasta_text(&records, 0, false).as_bytes())
        .unwrap();
    let files = [
        ("plain.fa", plain.clone().into_bytes()),
        ("one-line.fa.gz", gz.finish().unwrap()),
        ("blocks.fa.gz", bgzip(plain.as_bytes(), 65_280).0),
    ];
    for (name, data) in files {
        let path = dir.join(name);
        std::fs::write(&path, data).unwrap();
        for nthreads in [1, 3] {
            let got = scan_fasta(
                &motifs,
                &path,
                &ScanOptions {
                    nthreads,
                    ..opts.clone()
                },
                false,
            )
            .unwrap();
            let names: Vec<&str> = records.iter().map(|r| r.0.as_str()).collect();
            assert_eq!(got.names, names, "{name}");
            assert_eq!(
                got.lengths,
                seqs.iter().map(|s| s.len()).collect::<Vec<_>>()
            );
            assert_same_table(&want, &got.table, name);
            for i in 0..want.len() {
                assert_eq!(got.matched(i), want.matched(i, &seqs), "{name} hit {i}");
            }
        }
    }
}

/// Every window starting at `starts` that a gapped motif should hit, with
/// its strand and gap sizes, worked out from `score_match` and
/// `motif_pvalue` on the motif without its gaps.
fn naive_gapped(
    m: &Motif,
    seq: &str,
    starts: &[usize],
    cutoff: f64,
) -> Vec<(i64, bool, Vec<usize>)> {
    let plain = ungap(m, true);
    let g = &m.gap;
    let mut bounds = vec![0];
    bounds.extend(&g.gaploc);
    bounds.push(m.width());
    let widths: Vec<usize> = bounds.windows(2).map(|b| b[1] - b[0]).collect();
    let mut combos: Vec<Vec<usize>> = vec![Vec::new()];
    for (&lo, &hi) in g.mingap.iter().zip(&g.maxgap) {
        combos = combos
            .iter()
            .flat_map(|c| {
                (lo..=hi).map(move |s| {
                    let mut c = c.clone();
                    c.push(s);
                    c
                })
            })
            .collect();
    }
    let mut words = Vec::new();
    let mut keys = Vec::new();
    for &s in starts {
        for gaps in &combos {
            let w = m.width() + gaps.iter().sum::<usize>();
            if s + w > seq.len() {
                continue;
            }
            let window = &seq[s..s + w];
            if !window.bytes().all(|b| b"ACGTacgt".contains(&b)) {
                continue;
            }
            for minus in [false, true] {
                let oriented = if minus {
                    revcomp(window)
                } else {
                    window.to_ascii_uppercase()
                };
                let mut word = String::new();
                let mut at = 0;
                for (k, &pw) in widths.iter().enumerate() {
                    word.push_str(&oriented[at..at + pw]);
                    at += pw + gaps.get(k).copied().unwrap_or(0);
                }
                words.push(word);
                keys.push((s as i64, minus, gaps.clone()));
            }
        }
    }
    let scores = score_match(&plain, &words, false).unwrap();
    let q: Vec<(usize, f64)> = scores.iter().map(|&s| (0, s)).collect();
    let p = motif_pvalue(&[plain], &q, &PvalueOptions::default()).unwrap();
    let mut out: Vec<(i64, bool, Vec<usize>)> = keys
        .into_iter()
        .zip(p)
        .filter(|(_, p)| *p <= cutoff)
        .map(|(k, _)| k)
        .collect();
    out.sort();
    out
}

#[test]
fn gapped_scans_across_chunks_and_at_the_combination_limit() {
    let mut r = rng(33);
    let motifs = [
        add_gap(
            &random_motif(10, Alphabet::Dna, 330),
            Some(&[4]),
            &[0],
            &[40],
        )
        .unwrap(),
        add_gap(
            &random_motif(12, Alphabet::Dna, 331),
            Some(&[3, 8]),
            &[0, 1],
            &[12, 10],
        )
        .unwrap(),
        // Two gaps of 1 to 256 letters make 65,536 combinations, as many as
        // a hit can record.
        add_gap(
            &random_motif(16, Alphabet::Dna, 332),
            Some(&[5, 10]),
            &[1, 1],
            &[256, 256],
        )
        .unwrap(),
    ];
    let mut long = dna(&mut r, 3 * CHUNK + 700, false).into_bytes();
    // Best sites with their widest gaps, across each chunk boundary.
    let mut planted = Vec::new();
    for (k, m) in motifs.iter().enumerate() {
        let b = k + 1;
        let site = consensus(m);
        let g = &m.gap;
        let mut bounds = vec![0];
        bounds.extend(&g.gaploc);
        bounds.push(m.width());
        let mut text = String::new();
        for (j, piece) in bounds.windows(2).enumerate() {
            text.push_str(&site[piece[0]..piece[1]]);
            if let Some(&gap) = g.maxgap.get(j) {
                text.push_str(&"A".repeat(gap));
            }
        }
        let at = b * CHUNK - text.len() / 2;
        let minus = k == 1;
        let text = if minus { revcomp(&text) } else { text };
        long[at..at + text.len()].copy_from_slice(text.as_bytes());
        planted.push((k, at as i64, minus, g.maxgap.clone()));
    }
    let long = String::from_utf8(long).unwrap();
    let records: Vec<(String, String)> = vec![
        ("long".into(), long),
        ("messy".into(), dna(&mut r, 50_000, true)),
        ("short".into(), dna(&mut r, 30, false)),
    ];
    let seqs: Vec<&str> = records.iter().map(|r| r.1.as_str()).collect();
    let dir = TempDir::new("gapped");
    let path = dir.join("g.fa");
    std::fs::write(&path, fasta_text(&records, 60, false)).unwrap();

    // The motif with 65,536 combinations tries that many windows at each
    // start, so it gets a cutoff to match.
    let groups = [(0..2, 1e-5), (2..3, 1e-9)];
    for (range, cutoff) in groups {
        let group = &motifs[range.clone()];
        let opts = ScanOptions {
            pvalue: cutoff,
            nthreads: 4,
            correct_gaps: false,
            ..Default::default()
        };
        let t = scan_table(group, &seqs, &opts).unwrap();
        let t1 = scan_table(
            group,
            &seqs,
            &ScanOptions {
                nthreads: 1,
                ..opts.clone()
            },
        )
        .unwrap();
        assert_same_table(&t, &t1, "gapped, threads");
        let f = scan_fasta(group, &path, &opts, true).unwrap();
        assert_same_table(&t, &f.table, "gapped, scan_fasta");
        for i in 0..t.len() {
            assert_eq!(f.matched(i), t.matched(i, &seqs));
        }

        // Every hit scores what its letters outside the gaps score without
        // the gaps.
        for (mi, m) in group.iter().enumerate() {
            let plain = ungap(m, true);
            let rows: Vec<usize> = (0..t.len()).filter(|&i| t.motif[i] == mi as u32).collect();
            let words: Vec<String> = rows
                .iter()
                .map(|&i| t.matched(i, &seqs).replace('.', ""))
                .collect();
            let scores = score_match(&plain, &words, false).unwrap();
            for (k, &i) in rows.iter().enumerate() {
                assert_eq!(t.score[i], scores[k], "motif {mi}, hit {i}");
                let gaps = t.gap_sizes(i).unwrap();
                let w = m.width() + gaps.iter().sum::<usize>();
                assert_eq!((t.end[i] - t.start[i]) as usize, w);
                assert!(t.pvalue[i] <= cutoff);
            }
        }
        for (k, at, minus, gaps) in &planted {
            if !range.contains(k) {
                continue;
            }
            let mi = (k - range.start) as u32;
            assert!(
                (0..t.len()).any(|i| t.motif[i] == mi
                    && t.start[i] == *at
                    && t.minus[i] == *minus
                    && t.gap_sizes(i).unwrap() == *gaps),
                "planted site of motif {k}"
            );
        }
        if range.start > 0 {
            continue;
        }

        // Windows near the chunk boundaries, every combination tried, for
        // the motifs with fewer combinations.
        for (mi, m) in group.iter().enumerate() {
            let max_w = m.width() + m.gap.maxgap.iter().sum::<usize>();
            let starts: Vec<usize> = (1..=3)
                .flat_map(|b| b * CHUNK - max_w - 1..=b * CHUNK + 1)
                .collect();
            let want = naive_gapped(m, seqs[0], &starts, cutoff);
            let mut got: Vec<(i64, bool, Vec<usize>)> = (0..t.len())
                .filter(|&i| t.motif[i] == mi as u32 && t.sequence[i] == 0)
                .filter(|&i| starts.contains(&(t.start[i] as usize)))
                .map(|i| (t.start[i], t.minus[i], t.gap_sizes(i).unwrap()))
                .collect();
            got.sort();
            assert_eq!(got, want, "motif {mi}");
        }
    }
    // One more combination is too many.
    let over = add_gap(&motifs[2], Some(&[5, 10]), &[1, 0], &[256, 256]).unwrap();
    assert!(scan_table(&[over], &seqs, &ScanOptions::default()).is_err());
}

/// A two-letter PWM whose scores span about `2 * h` thousand grid points,
/// with a rare letter `b` that any hit must contain, so that its P-value
/// table above a cutoff of 0.003 is long however few windows hit.
fn long_tail_motif(h: f64, low: f64) -> Motif {
    let m = Matrix::from_columns(&[vec![-h, h], vec![0.0, low]]).unwrap();
    Motif::builder(m, Alphabet::custom(['a', 'b']).unwrap())
        .motif_type(Some(MotifType::Pwm))
        .bkg(Some(vec![0.999, 0.001]))
        .build()
        .unwrap()
}

#[test]
fn long_score_tails_are_batched_and_rebuilt() {
    // Forty motifs with tails of 140,000 to 220,000 grid points each: more than
    // the 4M points scan_encoded holds at once, so its motifs come in two
    // batches, and more than scan_fasta keeps, so it rebuilds the rest.
    let motifs: Vec<Motif> = (0..40)
        .map(|i| long_tail_motif(70.0 + i as f64, 0.001 * (1 + i % 3) as f64))
        .collect();
    let mut r = rng(32);
    let records: Vec<(String, String)> = (0..12)
        .map(|i| {
            let s: String = (0..r.random_range(500..3_000))
                .map(|_| if r.random_range(0..50) == 0 { 'b' } else { 'a' })
                .collect();
            (format!("s{i}"), s)
        })
        .collect();
    let seqs: Vec<&str> = records.iter().map(|r| r.1.as_str()).collect();
    let opts = ScanOptions {
        pvalue: 0.003,
        nthreads: 2,
        ..Default::default()
    };
    let t = scan_table(&motifs, &seqs, &opts).unwrap();
    // A window hits exactly when it holds a b.
    let windows = seqs
        .iter()
        .map(|s| {
            s.as_bytes()
                .windows(2)
                .filter(|w| w.contains(&b'b'))
                .count()
        })
        .sum::<usize>();
    assert_eq!(t.len(), windows * motifs.len());
    for (mi, m) in motifs.iter().enumerate() {
        let rows: Vec<usize> = (0..t.len()).filter(|&i| t.motif[i] == mi as u32).collect();
        let queries: Vec<(usize, f64)> = rows.iter().map(|&i| (0, t.score[i])).collect();
        let p = motif_pvalue(std::slice::from_ref(m), &queries, &PvalueOptions::default()).unwrap();
        for (k, &i) in rows.iter().enumerate() {
            assert_eq!(t.pvalue[i], p[k], "motif {mi}");
            let word = t.matched(i, &seqs);
            assert_eq!(t.score[i], score_match(m, &[&word], false).unwrap()[0]);
        }
    }
    let single = scan_table(
        &motifs,
        &seqs,
        &ScanOptions {
            nthreads: 1,
            ..opts.clone()
        },
    )
    .unwrap();
    assert_same_table(&t, &single, "threads");

    let dir = TempDir::new("tails");
    let path = dir.join("ab.fa");
    std::fs::write(&path, fasta_text(&records, 70, false)).unwrap();
    let f = scan_fasta(&motifs, &path, &opts, true).unwrap();
    assert_same_table(&t, &f.table, "scan_fasta");
}

// ---------------------------------------------------------------------------
// Genomes on disk

/// How the letters of a FASTA file written by `fasta_text` sit in it, to find
/// which letter a file offset falls on.
struct Layout {
    /// Offset of each record's first letter.
    starts: Vec<u64>,
    width: u64,
    nl: u64,
}

impl Layout {
    fn new(chroms: &[(String, String)], width: usize, crlf: bool) -> Layout {
        let nl = if crlf { 2 } else { 1 };
        let mut starts = Vec::new();
        let mut at = 0u64;
        for (name, seq) in chroms {
            at += name.len() as u64 + 1 + nl;
            starts.push(at);
            at += seq.len() as u64 + seq.len().div_ceil(width) as u64 * nl;
        }
        Layout {
            starts,
            width: width as u64,
            nl,
        }
    }

    /// The record and position at file offset `x`, or the next letter
    /// after it.
    fn position(&self, x: u64) -> Option<(usize, u64)> {
        let r = self.starts.partition_point(|&s| s <= x).checked_sub(1)?;
        let d = x - self.starts[r];
        let line = d / (self.width + self.nl);
        Some((
            r,
            line * self.width + (d % (self.width + self.nl)).min(self.width),
        ))
    }
}

fn genome_chroms(r: &mut StdRng) -> Vec<(String, String)> {
    let mut chroms = vec![
        ("chr1".to_string(), dna(r, 1_200_003, false)),
        ("chr2".to_string(), dna(r, 400_001, false)),
        ("chrM".to_string(), dna(r, 16_569, false)),
    ];
    for k in 0..5 {
        chroms.push((format!("chrUn{k}"), dna(r, k * 2 + 1, false)));
    }
    chroms.push(("chrN".to_string(), "N".repeat(70_000)));
    // Soft-masked stretches and runs of N, some long, on the big ones.
    for (_, s) in chroms.iter_mut().take(2) {
        let mut b = std::mem::take(s).into_bytes();
        for _ in 0..400 {
            let at = r.random_range(0..b.len() - 5_000);
            let len = r.random_range(1..5_000);
            if r.random_range(0..3) == 0 {
                b[at..at + len].fill(b'N');
            } else {
                b[at..at + len].make_ascii_lowercase();
            }
        }
        *s = String::from_utf8(b).unwrap();
    }
    chroms
}

fn check_genome(g: &Genome, chroms: &[(String, String)], extra: &[(usize, u64)], seed: u64) {
    let mut r = rng(seed);
    assert_eq!(g.len(), chroms.len());
    for (i, (name, seq)) in chroms.iter().enumerate() {
        assert_eq!(&g.names()[i], name);
        assert_eq!(g.lengths()[i], seq.len() as u64);
    }
    // Whole chromosomes.
    for (name, seq) in chroms {
        assert!(
            g.fetch(name, 0, seq.len() as u64).unwrap() == *seq,
            "{name}"
        );
    }
    let mut ranges = Vec::new();
    let mut want = Vec::new();
    let mut add = |i: usize, a: u64, b: u64| {
        ranges.push(GRange::new(i, a as i64, b as i64));
        want.push(chroms[i].1[a as usize..b as usize].to_string());
    };
    // Random pieces, mostly short and some long.
    for _ in 0..1_500 {
        let i = r.random_range(0..chroms.len());
        let len = chroms[i].1.len() as u64;
        let a = r.random_range(0..len);
        let max = if r.random_range(0..50) == 0 {
            200_000
        } else {
            400
        };
        let b = r.random_range(a..=len.min(a + max));
        add(i, a, b);
    }
    // Pieces around each position given, such as line and block ends.
    for &(i, p) in extra {
        let len = chroms[i].1.len() as u64;
        for (da, db) in [(0, 1), (1, 0), (1, 1), (2, 3), (70, 90), (0, 0)] {
            let (a, b) = (p.saturating_sub(da), (p + db).min(len));
            if a <= b {
                add(i, a, b);
            }
        }
    }
    let one = g.fetch_many(&ranges, 1).unwrap();
    let four = g.fetch_many(&ranges, 4).unwrap();
    for (k, w) in want.iter().enumerate() {
        assert!(one[k] == *w, "range {k}: {:?}", ranges[k]);
    }
    assert!(one == four);
    for (k, rg) in ranges.iter().enumerate().step_by(53) {
        let got = g
            .fetch_at(rg.chrom, rg.start as u64, rg.end as u64)
            .unwrap();
        assert!(got == want[k], "fetch_at {rg:?}");
    }
}

#[test]
fn genome_fetches_at_scale() {
    let dir = TempDir::new("genome");
    let mut r = rng(41);
    let chroms = genome_chroms(&mut r);
    for (k, (width, crlf)) in [(60, false), (80, true), (1_000_000, false)]
        .into_iter()
        .enumerate()
    {
        let text = fasta_text(&chroms, width, crlf);
        let layout = Layout::new(&chroms, width, crlf);
        // Line ends of the big chromosomes, a sample of them.
        let mut extra: Vec<(usize, u64)> = (0..2)
            .flat_map(|i| {
                (0..chroms[i].1.len() as u64)
                    .step_by(width * 97 + 1)
                    .map(move |p| (i, p - p % width as u64))
            })
            .collect();
        extra.extend((0..chroms.len()).map(|i| (i, chroms[i].1.len() as u64)));

        let path = dir.join(&format!("g{k}.fa"));
        std::fs::write(&path, &text).unwrap();
        check_genome(&Genome::open(&path).unwrap(), &chroms, &extra, 1);

        // bgzip, with blocks of the usual size and an odd one, each with and
        // without a .gzi index.
        let block = [65_280, 9_973, 0][k];
        if block > 0 {
            let (data, blocks) = bgzip(text.as_bytes(), block);
            let gz = dir.join(&format!("g{k}-{block}.fa.gz"));
            std::fs::write(&gz, &data).unwrap();
            let mut ends = extra.clone();
            ends.extend(blocks.iter().filter_map(|&(_, u)| layout.position(u)));
            check_genome(&Genome::open(&gz).unwrap(), &chroms, &ends, 2);
            std::fs::write(dir.join(&format!("g{k}-{block}.fa.gz.gzi")), gzi(&blocks)).unwrap();
            check_genome(&Genome::open(&gz).unwrap(), &chroms, &ends, 3);
        }
    }
    let tb = dir.join("g.2bit");
    let twobit_chroms: Vec<(String, String)> = chroms
        .iter()
        .map(|(n, s)| (n.clone(), s.replace('n', "N")))
        .collect();
    std::fs::write(&tb, twobit(&twobit_chroms)).unwrap();
    let extra: Vec<(usize, u64)> = (0..twobit_chroms.len())
        .map(|i| (i, twobit_chroms[i].1.len() as u64))
        .chain((0..4_000).map(|p| (0, p * 449)))
        .collect();
    check_genome(&Genome::open(&tb).unwrap(), &twobit_chroms, &extra, 4);
}

// ---------------------------------------------------------------------------
// Wide motifs and big alphabets

/// The P-value of a motif's best score, worked out directly: the chance of
/// drawing, at every position, a letter with that position's best entry on
/// the score grid.
fn best_score_pvalue(m: &Motif) -> f64 {
    let pwm = universalmotif::convert::convert_type(m, MotifType::Pwm, None, false).unwrap();
    pwm.matrix
        .columns()
        .map(|c| {
            let grid: Vec<i64> = c.iter().map(|&v| (v * 1000.0) as i64).collect();
            let best = *grid.iter().max().unwrap();
            grid.iter()
                .zip(&m.bkg)
                .filter(|(g, _)| **g == best)
                .map(|(_, b)| b)
                .sum::<f64>()
        })
        .product()
}

/// Check thresholds and P-values of a motif against each other, returning the
/// thresholds of the P-value ladder. Each call works out the score
/// distribution once, which for wide motifs is the slow part.
fn check_pvalue_consistency(m: &Motif, ladder: &[f64], what: &str) -> Vec<f64> {
    let ms = std::slice::from_ref(m);
    let o = PvalueOptions::default();
    let q: Vec<(usize, f64)> = ladder.iter().map(|&p| (0, p)).collect();
    let scores = pvalue_to_score(ms, &q, &o).unwrap();
    let min_p = min_pvalues(ms, &o).unwrap()[0];
    let (lo, hi) = universalmotif::scoring::motif_range(m, false).unwrap();
    // One batch of queries: each finite threshold and the grid point below
    // it, then scores across the whole range, then one above the best.
    let finite: Vec<usize> = (0..ladder.len())
        .filter(|&k| scores[k].is_finite())
        .collect();
    let mut queries: Vec<(usize, f64)> = finite
        .iter()
        .flat_map(|&k| [(0, scores[k]), (0, scores[k] - 0.001)])
        .collect();
    let n_grid = 400;
    queries.extend((0..=n_grid).map(|k| (0, lo + (hi - lo) * k as f64 / n_grid as f64)));
    queries.push((0, hi + 0.001));
    let p = motif_pvalue(ms, &queries, &o).unwrap();

    // Thresholds rise as the P-value falls, and each is the lowest grid score
    // whose P-value is within the cutoff.
    for k in 1..ladder.len() {
        assert!(scores[k] >= scores[k - 1], "{what}: {k}");
    }
    for k in 0..ladder.len() {
        if scores[k].is_infinite() {
            assert!(ladder[k] < min_p, "{what}: {} has no threshold", ladder[k]);
        }
    }
    for (j, &k) in finite.iter().enumerate() {
        let (at, below) = (p[2 * j], p[2 * j + 1]);
        assert!(
            scores[k] >= lo && scores[k] <= hi + 1e-9,
            "{what}: {}",
            scores[k]
        );
        assert!(at <= ladder[k], "{what}: {at} > {}", ladder[k]);
        if scores[k] > lo {
            assert!(below > ladder[k], "{what}: {below} <= {}", ladder[k]);
        }
    }
    // P-values fall as the score rises, across the whole range.
    let grid = &p[2 * finite.len()..2 * finite.len() + n_grid + 1];
    assert!((grid[0] - 1.0).abs() < 1e-9, "{what}: {}", grid[0]);
    assert!(grid.windows(2).all(|w| w[1] <= w[0]), "{what}");
    assert!(grid.iter().all(|&v| (0.0..=1.0).contains(&v)), "{what}");
    assert_eq!(*p.last().unwrap(), 0.0, "{what}");
    // The chance of the best score, checked directly.
    let best = best_score_pvalue(m);
    assert!(
        (min_p - best).abs() <= 1e-9 * best,
        "{what}: {min_p} against {best}"
    );
    scores
}

fn big_custom(n: usize) -> Alphabet {
    let letters: Vec<char> = ('!'..='~').filter(|c| *c != '>' && *c != ';').collect();
    Alphabet::custom(letters.into_iter().take(n)).unwrap()
}

#[test]
fn wide_motifs_and_big_alphabets_keep_pvalues_consistent() {
    let cases = [
        (random_motif(60, Alphabet::Dna, 501), "DNA width 60"),
        (random_motif(45, Alphabet::Rna, 502), "RNA width 45"),
        (random_motif(30, Alphabet::Aa, 503), "AA width 30"),
        (random_motif(12, big_custom(64), 504), "64 letters"),
        (random_motif(10, big_custom(65), 505), "65 letters"),
        (random_motif(8, big_custom(92), 506), "92 letters"),
    ];
    let ladder: Vec<f64> = (1..=40).map(|k| 10f64.powf(-0.5 * k as f64)).collect();
    let single: Vec<f64> = cases
        .iter()
        .flat_map(|(m, what)| check_pvalue_consistency(m, &ladder, what))
        .collect();
    // All the motifs at once, on more threads, give the same thresholds.
    let motifs: Vec<Motif> = cases.iter().map(|c| c.0.clone()).collect();
    let q: Vec<(usize, f64)> = (0..motifs.len())
        .flat_map(|i| ladder.iter().map(move |&p| (i, p)))
        .collect();
    let o = PvalueOptions {
        nthreads: 4,
        ..Default::default()
    };
    assert_eq!(pvalue_to_score(&motifs, &q, &o).unwrap(), single);
}

#[test]
fn wide_motif_scans_agree_with_score_match() {
    let mut r = rng(51);
    let cases: Vec<(Vec<Motif>, f64, Option<char>)> = vec![
        (
            vec![
                random_motif(30, Alphabet::Aa, 511),
                random_motif(3, Alphabet::Aa, 512),
            ],
            1e-3,
            Some('X'),
        ),
        (vec![random_motif(12, big_custom(64), 513)], 1e-3, Some(' ')),
        (vec![random_motif(9, big_custom(92), 514)], 1e-3, Some(' ')),
        (
            vec![
                random_motif(60, Alphabet::Dna, 515),
                random_motif(33, Alphabet::Dna, 516),
            ],
            1e-2,
            Some('N'),
        ),
    ];
    for (motifs, cutoff, odd) in cases {
        let alphabet = motifs[0].alphabet.clone();
        let letters = alphabet.letters().to_vec();
        let mut seqs: Vec<String> = (0..4)
            .map(|k| letters_of(&mut r, &letters, [CHUNK + 300, 2_000, 61, 7][k], odd))
            .collect();
        // A planted best site across the chunk boundary, so there is a hit.
        let site = consensus(&motifs[0]);
        let at = CHUNK - site.len() / 2;
        seqs[0].replace_range(at..at + site.len(), &site);
        let refs: Vec<&str> = seqs.iter().map(String::as_str).collect();
        let opts = ScanOptions {
            pvalue: cutoff,
            nthreads: 3,
            ..Default::default()
        };
        let t = scan_table(&motifs, &refs, &opts).unwrap();
        let t1 = scan_table(
            &motifs,
            &refs,
            &ScanOptions {
                nthreads: 1,
                ..opts.clone()
            },
        )
        .unwrap();
        assert_same_table(&t, &t1, "threads");
        assert!(t.len() > 1, "{alphabet}");
        check_hits(&motifs, &refs, &t, cutoff);
        assert!(
            (0..t.len()).any(|i| t.start[i] == at as i64 && t.motif[i] == 0),
            "{alphabet}: planted site"
        );
        // Every window of the shorter sequences, checked one by one, with
        // one P-value query per motif since each works out its whole score
        // distribution.
        let nucleotide = alphabet.is_nucleotide();
        for (mi, m) in motifs.iter().enumerate() {
            let w = m.width();
            let mut words = Vec::new();
            let mut keys = Vec::new();
            for (si, s) in refs.iter().enumerate().skip(1) {
                let chars: Vec<char> = s.chars().collect();
                for i in 0..(chars.len() + 1).saturating_sub(w) {
                    let word: String = chars[i..i + w].iter().collect();
                    if word.chars().any(|c| !letters.contains(&c)) {
                        continue;
                    }
                    if nucleotide {
                        words.push(revcomp(&word));
                        keys.push((si as i64, i as i64, true));
                    }
                    words.push(word);
                    keys.push((si as i64, i as i64, false));
                }
            }
            let scores = score_match(m, &words, false).unwrap();
            let q: Vec<(usize, f64)> = scores.iter().map(|&s| (0, s)).collect();
            let p = motif_pvalue(std::slice::from_ref(m), &q, &PvalueOptions::default()).unwrap();
            let mut want: Vec<(i64, i64, bool)> = keys
                .into_iter()
                .zip(p)
                .filter(|(_, p)| *p <= cutoff)
                .map(|(k, _)| k)
                .collect();
            want.sort();
            let got: Vec<(i64, i64, bool)> = (0..t.len())
                .filter(|&i| t.motif[i] == mi as u32 && t.sequence[i] > 0)
                .map(|i| (i64::from(t.sequence[i]), t.start[i], t.minus[i]))
                .collect();
            assert_eq!(got, want, "{alphabet}, motif {mi}");
        }
    }
}

// ---------------------------------------------------------------------------
// Many motifs

#[test]
fn many_motifs_compare_and_cluster_consistently() {
    let motifs: Vec<Motif> = (0..200)
        .map(|i| random_motif(6 + i % 10, Alphabet::Dna, 600 + i as u64))
        .collect();
    let opts = |nthreads| CompareOptions {
        min_overlap: 6,
        nthreads,
        ..Default::default()
    };
    let m1 = compare_matrix(&motifs, &opts(1), false).unwrap();
    let m4 = compare_matrix(&motifs, &opts(4), false).unwrap();
    let bits =
        |v: &Vec<Vec<f64>>| -> Vec<u64> { v.iter().flatten().map(|x| x.to_bits()).collect() };
    assert_eq!(bits(&m1.scores), bits(&m4.scores));
    for i in 0..motifs.len() {
        assert!((m1.scores[i][i] - 1.0).abs() < 1e-9, "self comparison {i}");
    }

    // A tree of all 200 is the clustering of the comparison scores.
    let tree = motif_tree(
        &motifs,
        &TreeOptions {
            nthreads: 4,
            ..Default::default()
        },
    )
    .unwrap();
    let dist: Vec<Vec<f64>> = m1
        .scores
        .iter()
        .map(|r| r.iter().map(|s| (1.0 - s) / 2.0).collect())
        .collect();
    assert_eq!(tree, hclust(&dist, Linkage::Complete).unwrap());
    assert_eq!(tree.merge.len(), motifs.len() - 1);
    let mut order = tree.order.clone();
    order.sort_unstable();
    assert_eq!(order, (0..motifs.len()).collect::<Vec<_>>());
    assert!(tree.height.windows(2).all(|w| w[0] <= w[1]));
    for k in [1, 2, 17, 200] {
        let cut = tree.cutree(k).unwrap();
        let mut ids = cut.clone();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), k);
    }

    // A few queries with P-values give the same rows as the full matrix with
    // P-values, on any number of threads.
    let queries = [0, 57, 199];
    let rows = compare_motifs(&motifs, &queries, &opts(1)).unwrap();
    assert!(rows == compare_motifs(&motifs, &queries, &opts(4)).unwrap());
    let full = compare_matrix(&motifs[..120], &opts(4), true).unwrap();
    let part = compare_motifs(&motifs[..120], &[3, 4], &opts(2)).unwrap();
    for c in &part {
        assert_eq!(c.score.to_bits(), full.scores[c.query][c.target].to_bits());
        assert_eq!(
            c.pvalue.to_bits(),
            full.pvalues[c.query][c.target].to_bits()
        );
        assert_eq!(
            c.qvalue.to_bits(),
            full.qvalues[c.query][c.target].to_bits()
        );
    }
    assert!(rows.iter().all(|c| (0.0..=1.0).contains(&c.pvalue)));
}

#[test]
fn merging_many_near_duplicates() {
    // Twenty distinct motifs, each with four noisy copies.
    let mut r = rng(71);
    let bases: Vec<Motif> = (0..20)
        .map(|i| random_motif(8 + i % 5, Alphabet::Dna, 700 + i as u64))
        .collect();
    let mut motifs = Vec::new();
    for b in &bases {
        motifs.push(b.clone());
        for _ in 0..4 {
            let cols: Vec<Vec<f64>> = b
                .matrix
                .columns()
                .map(|c| {
                    let noisy: Vec<f64> =
                        c.iter().map(|&v| v + r.random_range(0.0..0.05)).collect();
                    let s: f64 = noisy.iter().sum();
                    noisy.iter().map(|v| v / s).collect()
                })
                .collect();
            let m = Motif::builder(Matrix::from_columns(&cols).unwrap(), Alphabet::Dna)
                .motif_type(Some(MotifType::Ppm))
                .nsites(Some(100.0))
                .build()
                .unwrap();
            motifs.push(m);
        }
    }
    let opts = |nthreads| MergeSimilarOptions {
        nthreads,
        ..Default::default()
    };
    let one = merge_similar(&motifs, &opts(1)).unwrap();
    let four = merge_similar(&motifs, &opts(4)).unwrap();
    assert!(one == four);
    assert_eq!(one.clusters.len(), motifs.len());
    let n = one.motifs.len();
    // Numbered in order of first appearance.
    let mut next = 0;
    for &c in &one.clusters {
        assert!(c <= next);
        if c == next {
            next += 1;
        }
    }
    assert_eq!(next, n);
    // Clusters are the connected groups of motifs whose comparison q-value is
    // within the cutoff in either direction.
    let q = compare_matrix(
        &motifs,
        &CompareOptions {
            min_overlap: 5,
            nthreads: 4,
            ..Default::default()
        },
        true,
    )
    .unwrap()
    .qvalues;
    let mut label = vec![usize::MAX; motifs.len()];
    let mut n_groups = 0;
    for first in 0..motifs.len() {
        if label[first] != usize::MAX {
            continue;
        }
        let mut stack = vec![first];
        label[first] = n_groups;
        while let Some(i) = stack.pop() {
            for j in 0..motifs.len() {
                if label[j] == usize::MAX && q[i][j].min(q[j][i]) <= 0.01 {
                    label[j] = n_groups;
                    stack.push(j);
                }
            }
        }
        n_groups += 1;
    }
    assert_eq!(one.clusters, label);
    // Copies of a motif end up together, and there are at most as many
    // clusters as distinct motifs.
    for b in 0..bases.len() {
        let c = one.clusters[5 * b];
        assert!(
            (1..5).all(|k| one.clusters[5 * b + k] == c),
            "copies of {b}"
        );
    }
    assert!(n <= bases.len());
    for (c, m) in one.motifs.iter().enumerate() {
        let members = one.clusters.iter().filter(|&&x| x == c).count();
        if members == 1 {
            let i = one.clusters.iter().position(|&x| x == c).unwrap();
            assert!(*m == motifs[i]);
        }
    }
}

// ---------------------------------------------------------------------------
// Other parallel work

#[test]
fn sequence_work_does_not_depend_on_threads() {
    let mut r = rng(81);
    let mut seqs: Vec<String> = (0..400)
        .map(|_| {
            let len = r.random_range(150..600);
            dna(&mut r, len, false)
        })
        .collect();
    for len in [5_000, 6_250, 9_999, 20_000, 51_234, 400_000] {
        seqs.push(dna(&mut r, len, false));
    }
    let refs: Vec<&str> = seqs.iter().map(String::as_str).collect();
    // The short ones and one of 20,000 letters.
    let some: Vec<&str> = refs[..400].iter().chain(&refs[403..404]).copied().collect();

    // Shuffles, by every method, with and without windows.
    for method in [
        ShuffleMethod::Euler,
        ShuffleMethod::Markov,
        ShuffleMethod::Linear,
    ] {
        for window in [
            None,
            Some(WindowOptions {
                size: 100.0,
                overlap: 20.0,
            }),
        ] {
            let o = |nthreads| ShuffleOptions {
                k: 3,
                method,
                seed: 9,
                window,
                nthreads,
            };
            let a = shuffle_sequences(&some, &o(1)).unwrap();
            assert_eq!(a, shuffle_sequences(&some, &o(4)).unwrap(), "{method:?}");
            // Euler keeps k-mer counts, and Linear the letters.
            let k = if method == ShuffleMethod::Euler { 3 } else { 1 };
            let counts = |s: &str| count_kmers(s, k, &['A', 'C', 'G', 'T']);
            if method != ShuffleMethod::Markov && window.is_none() {
                assert!(
                    a.iter().zip(&some).all(|(x, y)| counts(x) == counts(y)),
                    "{method:?}"
                );
            }
        }
    }
    let c1 = create_sequences(&['A', 'C', 'G', 'T'], 50, 3_000, None, 3, 1).unwrap();
    assert_eq!(
        c1,
        create_sequences(&['A', 'C', 'G', 'T'], 50, 3_000, None, 3, 4).unwrap()
    );

    // k-mer counts, whole and in windows, against counting the windows one by
    // one.
    let bkg = |nthreads, window, merge| BkgOptions {
        k: vec![1, 2, 4],
        window,
        merge,
        nthreads,
        ..Default::default()
    };
    let a = get_bkg(&refs, &[], &bkg(1, None, false)).unwrap();
    assert!(a == get_bkg(&refs, &[], &bkg(4, None, false)).unwrap());
    // Windows only fit the long sequences.
    let long = &refs[refs.len() - 6..];
    let window = Some(WindowOptions {
        size: 5_000.0,
        overlap: 1_250.0,
    });
    let a = get_bkg(long, &[], &bkg(1, window, false)).unwrap();
    assert!(a == get_bkg(long, &[], &bkg(4, window, false)).unwrap());
    let alphabet = ['A', 'C', 'G', 'T'];
    let four = universalmotif::kmer::kmers(&alphabet, 4);
    let mut cache: std::collections::HashMap<(usize, usize), Vec<u64>> = Default::default();
    let mut n_rows = 0;
    for row in a.iter().filter(|r| r.kmer.len() == 4) {
        let (si, s, e) = (row.sequence.unwrap(), row.start.unwrap(), row.end.unwrap());
        assert!(e - s <= 5_000 && e <= long[si].len());
        let counts = cache
            .entry((si, s))
            .or_insert_with(|| count_kmers(&long[si][s..e], 4, &alphabet));
        let idx = four.iter().position(|k| *k == row.kmer).unwrap();
        assert_eq!(row.count, counts[idx]);
        n_rows += 1;
    }
    // Windows step by 3,750 letters to the end of each sequence.
    let n_windows: usize = long.iter().map(|s| (s.len() - 1_250).div_ceil(3_750)).sum();
    assert_eq!(n_rows, 256 * n_windows);
    let merged = get_bkg(&refs, &[], &bkg(4, None, true)).unwrap();
    assert!(merged == get_bkg(&refs, &[], &bkg(1, None, true)).unwrap());
    let total: u64 = merged
        .iter()
        .filter(|r| r.kmer.len() == 1)
        .map(|r| r.count)
        .sum();
    assert_eq!(total as usize, refs.iter().map(|s| s.len()).sum::<usize>());

    // Complexity in windows, by every method.
    for method in [
        ComplexityMethod::WoottonFederhen,
        ComplexityMethod::WoottonFederhenFast,
        ComplexityMethod::Trifonov,
        ComplexityMethod::TrifonovFast,
        ComplexityMethod::Dust,
    ] {
        let o = |nthreads| ComplexityOptions {
            method,
            window: 200.0,
            nthreads,
            ..Default::default()
        };
        let a = sequence_complexity(&some[some.len() - 30..], &o(1)).unwrap();
        assert!(a == sequence_complexity(&some[some.len() - 30..], &o(4)).unwrap());
        assert!(a.iter().all(|w| w.complexity.is_finite()), "{method:?}");
    }
}

#[test]
fn enrichment_and_hit_statistics_do_not_depend_on_threads() {
    let mut r = rng(91);
    let motifs = dna_motifs(&[6, 8, 10, 12], 900);
    let mut seqs: Vec<String> = (0..600).map(|_| dna(&mut r, 300, false)).collect();
    // Plant the first motif in a third of the sequences, near the middle.
    let site = consensus(&motifs[0]);
    for s in seqs.iter_mut().step_by(3) {
        let at = 140 + r.random_range(0..10);
        s.replace_range(at..at + site.len(), &site);
    }
    let refs: Vec<&str> = seqs.iter().map(String::as_str).collect();
    for test in [EnrichTest::Seqs, EnrichTest::Sites] {
        let o = |nthreads| EnrichOptions {
            pvalue: 1e-3,
            test,
            nthreads,
            ..Default::default()
        };
        let a = enrich_motifs(&motifs, &refs, None, &o(1)).unwrap();
        assert!(a == enrich_motifs(&motifs, &refs, None, &o(4)).unwrap());
        assert!(a[0].pvalue < 1e-5, "{test:?}: {}", a[0].pvalue);
    }

    // Co-occurrence and peaks from the hits of many sequences.
    let t = scan_table(
        &motifs,
        &refs,
        &ScanOptions {
            pvalue: 1e-3,
            nthreads: 4,
            ..Default::default()
        },
    )
    .unwrap();
    let coocc: Vec<CooccHit> = (0..t.len())
        .map(|i| CooccHit {
            motif: t.motif[i] as usize,
            sequence: t.sequence[i] as usize,
            start: t.start[i],
        })
        .collect();
    let co = |nthreads| CooccOptions {
        max_distance: Some(50),
        nthreads,
        ..Default::default()
    };
    let a = motif_coocc(&coocc, motifs.len(), refs.len(), &co(1)).unwrap();
    assert!(a == motif_coocc(&coocc, motifs.len(), refs.len(), &co(4)).unwrap());
    let peaks: Vec<PeakHit> = (0..t.len())
        .map(|i| PeakHit {
            motif: t.motif[i] as usize,
            sequence: t.sequence[i] as usize,
            centre: (t.start[i] + t.end[i] - 1) as f64 / 2.0,
            score: t.score[i],
        })
        .collect();
    let po = |nthreads| PeaksOptions {
        nthreads,
        ..Default::default()
    };
    let a = motif_peaks(&peaks, 300, &po(1)).unwrap();
    assert!(a == motif_peaks(&peaks, 300, &po(4)).unwrap());
    assert_eq!(a.len(), motifs.len());
}

#[test]
fn proximity_does_not_depend_on_threads() {
    use universalmotif::proximity::{
        ProximityHit, ProximityMethod, ProximityOptions, motif_proximity,
    };
    // Chromosomes of a few megabases, anchors across them, and the hits of
    // a few motifs, one of them drawn to the anchors.
    let mut r = rng(93);
    let lengths = [5_000_000i64, 3_000_000, 1_000_000];
    let universe: Vec<GRange> = lengths
        .iter()
        .enumerate()
        .map(|(c, &l)| GRange::new(c, 0, l))
        .collect();
    let anchors: Vec<GRange> = (0..1_500)
        .map(|_| {
            let c = r.random_range(0..3);
            let s = r.random_range(0..lengths[c] - 500);
            GRange::new(c, s, s + r.random_range(1..500))
        })
        .collect();
    let mut hits = Vec::new();
    for m in 0..5 {
        for _ in 0..300 {
            let (c, s) = if m == 0 && r.random_range(0..2) == 0 {
                let a = anchors[r.random_range(0..anchors.len())];
                (
                    a.chrom,
                    (a.start + r.random_range(-200..200)).clamp(0, lengths[a.chrom] - 20),
                )
            } else {
                let c = r.random_range(0..3);
                (c, r.random_range(0..lengths[c] - 20))
            };
            hits.push(ProximityHit {
                motif: m * 3,
                range: GRange::new(c, s, s + 12),
                score: r.random_range(5.0..15.0),
            });
        }
    }
    let chrom_lengths: Vec<Option<i64>> = lengths.iter().map(|&l| Some(l)).collect();
    for method in [
        ProximityMethod::Binomial,
        ProximityMethod::Permutation,
        ProximityMethod::Ks,
    ] {
        let o = |nthreads| ProximityOptions {
            method,
            n_perm: 50,
            nthreads,
            ..Default::default()
        };
        let a = motif_proximity(&hits, &anchors, &universe, &chrom_lengths, &o(1)).unwrap();
        let b = motif_proximity(&hits, &anchors, &universe, &chrom_lengths, &o(4)).unwrap();
        assert!(a == b, "{method:?}");
        assert_eq!(
            a.iter().map(|p| p.motif).collect::<Vec<_>>(),
            [0, 3, 6, 9, 12]
        );
        let best = a
            .iter()
            .min_by(|x, y| x.pvalue.total_cmp(&y.pvalue))
            .unwrap();
        assert_eq!(best.motif, 0, "{method:?}");
    }
}

#[test]
fn motif_finder_does_not_depend_on_threads() {
    let mut r = rng(95);
    let site = "TGACGTCATG";
    let seqs: Vec<String> = (0..150)
        .map(|i| {
            let mut s = dna(&mut r, 200, false);
            if i % 2 == 0 {
                let at = r.random_range(0..190);
                let w = if r.random_range(0..2) == 0 {
                    site.to_string()
                } else {
                    revcomp(site)
                };
                s.replace_range(at..at + 10, &w);
            }
            s
        })
        .collect();
    let refs: Vec<&str> = seqs.iter().map(String::as_str).collect();
    let o = |nthreads| FinderOptions {
        min_width: 6,
        max_width: 10,
        nmotifs: 3,
        nthreads,
        ..Default::default()
    };
    let a = motif_finder(&refs, None, &o(1)).unwrap();
    assert!(a == motif_finder(&refs, None, &o(4)).unwrap());
    assert!(!a.is_empty());
}

#[test]
fn matched_backgrounds_at_scale() {
    use universalmotif::background::{MatchOptions, match_bkg};
    let mut r = rng(97);
    let gc_seq = |r: &mut StdRng, len: usize, gc: f64| -> String {
        (0..len)
            .map(|_| {
                let x: f64 = r.random();
                let y = r.random_range(0..2);
                if x < gc {
                    b"GC"[y] as char
                } else {
                    b"AT"[y] as char
                }
            })
            .collect()
    };
    let targets: Vec<String> = (0..1_000)
        .map(|_| {
            let (len, gc) = (r.random_range(100..400), r.random_range(0.3..0.7));
            gc_seq(&mut r, len, gc)
        })
        .collect();
    let universe: Vec<String> = (0..20_000)
        .map(|_| {
            let (len, gc) = (r.random_range(80..450), r.random_range(0.2..0.8));
            gc_seq(&mut r, len, gc)
        })
        .collect();
    let opts = MatchOptions {
        n_per_target: 3,
        ..Default::default()
    };
    let a = match_bkg(&targets, &universe, &[], &opts).unwrap();
    assert!(a == match_bkg(&targets, &universe, &[], &opts).unwrap());
    assert_eq!(a.universe_index.len(), 3_000);
    assert_eq!(a.fallback_targets, 0);
    // Drawn without replacement across all targets.
    let mut seen = a.universe_index.clone();
    seen.sort_unstable();
    seen.dedup();
    assert_eq!(seen.len(), a.universe_index.len());
    // Draws are close in GC content to their targets.
    for (t, chunk) in a.universe_index.chunks(3).enumerate() {
        for &u in chunk {
            assert!((a.gc_targets[t] - a.gc_universe[u]).abs() < 0.16, "{t} {u}");
        }
    }
}
