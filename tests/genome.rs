//! Reading genomes a piece at a time from FASTA, bgzip and 2bit files.

use std::io::Write;
use std::path::Path;

use flate2::Compression;
use flate2::write::DeflateEncoder;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use universalmotif::genome::Genome;

fn random_chroms(rng: &mut StdRng) -> Vec<(String, String)> {
    [("chr1", 5000), ("chr2", 777), ("chrM", 61), ("tiny", 3)]
        .iter()
        .map(|&(name, len)| {
            let seq: String = (0..len)
                .map(|_| {
                    let c = b"ACGTacgtNN"[rng.random_range(0..10)];
                    c as char
                })
                .collect();
            (name.to_string(), seq)
        })
        .collect()
}

fn fasta_text(chroms: &[(String, String)], width: usize, crlf: bool) -> String {
    let nl = if crlf { "\r\n" } else { "\n" };
    let mut out = String::new();
    for (name, seq) in chroms {
        out.push_str(&format!(">{name} some description{nl}"));
        for chunk in seq.as_bytes().chunks(width) {
            out.push_str(std::str::from_utf8(chunk).unwrap());
            out.push_str(nl);
        }
    }
    out
}

fn check_fetches(genome: &Genome, chroms: &[(String, String)], rng: &mut StdRng) {
    assert_eq!(
        genome.names(),
        chroms.iter().map(|c| c.0.clone()).collect::<Vec<_>>()
    );
    // Many pieces at once, given out of order, come back in the order given.
    let mut ranges = Vec::new();
    let mut want = Vec::new();
    for _ in 0..200 {
        let i = rng.random_range(0..chroms.len());
        let seq = &chroms[i].1;
        let a = rng.random_range(0..seq.len());
        let b = rng.random_range(a..=seq.len().min(a + 300));
        ranges.push((genome.position(&chroms[i].0).unwrap(), a as u64, b as u64));
        want.push(seq[a..b].to_string());
    }
    assert_eq!(genome.fetch_many(&ranges, 2).unwrap(), want);
    for (name, seq) in chroms {
        let len = seq.len() as u64;
        assert_eq!(genome.lengths()[genome.position(name).unwrap()], len);
        assert_eq!(genome.fetch(name, 0, len).unwrap(), *seq);
        for _ in 0..50 {
            let a = rng.random_range(0..len);
            let b = rng.random_range(a..=len);
            assert_eq!(
                genome.fetch(name, a, b).unwrap(),
                seq[a as usize..b as usize],
                "{name} {a}..{b}"
            );
        }
    }
}

#[test]
fn plain_fasta() {
    let dir = tempdir();
    let mut rng = StdRng::seed_from_u64(1);
    let chroms = random_chroms(&mut rng);
    for (width, crlf) in [(60, false), (61, true), (10_000, false)] {
        let path = dir.join(format!("g{width}.fa"));
        std::fs::write(&path, fasta_text(&chroms, width, crlf)).unwrap();
        let built = Genome::open(&path).unwrap();
        check_fetches(&built, &chroms, &mut rng);
        // A .fai next to the file is used as it is.
        let mut fai = String::new();
        let nl = if crlf { 2 } else { 1 };
        let mut offset = 0;
        for (name, seq) in &chroms {
            offset += format!(">{name} some description").len() + nl;
            let bases = width.min(seq.len());
            fai.push_str(&format!(
                "{name}\t{}\t{offset}\t{bases}\t{}\n",
                seq.len(),
                bases + nl
            ));
            offset += seq.len() + seq.len().div_ceil(width) * nl;
        }
        std::fs::write(dir.join(format!("g{width}.fa.fai")), fai).unwrap();
        check_fetches(&Genome::open(&path).unwrap(), &chroms, &mut rng);
    }
}

#[test]
fn errors() {
    let dir = tempdir();
    let ragged = dir.join("ragged.fa");
    std::fs::write(&ragged, ">a\nACGT\nAC\nACGT\n").unwrap();
    let e = Genome::open(&ragged).unwrap_err().to_string();
    assert!(e.contains("differ in length"), "{e}");
    let dup = dir.join("dup.fa");
    std::fs::write(&dup, ">a\nACGT\n>a x\nAC\n").unwrap();
    assert!(
        Genome::open(&dup)
            .unwrap_err()
            .to_string()
            .contains("more than once")
    );
    let ok = dir.join("ok.fa");
    std::fs::write(&ok, ">a\nACGT\n").unwrap();
    let g = Genome::open(&ok).unwrap();
    assert!(g.fetch("a", 2, 5).is_err());
    assert!(g.fetch("b", 0, 1).is_err());
    assert_eq!(g.fetch("a", 2, 2).unwrap(), "");
    // Ordinary gzip cannot be read in pieces.
    let gz = dir.join("plain.fa.gz");
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(b">a\nACGT\n").unwrap();
    std::fs::write(&gz, enc.finish().unwrap()).unwrap();
    let e = Genome::open(&gz).unwrap_err().to_string();
    assert!(e.contains("bgzip"), "{e}");
}

/// Compress text in BGZF blocks of at most `block` bytes, as bgzip does,
/// returning the file and its blocks as (compressed, uncompressed) offsets.
fn bgzip(text: &[u8], block: usize) -> (Vec<u8>, Vec<(u64, u64)>) {
    let mut out = Vec::new();
    let mut blocks = Vec::new();
    let put = |data: &[u8], out: &mut Vec<u8>| {
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
        enc.write_all(data).unwrap();
        let deflated = enc.finish().unwrap();
        let bsize = (deflated.len() + 25) as u16;
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

#[test]
fn bgzip_fasta() {
    let dir = tempdir();
    let mut rng = StdRng::seed_from_u64(2);
    let chroms = random_chroms(&mut rng);
    let text = fasta_text(&chroms, 60, false);
    let (data, blocks) = bgzip(text.as_bytes(), 1000);
    let path = dir.join("g.fa.gz");
    std::fs::write(&path, &data).unwrap();
    // Without .gzi or .fai, both are built.
    check_fetches(&Genome::open(&path).unwrap(), &chroms, &mut rng);
    // A .gzi lists every block but the first.
    let mut gzi = Vec::new();
    gzi.extend_from_slice(&((blocks.len() - 1) as u64).to_le_bytes());
    for &(c, u) in &blocks[1..] {
        gzi.extend_from_slice(&c.to_le_bytes());
        gzi.extend_from_slice(&u.to_le_bytes());
    }
    std::fs::write(dir.join("g.fa.gz.gzi"), gzi).unwrap();
    check_fetches(&Genome::open(&path).unwrap(), &chroms, &mut rng);
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
                b'T' => 0,
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

#[test]
fn twobit_files() {
    let dir = tempdir();
    let mut rng = StdRng::seed_from_u64(3);
    // 2bit keeps N, and lower case only as soft-masking.
    let chroms: Vec<(String, String)> = random_chroms(&mut rng)
        .into_iter()
        .map(|(n, s)| (n, s.replace('n', "N")))
        .collect();
    let path = dir.join("g.2bit");
    std::fs::write(&path, twobit(&chroms)).unwrap();
    check_fetches(&Genome::open(&path).unwrap(), &chroms, &mut rng);
}

fn tempdir() -> std::path::PathBuf {
    static COUNT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let dir =
        std::env::temp_dir().join(format!("universalmotif-genome-{}-{n}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    assert!(Path::new(&dir).exists());
    dir
}
