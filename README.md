# universalmotif (Rust)

This crate is a Rust port of [universalmotif](https://bioconductor.org/packages/universalmotif/),
an R/Bioconductor package for working with sequence motifs. It covers most
of what the R package does: creating motifs and converting between count,
probability, weight, information content and contribution matrices,
reading and writing the common motif formats, exact P-values, scanning
sequences, comparing, merging and clustering motifs, enrichment testing,
and a fair amount of sequence handling such as k-lets, complexity,
shuffling and masking. The supported formats are MEME, JASPAR, HOMER,
TRANSFAC, CIS-BP, UniPROBE and plain matrices. It is also the engine
behind the Python version of the package,
[universalmotif-py](https://github.com/bjmt/universalmotif-py).

Results are checked against the R package itself. Where the two differ, it
is on purpose, usually because a bug in the R package was fixed here, and
[DIVERGENCES.md](DIVERGENCES.md) explains each case.

## A quick example

Let us make a motif from a few aligned sites, look at it as a log-odds
matrix, and use it to scan a sequence:

```rust
use universalmotif::MotifType;
use universalmotif::convert::convert_type;
use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::scan::{ScanOptions, scan_table};

fn main() -> Result<(), universalmotif::Error> {
    let opts = CreateOptions { nsites: Some(20.0), pseudocount: 1.0, ..Default::default() };
    let motif = from_strings(&["TGACTCA", "TGAGTCA", "TGACTAA"], None, &opts)?;
    println!("{}", convert_type(&motif, MotifType::Pwm, None, false)?);

    let seqs = ["CCCCTGACTCACCCCC"];
    let hits = scan_table(&[motif], &seqs, &ScanOptions { pvalue: 1e-3, ..Default::default() })?;
    for i in 0..hits.len() {
        println!("{}..{} {}", hits.start[i], hits.end[i], hits.matched(i, &seqs));
    }
    Ok(())
}
```

Both strands are scanned by default, so a site that is close to
palindromic, as this one is, can show up on both.

The [documentation](https://docs.rs/universalmotif) has a map of the
modules and more examples.

## A few things beyond the R package

- Scanning, comparison, merging, clustering and enrichment work for amino
  acid and custom alphabets, not just DNA and RNA.
- Gapped motifs are scanned on both strands, and a little over twice as
  fast as scanning every combination of gap sizes on its own.
- Very large sequence sets can be scanned straight from a FASTA file,
  gzipped or not, with `scan::scan_fasta`, which reads and scans a batch at
  a time, so a human genome takes a few hundred MB. Alternatively,
  sequences can be packed once at two bits per base with
  `io::fasta::read_encoded` and scanned as often as needed with
  `scan::scan_encoded`.
- Genomes can be read a piece at a time with `genome::Genome`, from FASTA
  files (plain or compressed with bgzip, using or building the same
  indexes as `samtools faidx`) and from 2bit files, so windows and peaks
  can be taken from a genome without loading it.

## Testing

```sh
cargo test
```

The tests compare the crate against reference values exported from the R
package into `tests/reference`. To regenerate them from a universalmotif
checkout:

```sh
Rscript scripts/make_reference.R ~/universalmotif
```

This needs R with the `pkgload` and `jsonlite` packages, plus the
dependencies of universalmotif.

## Benchmarks

`examples/bench.rs` times the main functions at realistic sizes, one task
per run, so that peak memory can be measured around it:

```sh
cargo build --release --example bench
UM_BENCH_MOTIFS=motifs.meme /usr/bin/time -l target/release/examples/bench scan_many 4
```

Without `UM_BENCH_MOTIFS` it uses 400 random motifs. The tasks are listed
in the file.

## Licence

MIT
