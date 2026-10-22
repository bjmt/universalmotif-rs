# universalmotif (Rust)

This is the Rust version of universalmotif, a library for working with sequence motifs. It provides a
motif type with validation, conversion between count, probability, weight,
information content and contribution matrices, and operations such as
reverse complementing and trimming. It is a port of the core of the
[universalmotif](https://bioconductor.org/packages/universalmotif/) R
package, and it is the engine behind the Python version
([universalmotif-py](../universalmotif-py)).

The library is at an early stage. Motif creation, type conversion, the
basic operations, reading and writing motif files (MEME, JASPAR, HOMER,
TRANSFAC, CIS-BP, UniPROBE, plain matrices and FASTA), exact P-values,
sequence scanning, motif comparison, merging and clustering, k-let
backgrounds, sequence complexity, shuffling, random sequences, implanting
sites, *de novo* motif discovery, enrichment testing, co-occurrence,
positional enrichment, proximity to genomic anchors and
composition-matched backgrounds are done. Scanning, comparison, merging,
clustering and enrichment work for amino acid and custom alphabets as well
as DNA and RNA.

## Example

```rust
use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::convert::convert_type;
use universalmotif::MotifType;

let opts = CreateOptions::default();
let motif = from_strings(&["TATAWAW"], None, &opts)?;
let pwm = convert_type(&motif, MotifType::Pwm, Some(1.0), false)?;
println!("{pwm}");
# Ok::<(), universalmotif::Error>(())
```

## Testing

```sh
cargo test
```

The tests compare the crate against reference values exported from the R
package. To regenerate them from a universalmotif checkout:

```sh
Rscript scripts/make_reference.R ~/universalmotif
```

This needs R with the `pkgload` and `jsonlite` packages, plus the
dependencies of universalmotif. The places where this crate deliberately
differs from the R package are listed in [DIVERGENCES.md](DIVERGENCES.md).

## Benchmarks

`examples/bench.rs` times the main functions at realistic sizes, one task
per run, so peak memory can be measured around it:

```sh
cargo build --release --example bench
UM_BENCH_MOTIFS=motifs.meme /usr/bin/time -l target/release/examples/bench scan_many 4
```

Without `UM_BENCH_MOTIFS` it uses 400 random motifs. The tasks are listed in
the file.

## License

MIT
