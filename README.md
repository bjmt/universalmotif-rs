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
TRANSFAC, CIS-BP, UniPROBE, plain matrices and FASTA), exact P-values and
sequence scanning are done; comparison, shuffling and motif discovery will
follow.

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

## License

MIT
