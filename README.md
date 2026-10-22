# umotif

umotif is a Rust library for working with sequence motifs. It provides a
motif type with validation, conversion between count, probability, weight,
information content and contribution matrices, and operations such as
reverse complementing and trimming. It is a port of the core of the
[universalmotif](https://bioconductor.org/packages/universalmotif/) R
package, and it is the engine behind the pyumotif Python package.

The library is at an early stage. Motif creation, type conversion, the
basic operations, and reading and writing MEME, JASPAR, HOMER, TRANSFAC,
CIS-BP, UniPROBE and plain matrix files are done; P-values, scanning,
comparison, shuffling and motif discovery will follow.

## Example

```rust
use umotif::create::{CreateOptions, from_strings};
use umotif::convert::convert_type;
use umotif::MotifType;

let opts = CreateOptions::default();
let motif = from_strings(&["TATAWAW"], None, &opts)?;
let pwm = convert_type(&motif, MotifType::Pwm, Some(1.0), false)?;
println!("{pwm}");
# Ok::<(), umotif::Error>(())
```

## Testing

```sh
cargo test
```

The tests compare umotif against reference values exported from
universalmotif. To regenerate them from a universalmotif checkout:

```sh
Rscript scripts/make_reference.R ~/universalmotif
```

This needs R with the `pkgload` and `jsonlite` packages, plus the
dependencies of universalmotif. The places where umotif deliberately
differs from universalmotif are listed in [DIVERGENCES.md](DIVERGENCES.md).

## License

MIT
