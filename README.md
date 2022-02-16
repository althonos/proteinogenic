# `proteinogenic` [![Star me](https://img.shields.io/github/stars/althonos/proteinogenic.svg?style=social&label=Star&maxAge=3600)](https://github.com/althonos/proteinogenic/stargazers)

*Chemical structure generation for protein sequences as [SMILES] string.*

[SMILES]: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system

[![Actions](https://img.shields.io/github/workflow/status/althonos/proteinogenic/Test?style=flat-square&maxAge=600)](https://github.com/althonos/proteinogenic/actions)
[![Codecov](https://img.shields.io/codecov/c/gh/althonos/proteinogenic/master.svg?style=flat-square&maxAge=600)](https://codecov.io/gh/althonos/proteinogenic)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/proteinogenic)
[![Crate](https://img.shields.io/crates/v/proteinogenic.svg?maxAge=600&style=flat-square)](https://crates.io/crates/proteinogenic)
[![Documentation](https://img.shields.io/badge/docs.rs-latest-4d76ae.svg?maxAge=2678400&style=flat-square)](https://docs.rs/proteinogenic)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/proteinogenic/blob/master/CHANGELOG.md)
[![GitHub issues](https://img.shields.io/github/issues/althonos/proteinogenic.svg?style=flat-square&maxAge=600)](https://github.com/althonos/proteinogenic/issues)


## 🔌 Usage

This crate builds on top of [`purr`](https://docs.rs/purr), a crate providing
primitives for reading and writing [SMILES].

Use the `AminoAcid` enum to encode the sequence residues, and build a SMILES
string with `proteinogenic::smiles`:

```rust
extern crate proteinogenic;

let sequence = "KGILGKLGVVQAGVDFVSGVWAGIKQSAKDHPNA";
let residues = sequence.chars()
  .map(|c| proteinogenic::AminoAcid::from_code1(c).unwrap());

let s = proteinogenic::smiles(residues);
```

This SMILES string can be used in conjunction with other cheminformatics toolkits,
for instance [OpenBabel](http://openbabel.org/wiki/Main_Page) which can generate a PNG figure:

![Skeletal formula of divergicin 750](https://raw.github.com/althonos/proteinogenic/master/static/divergicin.png)

Note that `proteinogenic` is not limited to building a SMILES string; it can
actually use any [`purr::walk::Follower`](https://docs.rs/purr/latest/purr/walk/trait.Follower.html)
implementor to generate an in-memory representation of a protein formula. If
your code is already compatible with `purr`, then you'll be able to use
protein sequences quite easily.

```rust
extern crate proteinogenic;
extern crate purr;

let sequence = "KGILGKLGVVQAGVDFVSGVWAGIKQSAKDHPNA";
let residues = sequence.chars()
  .map(|c| proteinogenic::AminoAcid::from_code1(c).unwrap());

let mut builder = purr::graph::Builder::new();
proteinogenic::visit(residues, &mut builder);

builder.build()
  .expect("failed to create a graph representation");
```

*The API is not yet stable, and may change to follow changes introduced by
`purr` or to improve the interface ergonomics.*

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/proteinogenic/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

<!-- ### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/proteinogenic/blob/main/CONTRIBUTING.md)
for more details. -->

## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/proteinogenic/blob/master/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.

## 🔍 See Also

If you're a bioinformatician and a Rustacean, you may be interested in these
other libraries:

- [`uniprot.rs`](https://github.com/althonos/uniprot.rs): Rust data structures
  for the UniProtKB databases.
- [`obofoundry.rs`](https://github.com/althonos/obofoundry.rs): Rust data
  structures for the OBO Foundry.
- [`fastobo`](https://github.com/fastobo/fastobo): Rust parser and abstract
  syntax tree for Open Biomedical Ontologies.
- [`pubchem.rs`](https://github.com/althonos/pubchem.rs): Rust data structures
  and API client for the PubChem API.

## 📜 License

This library is provided under the open-source
[MIT license](https://choosealicense.com/licenses/mit/).

*This project was developed by [Martin Larralde](https://github.com/althonos/)
during his PhD project at the [European Molecular Biology Laboratory](https://www.embl.de/)
in the [Zeller team](https://github.com/zellerlab).*
