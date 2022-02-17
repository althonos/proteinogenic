# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]

[Unreleased]: https://github.com/althonos/proteinogenic/compare/v0.2.0...HEAD


## [v0.2.0] - 2022-02-17

[v0.2.0]: https://github.com/althonos/proteinogenic/compare/v0.1.0...v0.2.0

### Fixed
- Kekulization of imidazole cycle of `L-histidine` residues.

### Changed
- Refactored API to make failible operations return a result.
- Renamed `AminoAcid::from_code1` to `AminoAcid::from_char`.
- Renamed `AminoAcid::from_code3` to `AminoAcid::from_code`.

### Added
- `AminoAcid::as_code` to view the 3-letter code of an `AminoAcid` variant.
- `L-pyrrolysine`, `dehydroalanine` and `(Z)-dehydrobutyrine` amino acids.
- Support for cross-link modifications like cystine or lanthionine.
- Support for head-to-tail homodetic cyclization.
- Dedicated error type for the new possible errors.


## [v0.1.0] - 2022-01-15

[v0.1.0]: https://github.com/althonos/proteinogenic/compare/dfa86e6...v0.1.0

Initial release.
