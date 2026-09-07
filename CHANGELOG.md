# Changelog

All notable changes to this project will be documented in this file.

The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and [Element](https://github.com/vector-im/element-android) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

[//]: # (Available sections in changelog)
[//]: # (### API changes warning ⚠️:)
[//]: # (### Added Features and Improvements 🙌:)
[//]: # (### Bugfix 🐛:)
[//]: # (### Other changes:)


## [Unreleased]


## [1.0.0] - 2026-09-07
### Added Features and Improvements 🙌:
- Add compatibility with current supported Python versions.
- Improve packaging configuration and dependency handling.

### Other changes:
- Improve type annotations and linting.
- Update documentation and continuous-integration configuration.
- Add Zenodo release metadata.


## [0.4.1] - 2025-11-24
### Bugfix 🐛:
- Remove unused deps to fix publishing to conda, #32
- Fix unreachable code in input method


## [0.4.0] - 2025-11-19
### Added Features:
- New `FeatureSet` class for loading multidimensional trajectory features.
- added `ForceEstimator.memory_kernel` to analyse force correlations at a fixed reference time point.
- Added tutorials for both new features
  
### API changes warning ⚠️:
- Drop python 3.8 support

### Added Features and Improvements 🙌:
- Add python 3.11-3.14 support🎉
- Add tutorial explaining the theory

### Other changes:
- Fix newlines in docs
- Add generalized smoothing of estimators
- Many minor bugfixes and improvements in the docs
- Changed syntax of plotting functions


## [0.3.0] - 2023-03-28
### Added Features and Improvements 🙌:
- Initial public release 🎉
- Added documentation, including tutorials 🎉
- Cleaned-up source code


## [0.2.1] - 2022-12-12
### Added Features and Improvements 🙌:
- Beta candidate


## [0.2.0] - 2022-05-23
### Added Features and Improvements 🙌:
- Alpha candidate


[Unreleased]: https://github.com/moldyn/dcTMD/compare/v1.0.0...main
[1.0.0]: https://github.com/moldyn/dcTMD/compare/v0.4.1...v1.0.0
[0.4.1]: https://github.com/moldyn/dcTMD/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/moldyn/dcTMD/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/moldyn/dcTMD/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/moldyn/dcTMD/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/moldyn/dcTMD/tree/v0.2.0
