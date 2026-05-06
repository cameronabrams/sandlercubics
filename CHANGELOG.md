# Changelog

All notable changes to sandlercubics will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.11.0]

### Added

- Output of enthalpy and entropy departures.

## [0.10.2] - 2026-02-07

### Fixed

- Minor bug fixes.

## [0.10.1] - 2026-02-04

### Changed

- Minor CLI and EOS refinements; quickstart documentation updated.

## [0.10.0] - 2026-02-04

### Added

- Integration with `sandlermisc`'s `ThermodynamicState` class.

## [0.9.2] - 2026-01-12

### Fixed

- Documentation configuration fixes.

## [0.9.1] - 2026-01-11

### Changed

- Documentation configuration and package metadata updates.

## [0.9.0] - 2026-01-11

### Added

- `.set_compound()` method now accepts compound name strings.
- Wider combinations of degrees of freedom.
- Output of two-phase saturated states.
- Direct calculation of absolute enthalpy and entropy (T ref 0 °C, P ref 0.1 MPa).
- Option to display both vapor and liquid phase properties simultaneously.

## [0.8.0] - 2026-01-10

### Changed

- Major refactoring of CLI and EOS implementations for improved extensibility.
- Added `__init__.py` public API surface.

## [0.6.0] - 2026-01-09

### Added

- Initial ReadTheDocs documentation setup.

## [0.5.1] - 2026-01-08

### Added

- Expanded CLI and EOS capabilities; initial test suite added.

## [0.4.0] - 2026-01-07

### Added

- Heat of vaporization (ΔH_vap) property calculation.
- Entropy of vaporization (ΔS_vap) property calculation.

## [0.3.0] - 2026-01-05

### Added

- Soave-Redlich-Kwong (SRK) equation of state.
- Comparative examples between vdW, PR, and SRK.

### Changed

- Refactored base EOS class for better extensibility.

## [0.2.1] - 2026-01-06

### Changed

- Redefined `CubicEOS` abstract base class; improved class hierarchy.

## [0.2.0] - 2025-12-31

### Added

- `StateReporter` class for formatted output.
- `delta` subcommand for property change calculations.
- Heat capacity integration for enthalpy and entropy changes.

## [0.1.1] - 2025-12-30

### Fixed

- Removed erroneous message in CLI output.

## [0.1.0] - 2025-12-30

### Added

- Initial release with van der Waals and Peng-Robinson equations of state.
- CLI `state` subcommand.
- Integration with sandlerprops database.
- Departure function calculations (Z, H_dep, S_dep).