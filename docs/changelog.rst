.. _changelog:

Changelog
=========

All notable changes to sandlercubics will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

[0.5.0] - 2020-02-02
--------------------

Added
~~~~~

* Output of two-phase saturated states now allowed
* Enhanced phase identification for subcritical conditions
* Option to display both vapor and liquid phase properties simultaneously

Changed
~~~~~~~

* Improved numerical stability near saturation curve
* Better handling of edge cases in two-phase region

Fixed
~~~~~

* Resolved issues with phase selection in CLI
* Corrected saturation property calculations

[0.4.0] - 2019-XX-XX
--------------------

Added
~~~~~

* Heat of vaporization (ΔH\ :sub:`vap`) property calculation
* Entropy of vaporization (ΔS\ :sub:`vap`) property calculation
* New methods for calculating phase transition properties

Enhanced
~~~~~~~~

* Documentation for phase equilibrium calculations
* Examples showing vaporization property usage

[0.3.0] - 2019-XX-XX
--------------------

Added
~~~~~

* Soave-Redlich-Kwong (SRK) equation of state implementation
* Support for third cubic EOS in CLI
* Comparative examples between vdW, PR, and SRK

Changed
~~~~~~~

* Refactored base EOS class for better extensibility
* Improved code organization for adding new equations

[0.2.1] - 2019-XX-XX
--------------------

Changed
~~~~~~~

* Redefined ``CubicEOS`` abstract base class
* Improved class hierarchy and inheritance structure
* Better separation of concerns between base and derived classes

Fixed
~~~~~

* Type hints and method signatures in base class
* Documentation strings for inherited methods

[0.2.0] - 2019-XX-XX
--------------------

Added
~~~~~

* ``StateReporter`` class for formatted output
* ``delta`` subcommand in CLI for property change calculations
* Heat capacity integration for enthalpy and entropy changes
* Support for custom heat capacity coefficients

Changed
~~~~~~~

* CLI now uses consistent formatting across all output
* Improved readability of state reports

[0.1.1] - 2019-XX-XX
--------------------

Fixed
~~~~~

* Removed erroneous thank-you message in CLI output
* Minor bug fixes in output formatting

[0.1.0] - 2019-XX-XX
--------------------

Added
~~~~~

* Initial release of sandlercubics
* van der Waals equation of state implementation
* Peng-Robinson equation of state implementation
* Command-line interface with ``state`` subcommand
* Integration with sandlerprops database
* Basic property calculations (Z, v, H\ :sub:`dep`, S\ :sub:`dep`)
* Python API for programmatic use

Features
~~~~~~~~

* Cubic equation solver using analytical methods
* Departure function calculations
* Support for both vapor and liquid phases
* Manual property input option (Tc, Pc, omega)

Documentation
~~~~~~~~~~~~~

* README with installation and usage examples
* MIT License
* Basic package metadata

[Unreleased]
------------

Planned features for future releases:

Considering
~~~~~~~~~~~

* Mixture support (VLE for binary and multicomponent systems)
* Additional equations of state (Redlich-Kwong, etc.)
* Fugacity coefficient calculations
* Activity coefficient models
* Flash calculations
* Enhanced saturation curve calculations
* Improved performance for batch calculations
* Python type stubs for better IDE support
* Additional output formats (JSON, CSV)
* Interactive plotting capabilities
* Jupyter notebook examples
* More comprehensive test suite
* Continuous integration setup

Breaking Changes Policy
-----------------------

sandlercubics follows semantic versioning:

* **MAJOR** version (X.0.0): Incompatible API changes
* **MINOR** version (0.X.0): New functionality, backwards-compatible
* **PATCH** version (0.0.X): Bug fixes, backwards-compatible

Deprecation Warnings
--------------------

No features are currently deprecated.

When features are deprecated, they will:

1. Remain functional in current MAJOR version
2. Issue ``DeprecationWarning`` when used
3. Include migration instructions in warning message
4. Be documented in this changelog
5. Be removed in next MAJOR version

Upgrade Guide
-------------

From 0.4.x to 0.5.0
~~~~~~~~~~~~~~~~~~~

No breaking changes. All existing code should continue to work.

New feature usage:

.. code-block:: python

   # Two-phase output now available in CLI
   sandlercubics state -T 180 -P 3.0 -eos pr -n methane --phase both

From 0.3.x to 0.4.0
~~~~~~~~~~~~~~~~~~~

No breaking changes. New vaporization properties are optional.

From 0.2.x to 0.3.0
~~~~~~~~~~~~~~~~~~~

New SRK equation available:

.. code-block:: python

   from sandlercubics.eos import SoaveRedlichKwongEOS
   
   eos = SoaveRedlichKwongEOS(Tc, Pc, omega)

From 0.1.x to 0.2.0
~~~~~~~~~~~~~~~~~~~

* CLI now requires explicit subcommand (``state`` or ``delta``)
* Old: ``sandlercubics -T 400 -P 0.5 ...``
* New: ``sandlercubics state -T 400 -P 0.5 ...``

Contributing to Changelog
--------------------------

When contributing changes:

1. Add entry under ``[Unreleased]`` section
2. Use appropriate subsection (Added/Changed/Deprecated/Removed/Fixed/Security)
3. Write concise, user-focused descriptions
4. Link to relevant issue/PR numbers
5. Maintainers will organize entries during release

Example entry::

   [Unreleased]
   ------------
   
   Added
   ~~~~~
   * New equation of state: Benedict-Webb-Rubin (#123)
   * Support for temperature-dependent viscosity calculations
   
   Fixed
   ~~~~~
   * Incorrect enthalpy calculation near critical point (#456)

See Also
--------

* :doc:`contributing` - How to contribute
* `GitHub Releases <https://github.com/cameronabrams/sandlercubics/releases>`_ - Release notes
* `GitHub Issues <https://github.com/cameronabrams/sandlercubics/issues>`_ - Bug reports and feature requests
