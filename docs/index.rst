.. sandlercubics documentation master file

sandlercubics
=============

.. image:: https://img.shields.io/pypi/v/sandlercubics.svg
   :target: https://pypi.org/project/sandlercubics/
   :alt: PyPI version

.. image:: https://img.shields.io/pypi/pyversions/sandlercubics.svg
   :target: https://pypi.org/project/sandlercubics/
   :alt: Python versions

Digitized cubic equations of state from Sandler's 5th edition
--------------------------------------------------------------

**sandlercubics** implements a Python interface to the cubic equations of state found in *Chemical, Biochemical, and Engineering Thermodynamics* (5th edition) by Stan Sandler (Wiley, USA). 

.. warning::
   This package should be used for **educational purposes only**.

Features
--------

* **Generalized van der Waals equation of state** for pure substances
* **Peng-Robinson equation of state** for pure substances
* **Soave-Redlich-Kwong equation of state** for pure substances
* Command-line interface for quick calculations
* Python API for integration into larger workflows
* Calculation of thermodynamic properties including:
  
  * Compressibility factor (Z)
  * Molar volume (v), molar enthalpy (h), molar entropy (s), and molar internal energy (u)
  * Enthalpy and entropy departure (H\ :sub:`dep`, S\ :sub:`dep`)
  * Vapor pressure and saturation temperature
  * Heat and entropy of vaporization

Quick Start
-----------

Installation::

   pip install sandlercubics

Basic usage from the command line::

   sandlercubics state -T 400 -P 0.5 -eos pr -n methane

Basic usage from Python:

.. code-block:: python

   from sandlercubics import PengRobinsonEOS
   eos = PengRobinsonEOS(T=400, P=0.5).set_compound('methane')
   print(f"Molar volume: {', '.join(f'{v: 6g}' for v in eos.v)} m³/mol")

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   
   installation
   quickstart
   cli
   examples
   theory

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   
   api/API.rst

.. toctree::
   :maxdepth: 1
   :caption: Developer Guide
   
   contributing
   changelog

License
=======

This project is licensed under the MIT License - see the LICENSE file for details.

Citation
========

If you use this package for academic work, please cite:

* Sandler, S. (2017). *Chemical, Biochemical, and Engineering Thermodynamics* (5th ed.). Wiley.

Contact
=======

Cameron F. Abrams - cfa22@drexel.edu

GitHub: https://github.com/cameronabrams/sandlercubics
