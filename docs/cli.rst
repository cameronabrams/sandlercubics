.. _cli:

Command-Line Interface
======================

sandlercubics provides a command-line interface (CLI) for quick thermodynamic calculations without writing Python code.

Overview
--------

The CLI has two main subcommands:

* ``state``: Calculate properties at a single thermodynamic state
* ``delta``: Calculate property changes between two states

Global Options
--------------

.. code-block:: bash

   sandlercubics --help

Shows general help and lists available subcommands. (Note that ``-h`` is reserved for enthalpy inputs.)  

.. code-block:: bash

   sandlercubics --version

Displays the installed version of sandlercubics.

``state`` Command
-----------------

Calculate thermodynamic properties for a pure substance at a specified state.

Syntax
~~~~~~

.. code-block:: bash

   sandlercubics state -eos EOS_TYPE -n COMPOUND [OPTIONS]

Required Arguments
~~~~~~~~~~~~~~~~~~

.. option:: -eos, --equation-of-state EOS_TYPE

   Equation of state to use. Options:
   
   * ``ideal``: Ideal Gas
   * ``vdw``: van der Waals
   * ``pr``: Peng-Robinson
   * ``srk``: Soave-Redlich-Kwong

.. option:: -n, --name COMPOUND

   Compound name (must be available in sandlerprops database)

Optional Arguments
~~~~~~~~~~~~~~~~~~

Two of the following state variables must be provided, at least one of which must be temperature or pressure:

.. option:: -T, --temperature TEMP

   Temperature (K)

.. option:: -P, --pressure PRESS

   Pressure (MPa)

.. option:: -u, --internal-energy U

   Internal energy (J/mol)

.. option:: -v, --molar-volume V

   Molar volume (m³/mol)

.. option:: -h, --enthalpy H

   Enthalpy (J/mol)

.. option:: -s, --entropy S

   Entropy (J/mol-K)

These are optional arguments:

.. option:: -pu , --pressure-unit PRESSURE_UNIT

   Pressure unit for input and output (default: MPa). Options: MPa, bar, atm

.. option:: -vu, --volume-unit VOLUME_UNIT

   Volume unit for output (default: m3). Options: m3, L, cm3

.. option:: --Tc CRITICAL_TEMP

   Override critical temperature (K) instead of using database value

.. option:: --Pc CRITICAL_PRESSURE

   Override critical pressure (MPa) instead of using database value

.. option:: --omega ACENTRIC_FACTOR

   Override acentric factor instead of using database value

.. option:: --Cp [CpA] [CpB] [CpC] [CpD]

   Heat capacity polynomial coefficients, used if the component is not specified or if you want to override database values. The ideal gas heat capacity is calculated as:
   
   .. math::
   
      C_p^{ig} = A + BT + CT^2 + DT^3

.. option:: --show-props
   
   Display critical constants and heat capacity used in calculations (default is to hide)

Examples
~~~~~~~~

**Basic calculation:**  Methane at 400 K and 0.5 MPa using the Peng-Robinson EOS:

.. code-block:: bash

   sandlercubics state -T 400 -P 0.5 -eos pr -n methane --show-props

Output::

   T     =   400 kelvin
   P     =   0.5 megapascal
   v     =  0.00662792 meter ** 3 / mole
   s     = -2.23817 joule / kelvin / mole
   h     =  3858.78 joule / mole
   u     =  544.821 joule / mole
   Pv    =  3313.96 joule / mole
   Z     =  0.996444
   Tsat  =  135.11 kelvin at 0.5 megapascal
   Tc    =  190.4 kelvin
   Pc    =    4.6 megapascal
   omega =  0.011
   CpA   =  19.25
   CpB   =  0.05213
   CpC   =  1.197e-05
   CpD   = -1.132e-08

**Using van der Waals equation:**

.. code-block:: bash

   sandlercubics state -T 300 -P 1.0 -eos vdw -n ethane

**Manual critical properties:**  Also with a constant pressure heat capacity of 32.0 J/mol-K and no temperature-dependence of Cp:

.. code-block:: bash

   sandlercubics state -T 350 -P 5.0 -eos pr --Tc 190.4 --Pc 4.6 --omega 0.011 --Cp 32.0 0.0 0.0 0.0 -n methane

**Two-phase calculation:** Methane at 180 K and 3.0 MPa using Peng-Robinson EOS; this state is in the two-phase region, so both vapor and liquid solutions are returned, vapor first.

.. code-block:: bash

   sandlercubics state -T 180 -P 3.0 -eos pr -n methane --show-props

Output::

   T     =   180 kelvin
   P     =     3 megapascal
   v     =  0.000312234 meter ** 3 / mole
   s     = -50.6319 joule / kelvin / mole
   h     = -5407.84 joule / mole
   u     = -6344.55 joule / mole
   Pv    =  936.703 joule / mole
   Z     =  0.625886
   Pvap  =  3.32679 megapascal at 180 kelvin
   Hvap  =  3686.6 joule / mole at 180 kelvin
   Svap  =  20.4811 joule / kelvin / mole at 180 kelvin
   Tsat  =  176.88 kelvin at 3 megapascal
   Tc    =  190.4 kelvin
   Pc    =    4.6 megapascal
   omega =  0.011
   CpA   =  19.25
   CpB   =  0.05213
   CpC   =  1.197e-05
   CpD   = -1.132e-08

``delta`` Command
-----------------

Calculate changes in thermodynamic properties between two states.

Syntax
~~~~~~

.. code-block:: bash

   sandlercubics delta -eos EOS_TYPE -n COMPOUND [OPTIONS]

Required Arguments
~~~~~~~~~~~~~~~~~~

.. option:: -eos, --equation-of-state EOS_TYPE

   Equation of state (``vdw``, ``pr``, or ``srk``)

.. option:: -n, --name COMPOUND

   Compound name

Optional Arguments
~~~~~~~~~~~~~~~~~~

.. option:: --show-states

   Display full state information for both states in addition to property differences

.. option:: --Tc CRITICAL_TEMP

   Override critical temperature (K)

.. option:: --Pc CRITICAL_PRESSURE

   Override critical pressure (MPa)

.. option:: --omega ACENTRIC_FACTOR

   Override acentric factor

.. option:: --Cp [CpA] [CpB] [CpC] [CpD]

   Heat capacity polynomial coefficients, if the component is not specified or if you want to override database values.

Examples
~~~~~~~~

**Basic state change calculations:**

.. code-block:: bash

   sandlercubics delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr

Output::

State-change calculations for methane using Peng-Robinson Equation of State:
ΔT  =     50 kelvin
ΔP  =      8 megapascal
Δh  =  1571.86 joule / mole
Δs  = -1.44983 joule / kelvin / mole
Δu  =  1104.74 joule / mole
Δv  = -0.000155344 meter ** 3 / mole
ΔPv =  467.124 joule / mole
ΔZ  =  0.0246816

**With state details:**

.. code-block:: bash

   sandlercubics delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr --show-states

Output::

   State-change calculations for methane using Peng-Robinson Equation of State:

   State 1:                                 State 2:
   T  =   350 kelvin                        T  =   400 kelvin
   P  =   7.5 megapascal                    P  =  15.5 megapascal
   v  =  0.000359369 meter ** 3 / mole      v  =  0.000204025 meter ** 3 / mole
   s  = -32.095 joule / kelvin / mole       s  = -33.5449 joule / kelvin / mole
   h  =  929.35 joule / mole                h  =  2501.21 joule / mole
   u  = -1765.92 joule / mole               u  = -661.182 joule / mole
   Pv =  2695.27 joule / mole               Pv =  3162.39 joule / mole
   Z  =  0.92619                            Z  =  0.950871

   Property changes:
   ΔT  =     50 kelvin
   ΔP  =      8 megapascal
   Δh  =  1571.86 joule / mole
   Δs  = -1.44983 joule / kelvin / mole
   Δu  =  1104.74 joule / mole
   Δv  = -0.000155344 meter ** 3 / mole
   ΔPv =  467.124 joule / mole
   ΔZ  =  0.0246816

Units Reference
---------------

The CLI uses SI units consistently unless otherwise specified, implemented by the pint package. The following table summarizes the units used for input and output:

=================== ==================
Property            Unit
=================== ==================
Temperature (T)     Kelvin (K)
Pressure (P)        Megapascal (MPa)
Molar volume (v)    m³/mol
Enthalpy (H)        J/mol
Entropy (S)         J/(mol-K)
Internal energy (U) J/mol
=================== ==================

