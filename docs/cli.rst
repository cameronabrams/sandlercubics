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

Examples
~~~~~~~~

**Basic calculation:**  Methane at 400 K and 0.5 MPa using the Peng-Robinson EOS:

.. code-block:: bash

   sandlercubics state -T 400 -P 0.5 -eos pr -n methane

Output::

   State report for methane using Peng-Robinson Equation of State:
   T    = 400.00 K
   P    = 0.50 MPa
   Z    = 0.9964
   v    = 0.00662792 m3/mol
   h    = 3858.78 J/mol
   s    = -2.24 J/mol-K
   hdep = -54.75 J/mol
   sdep = -0.11 J/mol-K

   Constants used for calculations:
   Tc    = 190.40 K
   Pc    = 4.60 MPa
   omega = 0.011
   Tref  = 298.15 K
   Pref  = 0.10 MPa
   CpA   = 19.25 J/mol-K
   CpB   = 5.213e-02 J/mol-K^2
   CpC   = 1.197e-05 J/mol-K^3
   CpD   = -1.132e-08 J/mol-K^4

**Using van der Waals equation:**

.. code-block:: bash

   sandlercubics state -T 300 -P 1.0 -eos vdw -n ethane

**Manual critical properties:**  Also with a constant pressure heat capacity of 32.0 J/mol-K and no temperature-dependence of Cp:

.. code-block:: bash

   sandlercubics state -T 350 -P 5.0 -eos pr --Tc 190.4 --Pc 4.6 --omega 0.011 --Cp 32.0 0.0 0.0 0.0 -n methane

**Two-phase calculation:** Methane at 180 K and 3.0 MPa using Peng-Robinson EOS; this state is in the two-phase region, so both vapor and liquid solutions are returned, vapor first.

.. code-block:: bash

   sandlercubics state -T 180 -P 3.0 -eos pr -n methane

Output::

   State report for methane using Peng-Robinson Equation of State:
   T              = 180.00 K
   P              = 3.00 MPa
   Z              = 0.6259, 0.1244
   v              = 0.000312234, 6.20463e-05 m3/mol
   h              = -5407.84, -9334.59 J/mol
   s              = -50.63, -72.85 J/mol-K
   hdep           = -1597.87, -5524.61 J/mol
   sdep           = -6.22, -28.43 J/mol-K
   Pvap(180.00 K) = 3.33 MPa
   Hvap(180.00 K) = 3686.60 J/mol
   Svap(180.00 K) = 20.4811 J/mol-K
   Tsat(3.00 MPa) = 176.88 K

   Constants used for calculations:
   Tc    = 190.40 K
   Pc    = 4.60 MPa
   omega = 0.011
   Tref  = 298.15 K
   Pref  = 0.10 MPa
   CpA   = 19.25 J/mol-K
   CpB   = 5.213e-02 J/mol-K^2
   CpC   = 1.197e-05 J/mol-K^3
   CpD   = -1.132e-08 J/mol-K^4

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

**Basic state change:**

.. code-block:: bash

   sandlercubics delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr

Output::

   Property differences for methane using Peng-Robinson Equation of State:
   Delta H = 1571.86 J/mol
   Delta S = -1.45 J/mol-K
   Delta U = 1104.74 J/mol

   Constants used for calculations:
   Tc    = 190.40 K
   Pc    = 4.60 MPa
   omega = 0.011
   Tref  = 298.15 K
   Pref  = 0.10 MPa
   CpA   = 19.25 J/mol-K
   CpB   = 5.213e-02 J/mol-K^2
   CpC   = 1.197e-05 J/mol-K^3
   CpD   = -1.132e-08 J/mol-K^4

**With state details:**

.. code-block:: bash

   sandlercubics delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr --show-states

Output::

   State 1:
   T    = 350.00 K
   P    = 7.50 MPa
   Z    = 0.9262
   v    = 0.000359369 m3/mol
   h    = 929.35 J/mol
   s    = -32.10 J/mol-K
   hdep = -989.93 J/mol
   sdep = -2.13 J/mol-K

   State 2:
   T    = 400.00 K
   P    = 15.50 MPa
   Z    = 0.9509
   v    = 0.000204025 m3/mol
   h    = 2501.21 J/mol
   s    = -33.54 J/mol-K
   hdep = -1412.32 J/mol
   sdep = -2.86 J/mol-K


   Property differences for methane using Peng-Robinson Equation of State:
   Delta H = 1571.86 J/mol
   Delta S = -1.45 J/mol-K
   Delta U = 1104.74 J/mol

   Constants used for calculations:
   Tc    = 190.40 K
   Pc    = 4.60 MPa
   omega = 0.011
   Tref  = 298.15 K
   Pref  = 0.10 MPa
   CpA   = 19.25 J/mol-K
   CpB   = 5.213e-02 J/mol-K^2
   CpC   = 1.197e-05 J/mol-K^3
   CpD   = -1.132e-08 J/mol-K^4

Units Reference
---------------

The CLI uses SI units consistently unless otherwise specified. The following table summarizes the units used for input and output:

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

.. note::
   When using the Python API, pressure for EOS initialization is typically in **bar**, 
   not MPa. The CLI handles this conversion automatically.  The user can specify pressure units with the ``-pu`` option, the volume output units with the ``-vu`` option.

Tips and Best Practices
-----------------------

1. **Compound names**: Use lowercase names as they appear in the sandlerprops database (e.g., ``methane``, ``ethane``, ``water``)

2. **EOS selection**: 
   
   * Use ``ideal`` for low-pressure gases where non-ideal behavior is negligible
   * Use ``pr`` (Peng-Robinson) for most hydrocarbons and petroleum applications
   * Use ``srk`` (Soave-Redlich-Kwong) for gas processing
   * Use ``vdw`` (van der Waals) for educational purposes or as a baseline

3. **Phase specifications**: Two-phase conditions can have multiple solutions.  When two values are reported for Z, v, h, or s, the first corresponds to the vapor phase and the second to the liquid phase.

4. **Custom properties**: If a compound isn't in the database, you can manually specify Tc, Pc, omega, and Cp from any reliable source.

Error Messages
--------------

Common CLI errors and their solutions:

**"Compound not found"**
   The compound name isn't in the sandlerprops database. Either check the spelling or use manual property input (--Tc, --Pc, --omega, --Cp).

**"Convergence error"**
   The numerical solver couldn't find a solution. This can happen near critical points or at extreme conditions. Try slightly different T/P values.

