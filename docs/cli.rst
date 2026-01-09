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

Shows general help and lists available subcommands.

.. code-block:: bash

   sandlercubics --version

Displays the installed version of sandlercubics.

``state`` Command
-----------------

Calculate thermodynamic properties for a pure substance at a specified state.

Syntax
~~~~~~

.. code-block:: bash

   sandlercubics state -T TEMP -P PRESSURE -eos EOS_TYPE -n COMPOUND [OPTIONS]

Required Arguments
~~~~~~~~~~~~~~~~~~

.. option:: -T, --temperature TEMP

   Temperature in Kelvin (K)

.. option:: -P, --pressure PRESSURE

   Pressure in megapascals (MPa)

.. option:: -eos, --equation-of-state EOS_TYPE

   Equation of state to use. Options:
   
   * ``vdw``: van der Waals
   * ``pr``: Peng-Robinson
   * ``srk``: Soave-Redlich-Kwong

.. option:: -n, --name COMPOUND

   Compound name (must be available in sandlerprops database)

Optional Arguments
~~~~~~~~~~~~~~~~~~

.. option:: --Tc CRITICAL_TEMP

   Override critical temperature (K) instead of using database value

.. option:: --Pc CRITICAL_PRESSURE

   Override critical pressure (MPa) instead of using database value

.. option:: --omega ACENTRIC_FACTOR

   Override acentric factor instead of using database value

.. option:: --phase PHASE

   Specify phase for two-phase conditions. Options:
   
   * ``vapor`` or ``v``: Vapor phase
   * ``liquid`` or ``l``: Liquid phase
   * ``both``: Show both phases (for saturated states)

Examples
~~~~~~~~

**Basic calculation:**

.. code-block:: bash

   sandlercubics state -T 400 -P 0.5 -eos pr -n methane

Output::

   EOS  = pr
   T    = 400.00 K
   P    = 0.50 MPa
   Z    = 1.00
   v    = 0.006628 m³/mol
   Hdep = -54.75 J/mol
   Sdep = -0.11 J/mol-K
   Tc   = 190.40 K
   Pc   = 4.60 MPa
   omega = 0.011

**Using van der Waals equation:**

.. code-block:: bash

   sandlercubics state -T 300 -P 1.0 -eos vdw -n ethane

**Manual critical properties:**

.. code-block:: bash

   sandlercubics state -T 350 -P 5.0 -eos pr --Tc 190.4 --Pc 4.6 --omega 0.011

**Two-phase calculation:**

.. code-block:: bash

   sandlercubics state -T 180 -P 3.0 -eos pr -n methane --phase both

``delta`` Command
-----------------

Calculate changes in thermodynamic properties between two states.

Syntax
~~~~~~

.. code-block:: bash

   sandlercubics delta -T1 TEMP1 -P1 PRESS1 -T2 TEMP2 -P2 PRESS2 -eos EOS_TYPE -n COMPOUND [OPTIONS]

Required Arguments
~~~~~~~~~~~~~~~~~~

.. option:: -T1, --temperature1 TEMP1

   Initial temperature (K)

.. option:: -P1, --pressure1 PRESS1

   Initial pressure (MPa)

.. option:: -T2, --temperature2 TEMP2

   Final temperature (K)

.. option:: -P2, --pressure2 PRESS2

   Final pressure (MPa)

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

.. option:: --CpA, --CpB, --CpC, --CpD COEFFICIENTS

   Heat capacity polynomial coefficients. The ideal gas heat capacity is calculated as:
   
   .. math::
   
      C_p^{ig} = A + BT + CT^2 + DT^3

Examples
~~~~~~~~

**Basic state change:**

.. code-block:: bash

   sandlercubics delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr

Output::

   Property differences:
   Delta H = 2416.63 J/mol
   Delta S = 0.02 J/mol-K
   Delta U = 2883.75 J/mol

**With state details:**

.. code-block:: bash

   sandlercubics delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr --show-states

Output::

   State 1:
   EOS  = pr
   T    = 350.00 K
   P    = 7.50 MPa
   Z    = 0.93
   v    = 0.000359 m³/mol
   Hdep = -989.93 J/mol
   Sdep = -2.13 J/mol-K

   State 2:
   EOS  = pr
   T    = 400.00 K
   P    = 15.50 MPa
   Z    = 0.95
   v    = 0.000204 m³/mol
   Hdep = -1412.32 J/mol
   Sdep = -2.86 J/mol-K

   Property differences:
   Delta H = 2416.63 J/mol
   Delta S = 0.02 J/mol-K
   Delta U = 2883.75 J/mol

   Constants used for calculations:
   Tc    = 190.40 K
   Pc    = 4.60 MPa
   omega = 0.011
   CpA   = 19.25 J/mol-K
   CpB   = 5.213e-02 J/mol-K²
   CpC   = 1.197e-05 J/mol-K³
   CpD   = -1.132e-08 J/mol-K⁴

Units Reference
---------------

The CLI uses SI units consistently:

==================  ============
Property            Unit
==================  ============
Temperature (T)     Kelvin (K)
Pressure (P)        Megapascal (MPa)
Molar volume (v)    m³/mol
Enthalpy (H)        J/mol
Entropy (S)         J/(mol·K)
Internal energy (U) J/mol
==================  ============

.. note::
   When using the Python API, pressure for EOS initialization is typically in **bar**, 
   not MPa. The CLI handles this conversion automatically.

Tips and Best Practices
-----------------------

1. **Compound names**: Use lowercase names as they appear in the sandlerprops database (e.g., ``methane``, ``ethane``, ``water``)

2. **EOS selection**: 
   
   * Use ``pr`` (Peng-Robinson) for most hydrocarbons and petroleum applications
   * Use ``srk`` (Soave-Redlich-Kwong) for gas processing
   * Use ``vdw`` (van der Waals) for educational purposes or as a baseline

3. **Phase specifications**: Two-phase conditions can have multiple solutions. Use ``--phase`` to specify which phase you want when near saturation conditions.

4. **Custom properties**: If a compound isn't in the database, you can manually specify Tc, Pc, and omega from any reliable source.

Error Messages
--------------

Common CLI errors and their solutions:

**"Compound not found"**
   The compound name isn't in the sandlerprops database. Either check the spelling or use manual property input (--Tc, --Pc, --omega).

**"Convergence error"**
   The numerical solver couldn't find a solution. This can happen near critical points or at extreme conditions. Try slightly different T/P values.

**"Two-phase region detected"**
   The state point is in the two-phase region. Use ``--phase`` to specify liquid or vapor.

Scripting with the CLI
----------------------

The CLI output is designed to be parsable. You can use standard Unix tools to extract values:

.. code-block:: bash

   # Extract just the molar volume
   sandlercubics state -T 400 -P 0.5 -eos pr -n methane | grep "^v" | awk '{print $3}'
   
   # Run calculations for multiple pressures
   for P in 0.5 1.0 2.0 5.0 10.0; do
       echo "P = $P MPa:"
       sandlercubics state -T 400 -P $P -eos pr -n methane | grep "^Z"
   done

For more complex workflows, consider using the Python API instead.
