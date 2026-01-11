.. _quickstart:

Quick Start Guide
=================

This guide will get you up and running with sandlercubics in just a few minutes.

Basic Concepts
--------------

sandlercubics provides implementations of pure-species cubic equations of state (EOS) presented in Sandler's textbook:

* **van der Waals (vdW)**: The original cubic EOS, useful for understanding basic concepts
* **Peng-Robinson (PR)**: Industry-standard EOS with good accuracy for hydrocarbons
* **Soave-Redlich-Kwong (SRK)**: Alternative cubic EOS popular in gas processing

You can set up a state using temperature (T) and pressure (P) as inputs, along with a compound name from the sandlerprops database or by providing critical properties directly.  Instead of T and P, you can also specify other combinations of state variables, such as T and molar volume (v), or P and molar enthalpy (h).  Any two independent state variables are sufficient to define the thermodynamic state of a pure substance.  However, at least one of the two variables must be temperature or pressure (for now).

The following computed properties are always represented by numpy arrays, even if they contain only a single value:

* Compressibility factor (Z)
* Molar enthalpy (h), entropy (s), volume (v), and internal energy (u) (with respect to reference temperature of 298.15 K and reference pressure of 0.1 MPa)
* Enthalpy and entropy departure (H\ :sub:`dep`, S\ :sub:`dep`)

The following properties are available when applicable as scalars:

* Heat of vaporization (ΔH\ :sub:`vap`) at state temperature and saturation pressure
* Entropy of vaporization (ΔS\ :sub:`vap`) at state temperature and saturation pressure
* Vapor pressure (P\ :sup:`vap`) (aka, saturation pressure) at state temperature
* Saturation temperature (T\ :sup:`sat`) at state pressure
* Vapor fraction (x\ :sub:`vap`) for two-phase states

Temperature is in Kelvin (K) and pressure is in megapascals (MPa) by default.  Other units of pressure are supported.  See the :doc:`API` documentation for details.

First Calculation
-----------------

Let's calculate the molar volume of methane at 400 K and 0.5 MPa using the Peng-Robinson equation:

From the Command Line
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   sandlercubics state -T 400 -P 0.5 -eos pr -n methane

Output::

   State report for methane using Peng-Robinson Equation of State:
   T    =  400.00 K
   P    =  0.50 MPa
   Z    =  0.996444
   v    =  0.00662792 m3/mol
   h    =  3858.78 J/mol
   s    = -2.23817 J/mol-K
   hdep = -54.7512 J/mol
   sdep = -0.107042 J/mol-K

From Python
~~~~~~~~~~~

.. code-block:: python

   from sandlercubics import PengRobinsonEOS
   # Create EOS object and set state
   pr_methane = PengRobinsonEOS(T=400, P=0.5).set_compound('methane')
   
   # Access calculated properties; remember these are numpy arrays
   print(f"Compressibility factor: {', '.join([f'{z: 4g}' for z in pr_methane.Z])}")
   print(f"Molar volume: {', '.join([f'{v: 6g}' for v in pr_methane.v])} m³/mol")
   print(f"Enthalpy: {', '.join([f'{h: 6g}' for h in pr_methane.h])} J/mol")
   print(f"Entropy: {', '.join([f'{s: 6g}' for s in pr_methane.s])} J/mol-K")
   print(f"Enthalpy departure: {', '.join([f'{hdep: 6g}' for hdep in pr_methane.h_departure])} J/mol")
   print(f"Entropy departure: {', '.join([f'{sdep: 6g}' for sdep in pr_methane.s_departure])} J/mol-K")

You need not set the state variables (T, P) during initialization; you can also just directly assign them later:

.. code-block:: python

   from sandlercubics import PengRobinsonEOS
   # Create EOS object without state and without compound
   pr = PengRobinsonEOS()
   
   # Set compound later; .set_compound returns self for convenience
   pr_methane = pr.set_compound('methane')
   
   # Set state later
   pr_methane.T = 400
   pr_methane.P = 0.5
   
   # Now access properties as before
   print(f"Molar volume: {', '.join([f'{v: 6g}' for v in pr_methane.v])} m³/mol")

State-Change Calculations
-------------------------

Calculate property changes between two thermodynamic states:

From the Command Line
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   sandlercubics delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr --show-states

Output::

   State-change calculations for methane using Peng-Robinson Equation of State:

   State 1:                       State 2:
   T    =  350.00 K               T    =  400.00 K
   P    =  7.50 MPa               P    =  15.50 MPa
   Z    =  0.92619                Z    =  0.950871
   v    =  0.000359369 m3/mol     v    =  0.000204025 m3/mol
   h    =  929.35 J/mol           h    =  2501.21 J/mol
   s    = -32.095 J/mol-K         s    = -33.5449 J/mol-K
   hdep = -989.935 J/mol          hdep = -1412.32 J/mol
   sdep = -2.12621 J/mol-K        sdep = -2.86197 J/mol-K

   Property changes:
   Δh =  1571.86 J/mol
   Δs = -1.44983 J/mol-K
   Δu =  1104.74 J/mol

From Python
~~~~~~~~~~~

State-change calculations can be performed by creating two EOS objects representing the initial and final states, then using the built-in methods to compute property differences:

.. code-block:: python

   from sandlercubics import PengRobinsonEOS
   
   # State 1
   pr_methane1 = PengRobinsonEOS(T=350, P=7.5).set_compound('methane')
   
   # State 2
   pr_methane2 = PengRobinsonEOS(T=400, P=15.5).set_compound('methane')
   
   # Calculate property differences using built-in methods
   delta_h = pr_methane1.delta_h(pr_methane2)
   delta_s = pr_methane1.delta_s(pr_methane2)
   delta_u = pr_methane1.delta_u(pr_methane2)
   
   print(f"Δh = {', '.join([f'{dh: 7g}' for dh in delta_h])} J/mol")
   print(f"Δs = {', '.join([f'{ds: 7g}' for ds in delta_s])} J/mol-K")
   print(f"Δu = {', '.join([f'{du: 7g}' for du in delta_u])} J/mol")

Available Equations of State
-----------------------------

Ideal gas
~~~~~~~~~~~~~

.. code-block:: python

   from sandlercubics import IdealGasEOS
   
   ig = IdealGasEOS()

van der Waals
~~~~~~~~~~~~~

.. code-block:: python

   from sandlercubics import VanDerWaalsEOS
   
   vdw_methane = VanDerWaalsEOS().set_compound('methane')

Peng-Robinson
~~~~~~~~~~~~~

.. code-block:: python

   from sandlercubics import PengRobinsonEOS
   
   pr_benzene = PengRobinsonEOS().set_compound('benzene')

Soave-Redlich-Kwong
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from sandlercubics import SoaveRedlichKwongEOS
   
   srk_ethanol = SoaveRedlichKwongEOS().set_compound('ethanol')

Next Steps
----------

* Learn more about the :doc:`cli` for advanced command-line usage
* Explore :doc:`examples` for more complex calculations
* Read about the :doc:`theory` behind cubic equations of state
* Check the :doc:`api/API` for complete API documentation
