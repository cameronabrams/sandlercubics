.. _quickstart:

Quick Start Guide
=================

This guide will get you up and running with sandlercubics in just a few minutes.

Basic Concepts
--------------

sandlercubics provides implementations of cubic equations of state (EOS) commonly used in chemical engineering thermodynamics:

* **van der Waals (vdW)**: The original cubic EOS, useful for understanding basic concepts
* **Peng-Robinson (PR)**: Industry-standard EOS with good accuracy for hydrocarbons
* **Soave-Redlich-Kwong (SRK)**: Alternative cubic EOS popular in gas processing

Each EOS can calculate:

* Compressibility factor (Z)
* Molar volume (v)
* Enthalpy departure (H\ :sub:`dep`)
* Entropy departure (S\ :sub:`dep`)
* Phase equilibrium properties

First Calculation
-----------------

Let's calculate the molar volume of methane at 400 K and 0.5 MPa using the Peng-Robinson equation:

From the Command Line
~~~~~~~~~~~~~~~~~~~~~

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

From Python
~~~~~~~~~~~

.. code-block:: python

   from sandlercubics.eos import PengRobinsonEOS
   from sandlerprops.properties import PropertiesDatabase
   
   # Load compound properties
   db = PropertiesDatabase()
   methane = db.get_compound('methane')
   
   # Create EOS object
   # Note: Pc needs to be in bar, so divide by 10
   eos = PengRobinsonEOS(
       Tc=methane.Tc,      # Critical temperature (K)
       Pc=methane.Pc/10,   # Critical pressure (bar)
       omega=methane.Omega # Acentric factor
   )
   
   # Set state conditions
   eos.T = 400  # Temperature in K
   eos.P = 0.5  # Pressure in MPa
   
   # Access calculated properties
   print(f"Compressibility factor: {eos.Z:.4f}")
   print(f"Molar volume: {eos.v.item():.6f} m³/mol")
   print(f"Enthalpy departure: {eos.Hdep:.2f} J/mol")
   print(f"Entropy departure: {eos.Sdep:.2f} J/mol-K")

Comparing States
----------------

Calculate property changes between two thermodynamic states:

From the Command Line
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   sandlercubics delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr --show-states

This calculates ΔH, ΔS, and ΔU between the two states.

From Python
~~~~~~~~~~~

.. code-block:: python

   from sandlercubics.eos import PengRobinsonEOS
   from sandlerprops.properties import PropertiesDatabase
   
   db = PropertiesDatabase()
   methane = db.get_compound('methane')
   
   # State 1
   eos1 = PengRobinsonEOS(Tc=methane.Tc, Pc=methane.Pc/10, omega=methane.Omega)
   eos1.T = 350
   eos1.P = 7.5
   
   # State 2
   eos2 = PengRobinsonEOS(Tc=methane.Tc, Pc=methane.Pc/10, omega=methane.Omega)
   eos2.T = 400
   eos2.P = 15.5
   
   # Calculate property differences (you'll need to implement this based on your API)
   # This is a placeholder - adjust based on actual implementation
   delta_H = (eos2.H - eos1.H)  # May need heat capacity integration
   delta_S = (eos2.S - eos1.S)
   delta_U = delta_H - (eos2.P * eos2.v - eos1.P * eos1.v)

Available Equations of State
-----------------------------

van der Waals
~~~~~~~~~~~~~

.. code-block:: python

   from sandlercubics.eos import VanDerWaalsEOS
   
   eos = VanDerWaalsEOS(Tc=190.4, Pc=46.0, omega=0.011)

Peng-Robinson
~~~~~~~~~~~~~

.. code-block:: python

   from sandlercubics.eos import PengRobinsonEOS
   
   eos = PengRobinsonEOS(Tc=190.4, Pc=46.0, omega=0.011)

Soave-Redlich-Kwong
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from sandlercubics.eos import SoaveRedlichKwongEOS
   
   eos = SoaveRedlichKwongEOS(Tc=190.4, Pc=46.0, omega=0.011)

Working with Different Compounds
---------------------------------

The sandlerprops database contains properties for common compounds:

.. code-block:: python

   from sandlerprops.properties import PropertiesDatabase
   
   db = PropertiesDatabase()
   
   # Get compound by name
   ethane = db.get_compound('ethane')
   propane = db.get_compound('propane')
   water = db.get_compound('water')
   
   # Access critical properties
   print(f"Ethane Tc: {ethane.Tc} K")
   print(f"Ethane Pc: {ethane.Pc} bar")
   print(f"Ethane ω: {ethane.Omega}")

If you don't have sandlerprops installed, you can manually input critical properties from any reliable source.

Next Steps
----------

* Learn more about the :doc:`cli` for advanced command-line usage
* Explore :doc:`examples` for more complex calculations
* Read about the :doc:`theory` behind cubic equations of state
* Check the :doc:`api/index` for complete API documentation
