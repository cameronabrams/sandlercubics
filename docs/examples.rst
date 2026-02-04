.. _examples:

Examples
========

This page provides practical examples of using sandlercubics for common thermodynamic calculations.

Basic Property Calculations
----------------------------

Example 1: Single State Point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate properties of propane at 350 K and 2.0 MPa using Peng-Robinson:

.. code-block:: python

   from sandlercubics import PengRobinsonEOS
   
   # Initialize EOS and set state
   eos = PengRobinsonEOS(T=350, P=2.0).set_compound('propane')
   
   # Display results
   print(f"State conditions:")
   print(f"  T = {eos.T} K")
   print(f"  P = {eos.P} MPa")
   print(f"\nCalculated properties:")
   print(f"  Z    = {eos.Z:.4f}")
   print(f"  v    = {eos.v.item():.6f} m³/mol")
   print(f"  Hdep = {eos.Hdep:.2f} J/mol")
   print(f"  Sdep = {eos.Sdep:.4f} J/(mol·K)")

Example 2: Comparing Equations of State
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare different equations for the same state:

.. code-block:: python

   from sandlercubics.eos import VanDerWaalsEOS, PengRobinsonEOS, SoaveRedlichKwongEOS
   from sandlerprops.properties import PropertiesDatabase
   
   db = PropertiesDatabase()
   ethane = db.get_compound('ethane')
   
   # Same critical properties for all equations
   T, P = 300, 5.0  # K, MPa
   Tc, Pc, omega = ethane.Tc, ethane.Pc / 10, ethane.Omega
   
   # Create EOS objects
   equations = {
       'van der Waals': VanDerWaalsEOS(Tc, Pc, omega),
       'Peng-Robinson': PengRobinsonEOS(Tc, Pc, omega),
       'Soave-RK': SoaveRedlichKwongEOS(Tc, Pc, omega)
   }
   
   # Compare results
   print(f"Ethane at T = {T} K, P = {P} MPa\n")
   print(f"{'Equation':<20} {'Z':>10} {'v (m³/mol)':>15} {'Hdep (J/mol)':>15}")
   print("-" * 65)
   
   for name, eos in equations.items():
       eos.T = T
       eos.P = P
       print(f"{name:<20} {eos.Z:>10.4f} {eos.v.item():>15.6f} {eos.Hdep:>15.2f}")

Phase Behavior
--------------

Example 3: Vapor-Liquid Equilibrium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Identify phases in a two-phase region:

.. code-block:: python

   from sandlercubics.eos import PengRobinsonEOS
   from sandlerprops.properties import PropertiesDatabase
   
   db = PropertiesDatabase()
   propane = db.get_compound('propane')
   
   eos = PengRobinsonEOS(
       Tc=propane.Tc,
       Pc=propane.Pc / 10,
       omega=propane.Omega
   )
   
   # Subcritical conditions (two-phase possible)
   eos.T = 300  # K (Tc of propane is ~370 K)
   eos.P = 1.0  # MPa
   
   # Check if we're in two-phase region
   # You may need to implement phase checking based on your API
   if hasattr(eos, 'phases'):
       for phase in eos.phases:
           print(f"{phase.name} phase:")
           print(f"  Z = {phase.Z:.4f}")
           print(f"  v = {phase.v:.6f} m³/mol")
           print(f"  Hdep = {phase.Hdep:.2f} J/mol")
           print()


Property Changes
----------------

Example 5: Isothermal Compression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate work and property changes during isothermal compression:

.. code-block:: python

   import numpy as np
   from sandlercubics.eos import PengRobinsonEOS
   from sandlerprops.properties import PropertiesDatabase
   
   db = PropertiesDatabase()
   nitrogen = db.get_compound('nitrogen')
   
   eos = PengRobinsonEOS(
       Tc=nitrogen.Tc,
       Pc=nitrogen.Pc / 10,
       omega=nitrogen.Omega
   )
   
   # Isothermal process
   T = 300  # K (constant)
   P_initial = 1.0  # MPa
   P_final = 10.0  # MPa
   
   # Initial state
   eos.T = T
   eos.P = P_initial
   Z1, v1 = eos.Z, eos.v.item()
   H1, S1 = eos.Hdep, eos.Sdep
   
   # Final state
   eos.P = P_final
   Z2, v2 = eos.Z, eos.v.item()
   H2, S2 = eos.Hdep, eos.Sdep
   
   # Property changes (ideal gas contribution + departure function changes)
   R = 8.314  # J/(mol·K)
   Delta_H_dep = H2 - H1
   Delta_S_dep = S2 - S1
   Delta_S_ig = -R * np.log(P_final / P_initial)
   
   print(f"Isothermal Compression of Nitrogen")
   print(f"T = {T} K (constant)")
   print(f"P: {P_initial} → {P_final} MPa\n")
   print(f"State Changes:")
   print(f"  Z: {Z1:.4f} → {Z2:.4f}")
   print(f"  v: {v1:.6f} → {v2:.6f} m³/mol\n")
   print(f"Property Changes:")
   print(f"  ΔH_dep = {Delta_H_dep:.2f} J/mol")
   print(f"  ΔS_dep = {Delta_S_dep:.4f} J/(mol·K)")
   print(f"  ΔS_ig  = {Delta_S_ig:.4f} J/(mol·K)")
   print(f"  ΔS_total = {Delta_S_dep + Delta_S_ig:.4f} J/(mol·K)")

Example 6: Isobaric Heating
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate enthalpy change during constant pressure heating:

.. code-block:: python

   from sandlercubics.eos import PengRobinsonEOS
   from sandlerprops.properties import PropertiesDatabase
   import numpy as np
   
   db = PropertiesDatabase()
   methane = db.get_compound('methane')
   
   eos = PengRobinsonEOS(
       Tc=methane.Tc,
       Pc=methane.Pc / 10,
       omega=methane.Omega
   )
   
   # Process conditions
   P = 5.0  # MPa (constant)
   T1, T2 = 300, 400  # K
   
   # Heat capacity coefficients (from database or literature)
   # Cp_ig = A + B*T + C*T^2 + D*T^3
   CpA = 19.25  # J/(mol·K)
   CpB = 5.213e-2  # J/(mol·K²)
   CpC = 1.197e-5  # J/(mol·K³)
   CpD = -1.132e-8  # J/(mol·K⁴)
   
   # State 1
   eos.T = T1
   eos.P = P
   Hdep1 = eos.Hdep
   
   # State 2
   eos.T = T2
   eos.P = P
   Hdep2 = eos.Hdep
   
   # Ideal gas enthalpy change (integrate Cp)
   def cp_integral(T):
       return CpA * T + CpB * T**2 / 2 + CpC * T**3 / 3 + CpD * T**4 / 4
   
   Delta_H_ig = cp_integral(T2) - cp_integral(T1)
   Delta_H_dep = Hdep2 - Hdep1
   Delta_H_total = Delta_H_ig + Delta_H_dep
   
   print(f"Isobaric Heating of Methane")
   print(f"P = {P} MPa (constant)")
   print(f"T: {T1} → {T2} K\n")
   print(f"Enthalpy Changes:")
   print(f"  ΔH_ig  = {Delta_H_ig:.2f} J/mol")
   print(f"  ΔH_dep = {Delta_H_dep:.2f} J/mol")
   print(f"  ΔH_total = {Delta_H_total:.2f} J/mol")

