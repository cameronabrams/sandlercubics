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

   from sandlercubics.eos import PengRobinsonEOS
   from sandlerprops.properties import PropertiesDatabase
   
   # Load compound data
   db = PropertiesDatabase()
   propane = db.get_compound('propane')
   
   # Initialize EOS
   eos = PengRobinsonEOS(
       Tc=propane.Tc,
       Pc=propane.Pc / 10,  # Convert MPa to bar
       omega=propane.Omega
   )
   
   # Set state
   eos.T = 350  # K
   eos.P = 2.0  # MPa
   
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

Example 4: Saturation Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate properties along the saturation curve:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from sandlercubics.eos import PengRobinsonEOS
   from sandlerprops.properties import PropertiesDatabase
   
   db = PropertiesDatabase()
   methane = db.get_compound('methane')
   
   eos = PengRobinsonEOS(
       Tc=methane.Tc,
       Pc=methane.Pc / 10,
       omega=methane.Omega
   )
   
   # Temperature range (subcritical)
   T_range = np.linspace(120, 185, 20)  # K
   
   # Calculate saturation pressures (this is a placeholder)
   # You'll need to implement or use your package's saturation solver
   P_sat = []
   v_liquid = []
   v_vapor = []
   
   for T in T_range:
       # Calculate saturation pressure for this temperature
       # This is pseudocode - adjust based on your API
       # Psat, vl, vv = eos.saturation(T)
       # P_sat.append(Psat)
       # v_liquid.append(vl)
       # v_vapor.append(vv)
       pass
   
   # Plot vapor pressure curve
   # plt.figure(figsize=(10, 6))
   # plt.plot(T_range, P_sat, 'b-', linewidth=2)
   # plt.xlabel('Temperature (K)')
   # plt.ylabel('Pressure (MPa)')
   # plt.title('Vapor Pressure Curve for Methane')
   # plt.grid(True, alpha=0.3)
   # plt.show()

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

Advanced Applications
---------------------

Example 7: PV Diagram
~~~~~~~~~~~~~~~~~~~~~

Generate a pressure-volume diagram showing isotherms:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from sandlercubics.eos import PengRobinsonEOS
   from sandlerprops.properties import PropertiesDatabase
   
   db = PropertiesDatabase()
   co2 = db.get_compound('carbon dioxide')
   
   eos = PengRobinsonEOS(
       Tc=co2.Tc,
       Pc=co2.Pc / 10,
       omega=co2.Omega
   )
   
   # Create volume range
   v_range = np.logspace(-4, -1, 200)  # m³/mol
   
   # Calculate isotherms
   temperatures = [280, 300, 304.2, 320, 350]  # K (Tc ≈ 304.2 K)
   
   plt.figure(figsize=(12, 8))
   
   for T in temperatures:
       eos.T = T
       pressures = []
       
       for v in v_range:
           # Calculate pressure from volume (may need special method)
           # This is pseudocode
           # P = eos.P_from_v(v)
           # pressures.append(P)
           pass
       
       # label = f'T = {T} K'
       # if abs(T - 304.2) < 1:
       #     label += ' (critical)'
       # plt.plot(v_range, pressures, linewidth=2, label=label)
   
   plt.xlabel('Molar Volume (m³/mol)')
   plt.ylabel('Pressure (MPa)')
   plt.title('PV Diagram for CO₂ (Peng-Robinson EOS)')
   plt.xscale('log')
   plt.yscale('log')
   plt.grid(True, alpha=0.3, which='both')
   plt.legend()
   plt.tight_layout()
   # plt.show()

Example 8: Reduced Property Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Plot compressibility factor vs reduced pressure:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   from sandlercubics.eos import PengRobinsonEOS
   
   # Use nitrogen properties
   Tc, Pc, omega = 126.2, 3.39, 0.037  # K, MPa, dimensionless
   
   eos = PengRobinsonEOS(Tc, Pc / 10, omega)
   
   # Reduced conditions
   T_r_values = [1.0, 1.2, 1.5, 2.0]  # Tr = T/Tc
   P_r_range = np.linspace(0.1, 10, 100)  # Pr = P/Pc
   
   plt.figure(figsize=(12, 8))
   
   for T_r in T_r_values:
       T = T_r * Tc
       eos.T = T
       Z_values = []
       
       for P_r in P_r_range:
           P = P_r * Pc
           eos.P = P
           Z_values.append(eos.Z)
       
       plt.plot(P_r_range, Z_values, linewidth=2, label=f'Tr = {T_r}')
   
   # Ideal gas line
   plt.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='Ideal gas')
   
   plt.xlabel('Reduced Pressure (Pr)')
   plt.ylabel('Compressibility Factor (Z)')
   plt.title('Generalized Compressibility Chart (Peng-Robinson)')
   plt.grid(True, alpha=0.3)
   plt.legend()
   plt.xlim([0, 10])
   plt.ylim([0, 2])
   plt.tight_layout()
   # plt.show()

Example 9: Custom Compound
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Work with a compound not in the database:

.. code-block:: python

   from sandlercubics.eos import PengRobinsonEOS
   
   # Custom compound properties (e.g., from literature)
   # Example: n-hexane
   Tc = 507.6  # K
   Pc = 30.25  # bar (not MPa for EOS initialization!)
   omega = 0.301
   
   eos = PengRobinsonEOS(Tc, Pc, omega)
   
   # Calculate properties
   eos.T = 400  # K
   eos.P = 2.0  # MPa
   
   print(f"Custom Compound Properties:")
   print(f"  Critical temperature: {Tc} K")
   print(f"  Critical pressure: {Pc} bar")
   print(f"  Acentric factor: {omega}")
   print(f"\nState at T = {eos.T} K, P = {eos.P} MPa:")
   print(f"  Z = {eos.Z:.4f}")
   print(f"  v = {eos.v.item():.6f} m³/mol")

Example 10: Batch Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Process multiple compounds efficiently:

.. code-block:: python

   from sandlercubics.eos import PengRobinsonEOS
   from sandlerprops.properties import PropertiesDatabase
   import pandas as pd
   
   db = PropertiesDatabase()
   
   # Compounds to analyze
   compounds = ['methane', 'ethane', 'propane', 'n-butane']
   
   # Fixed conditions
   T, P = 350, 5.0  # K, MPa
   
   results = []
   
   for name in compounds:
       compound = db.get_compound(name)
       eos = PengRobinsonEOS(compound.Tc, compound.Pc / 10, compound.Omega)
       eos.T = T
       eos.P = P
       
       results.append({
           'Compound': name,
           'Tc (K)': compound.Tc,
           'Pc (MPa)': compound.Pc,
           'ω': compound.Omega,
           'Z': eos.Z,
           'v (m³/mol)': eos.v.item(),
           'Hdep (J/mol)': eos.Hdep,
           'Sdep (J/(mol·K))': eos.Sdep
       })
   
   # Create DataFrame for nice display
   df = pd.DataFrame(results)
   print(f"\nResults at T = {T} K, P = {P} MPa:\n")
   print(df.to_string(index=False))

Notes on Examples
-----------------

**Placeholder Code**

Some examples include placeholder code (marked with comments) for functionality that depends on specific implementation details in your package. You'll need to:

1. Replace pseudocode with actual API calls
2. Implement any missing methods (e.g., saturation calculations)
3. Adjust property access patterns to match your class structure

**Unit Consistency**

Always check units carefully:

* EOS initialization typically uses **bar** for pressure
* CLI and most calculations use **MPa**
* Temperature is always in **Kelvin**
* Property outputs are in SI units

**Performance Considerations**

For large-scale calculations:

* Reuse EOS objects when possible
* Batch similar calculations together
* Consider vectorizing operations if your implementation supports it
