.. _theory:

Theoretical Background
======================

This section provides the theoretical foundation for cubic equations of state and their implementation in sandlercubics.

Cubic Equations of State
-------------------------

Cubic equations of state are mathematical models that relate pressure, temperature, and molar volume of a pure substance or mixture. They are called "cubic" because they result in a cubic equation when solved for molar volume.

The general form can be written as:

.. math::

   P = \frac{RT}{v-b} - \frac{a(T)}{v^2 + ubv + wb^2}

where:

* :math:`P` = pressure
* :math:`T` = temperature
* :math:`v` = molar volume
* :math:`R` = universal gas constant (8.314 J/(mol·K))
* :math:`a(T)` = attraction parameter (temperature-dependent)
* :math:`b` = covolume (repulsion parameter)
* :math:`u, w` = equation-specific constants

Van der Waals Equation
-----------------------

The van der Waals equation (1873) is the original and simplest cubic equation of state:

.. math::

   P = \frac{RT}{v-b} - \frac{a}{v^2}

Parameters
~~~~~~~~~~

.. math::

   a = \frac{27R^2T_c^2}{64P_c}

.. math::

   b = \frac{RT_c}{8P_c}

where :math:`T_c` and :math:`P_c` are the critical temperature and pressure.

Characteristics
~~~~~~~~~~~~~~~

* Simple and intuitive: :math:`b` accounts for molecular volume, :math:`a` for attractive forces
* Historically important but less accurate than modern equations
* Good for educational purposes and understanding basic concepts
* Tends to overestimate pressures at high densities

Peng-Robinson Equation
-----------------------

The Peng-Robinson equation (1976) is widely used in the petroleum and gas industries:

.. math::

   P = \frac{RT}{v-b} - \frac{a\alpha(T)}{v(v+b) + b(v-b)}

This can be rewritten as:

.. math::

   P = \frac{RT}{v-b} - \frac{a\alpha(T)}{v^2 + 2bv - b^2}

Parameters
~~~~~~~~~~

.. math::

   a = 0.45724\frac{R^2T_c^2}{P_c}

.. math::

   b = 0.07780\frac{RT_c}{P_c}

.. math::

   \alpha(T) = \left[1 + \kappa\left(1 - \sqrt{T_r}\right)\right]^2

.. math::

   \kappa = 0.37464 + 1.54226\omega - 0.26992\omega^2

where:

* :math:`T_r = T/T_c` is the reduced temperature
* :math:`\omega` is the acentric factor

Characteristics
~~~~~~~~~~~~~~~

* Excellent accuracy for liquid densities
* Good performance near critical point
* Standard in petroleum industry
* Reliable for vapor-liquid equilibrium calculations
* Performs well for nonpolar and slightly polar substances

Soave-Redlich-Kwong Equation
-----------------------------

The Soave-Redlich-Kwong (SRK) equation (1972) modifies the original Redlich-Kwong equation:

.. math::

   P = \frac{RT}{v-b} - \frac{a\alpha(T)}{v(v+b)}

This corresponds to :math:`u=1, w=0` in the general form.

Parameters
~~~~~~~~~~

.. math::

   a = 0.42748\frac{R^2T_c^2}{P_c}

.. math::

   b = 0.08664\frac{RT_c}{P_c}

.. math::

   \alpha(T) = \left[1 + m\left(1 - \sqrt{T_r}\right)\right]^2

.. math::

   m = 0.480 + 1.574\omega - 0.176\omega^2

Characteristics
~~~~~~~~~~~~~~~

* Widely used in gas processing
* Better for vapor pressures than liquid densities
* Simpler mathematical form than Peng-Robinson
* Popular in process simulators
* Good for hydrocarbons and light gases

Compressibility Factor
----------------------

The compressibility factor :math:`Z` measures the deviation from ideal gas behavior:

.. math::

   Z = \frac{Pv}{RT}

For an ideal gas, :math:`Z = 1`. Real gases have :math:`Z \neq 1`, with:

* :math:`Z < 1`: attractive forces dominate (typical at moderate pressures)
* :math:`Z > 1`: repulsive forces dominate (typical at very high pressures)

Each cubic equation can be rewritten in terms of :math:`Z`:

.. math::

   Z^3 + pZ^2 + qZ + r = 0

where :math:`p, q, r` are equation-specific functions of temperature and pressure. The three roots correspond to:

1. Vapor phase (largest Z)
2. Liquid phase (smallest positive Z)  
3. Unphysical root (may be negative or complex)

Departure Functions
-------------------

Departure functions quantify the difference between real and ideal gas properties.

Enthalpy Departure
~~~~~~~~~~~~~~~~~~

.. math::

   H^{dep} = H^{real} - H^{ig} = \int_\infty^v \left[T\left(\frac{\partial P}{\partial T}\right)_v - P\right] dv

This integral can be evaluated analytically for each cubic equation.

Entropy Departure
~~~~~~~~~~~~~~~~~

.. math::

   S^{dep} = S^{real} - S^{ig} = \int_\infty^v \left[\left(\frac{\partial P}{\partial T}\right)_v - \frac{R}{v}\right] dv

Phase Equilibrium
-----------------

At vapor-liquid equilibrium, the cubic equation has three real roots. The equilibrium condition is:

.. math::

   f^V = f^L

where :math:`f` is the fugacity. For cubic equations:

.. math::

   \ln \phi = \frac{1}{RT}\int_\infty^v \left[P - \frac{RT}{v}\right] dv - \ln Z

where :math:`\phi = f/P` is the fugacity coefficient.

The Maxwell equal-area rule provides an alternative equilibrium criterion: the area between the isotherm and the equilibrium pressure must be equal for the liquid and vapor regions.

Implementation Notes
--------------------

Numerical Methods
~~~~~~~~~~~~~~~~~

sandlercubics uses several numerical techniques:

1. **Root finding**: Analytical solution of cubic equation or iterative methods
2. **Phase identification**: Gibbs energy minimization to select physical root
3. **Saturation calculations**: Iterative solution of fugacity equality
4. **Property integration**: Analytical integration where possible

Convergence
~~~~~~~~~~~

Near the critical point, cubic equations become numerically challenging:

* Multiple roots converge
* Small changes in T or P cause large changes in properties
* Derivatives become very large

The package includes safeguards for these conditions but may still encounter convergence issues very close to the critical point.

Limitations
-----------

Users should be aware of fundamental limitations:

1. **Polar compounds**: Cubic equations are less accurate for highly polar substances (water, alcohols, acids)
2. **Associating fluids**: Hydrogen bonding effects are not captured
3. **Critical region**: Accuracy decreases near the critical point
4. **Quantum effects**: Not applicable to cryogenic helium or hydrogen at very low temperatures
5. **Mixtures**: Current implementation is for pure substances only

Accuracy Expectations
---------------------

Typical accuracy for well-characterized compounds:

==================  ==================  ==================
Property            Peng-Robinson       Soave-Redlich-Kwong
==================  ==================  ==================
Vapor pressure      ±5-10%              ±5-15%
Liquid density      ±2-5%               ±5-10%
Vapor density       ±1-3%               ±1-3%
Critical region     ±10-30%             ±10-30%
==================  ==================  ==================

Errors are generally smaller for hydrocarbons and increase for polar compounds.

References
----------

1. Sandler, S. I. (2017). *Chemical, Biochemical, and Engineering Thermodynamics* (5th ed.). Wiley.

2. van der Waals, J. D. (1873). "Over de Continuiteit van den Gas- en Vloeistoftoestand". PhD thesis, Leiden University.

3. Peng, D. Y., & Robinson, D. B. (1976). "A New Two-Constant Equation of State". *Industrial & Engineering Chemistry Fundamentals*, 15(1), 59-64.

4. Soave, G. (1972). "Equilibrium constants from a modified Redlich-Kwong equation of state". *Chemical Engineering Science*, 27(6), 1197-1203.

5. Redlich, O., & Kwong, J. N. S. (1949). "On the Thermodynamics of Solutions. V. An Equation of State. Fugacities of Gaseous Solutions". *Chemical Reviews*, 44(1), 233-244.

Further Reading
---------------

For deeper understanding of thermodynamic theory and cubic equations:

* Reid, R. C., Prausnitz, J. M., & Poling, B. E. (1987). *The Properties of Gases and Liquids* (4th ed.). McGraw-Hill.
* Elliott, J. R., & Lira, C. T. (2012). *Introductory Chemical Engineering Thermodynamics* (2nd ed.). Prentice Hall.
* Poling, B. E., Prausnitz, J. M., & O'Connell, J. P. (2001). *The Properties of Gases and Liquids* (5th ed.). McGraw-Hill.
