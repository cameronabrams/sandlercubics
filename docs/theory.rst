.. _theory:

Theoretical Background
======================

This section provides the theoretical foundation for cubic equations of state and their implementation in sandlercubics.  The presentation follows closely the treatment found in *Chemical, Biochemical, and Engineering Thermodynamics* by Stan Sandler [1]_.

Cubic Equations of State
-------------------------

Cubic equations of state are mathematical models that relate pressure, temperature, and molar volume of a pure substance or mixture. They are called "cubic" because they result in a cubic equation expressed in terms of compressibility factor :math:`Z`:

.. math::

   Z \equiv \frac{Pv}{RT}

Each cubic equation can be rewritten in terms of :math:`Z`:

.. math::

   Z^3 + pZ^2 + qZ + r = 0

where :math:`p, q, r` are equation-specific functions of temperature and pressure. The three possible classes of solutions
to the cubic equation depend on the thermodynamic state:

1. Single real root (supercritical or single-phase states)
2. Three real roots (two-phase region below critical point)

When three real roots exist, they correspond to:

1. Vapor phase (largest Z)
2. Liquid phase (smallest positive Z)  
3. Non-physical root (intermediate Z, discarded)


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

where :math:`T_c` and :math:`P_c` are the critical temperature and pressure.  When cast as a cubic in :math:`Z`, the van der Waals equation becomes:

.. math::

   Z^3 - (1+B)Z^2 + AZ - AB = 0

where:

.. math::

   A = \frac{aP}{R^2T^2}

and

.. math::

   B = \frac{bP}{RT}

By convention, all cubics use the same formulations of :math:`A` and :math:`B` for consistency.

Peng-Robinson Equation
-----------------------

The Peng-Robinson equation (1976) is widely used in the petroleum and gas industries:

.. math::

   P = \frac{RT}{v-b} - \frac{a\alpha(T)}{v(v+b) + b(v-b)}

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

When expressed as a cubic in :math:`Z`, the Peng-Robinson equation becomes:

.. math::

   Z^3 - (1-B)Z^2 + (A - 3B^2 - 2B)Z - (AB - B^2 - B^3) = 0

Soave-Redlich-Kwong Equation
-----------------------------

The Soave-Redlich-Kwong (SRK) equation (1972) modifies the original Redlich-Kwong equation:

.. math::

   P = \frac{RT}{v-b} - \frac{a\alpha(T)}{v(v+b)}

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

When cast as a cubic in :math:`Z`, the SRK equation becomes:

.. math::

   Z^3 - (1-B)Z^2 + (A - B - B^2)Z - AB = 0


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

Implementation Notes
--------------------

Numerical Methods
~~~~~~~~~~~~~~~~~

sandlercubics uses two main numerical techniques:

1. **Root finding**: Analytical solution of cubic equation or iterative methods
2. **Saturation calculations**: Iterative solution of fugacity equality as explained in Chapter 7 of Sandler [1]_

Convergence
~~~~~~~~~~~

Near the critical point, cubic equations become numerically challenging:

* Multiple roots converge
* Small changes in T or P cause large changes in properties
* Derivatives become very large

Limitations
-----------

Users should be aware of fundamental limitations:

1. **Polar compounds**: Cubic equations are less accurate for highly polar substances (water, alcohols, acids)
2. **Associating fluids**: Hydrogen bonding effects are not captured
3. **Critical region**: Accuracy decreases near the critical point
4. **Quantum effects**: Not applicable to cryogenic helium or hydrogen at very low temperatures
5. **Mixtures**: Current implementation is for pure substances only

References
----------
.. _1:

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
