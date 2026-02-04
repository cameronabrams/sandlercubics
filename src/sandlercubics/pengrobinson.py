# Author: Cameron F. Abrams, <cfa22@drexel.edu>

from __future__ import annotations
from .eos import CubicEOS, sqrt_2
from dataclasses import dataclass
import numpy as np
import logging
from sandlermisc import ureg, R

logger = logging.getLogger(__name__)

@dataclass
class PengRobinsonEOS(CubicEOS):
    """
    Pure-component Peng-Robinson equation of state
    """

    name: str = "Peng-Robinson Equation of State"
    description: str = "Peng-Robinson Equation of State"
    
    def _calc_kappa(self) -> float:
        """
        kappa parameter for Peng-Robinson EOS
        """
        return 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2

    @property
    def kappa(self) -> float:
        """
        kappa parameter for Peng-Robinson EOS
        """
        return self._calc_kappa()
    
    def _calc_alpha(self) -> float:
        """
        alpha parameter for Peng-Robinson EOS
        """
        return (1 + self.kappa * (1 - np.sqrt(self.T / self.Tc)))**2

    @property
    def alpha(self) -> float:
        """
        alpha parameter for Peng-Robinson EOS
        """
        return self._calc_alpha()

    def _calc_a(self) -> float:
        """
        Calculates parameter a for Peng-Robinson EOS
        """
        return 0.45724 * R**2 * self.Tc**2 / self.Pc * self.alpha
    
    def _calc_b(self) -> float:
        """
        Calculates parameter b for Peng-Robinson EOS
        """
        return 0.07780 * R * self.Tc / self.Pc

    def _calc_da_dT(self) -> float:
        """
        Temperature derivative of parameter a for Peng-Robinson EOS
        """
        return -self.a * self.kappa / np.sqrt(self.alpha * self.T * self.Tc)
    
    def _calc_P(self) -> float:
        """
        Calculates pressure from Peng-Robinson EOS
        """
        v = self.v
        a = self.a
        b = self.b
        T = self.T
        return R * T / (v - b) - a / (v**2 + 2 * b * v - b**2)

    @property
    def cubic_coeff(self) -> np.ndarray:
        """
        cubic coefficients for Peng-Robinson EOS
        """
        A, B = self.A, self.B
        return np.array([1.0, -1 + B, A - 3 * B**2 - 2*B, -A * B + B**2 + B**3])

    @property
    def lrfrac(self) -> float:
        """
        Helper property for logarithmic fraction term in Peng-Robinson EOS
        """
        z = self.Z
        B = self.B
        num_arg = z + (1 + sqrt_2) * B
        den_arg = z + (1 - sqrt_2) * B
        if num_arg * den_arg <= 0:
            logger.warning(f'lrfrac: Invalid argument for logarithm: num_arg={num_arg}, den_arg={den_arg}. Returning NaN.')
            logger.warning(f'lrfrac: z={z}, B={B} 1+sqrt(2)={1 + sqrt_2}, 1 - sqrt(2)={1 - sqrt_2}')
            raise ValueError('Invalid argument for logarithm in lrfrac calculation.')
        return np.log(num_arg / den_arg)
        
    def _calc_h_departure(self) -> float:
        """
        Enthalpy departure for Peng-Robinson EOS
        """
        z = self.Z
        return R * self.T * (z - 1) + (self.T * self.da_dT - self.a) / (2 * sqrt_2 * self.b) * self.lrfrac

    def _calc_s_departure(self) -> float:
        """
        Entropy departure for Peng-Robinson EOS
        """
        z = self.Z
        return R * np.log(z - self.B) + self.da_dT/(2 * sqrt_2 * self.b) * self.lrfrac
    
    def _calc_logphi(self) -> float:
        """
        natural log of fugacity coefficient for Peng-Robinson EOS
        """
        z = self.Z
        return z - 1 - np.log(z - self.B) - self.A / (2 * np.sqrt(2) * self.B) * self.lrfrac

    