# Author: Cameron F. Abrams, <cfa22@drexel.edu>

from .eos import CubicEOS, sqrt_2
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger(__name__)

@dataclass
class SoaveRedlichKwongEOS(CubicEOS):
    """
    Pure-component Soave-Redlich-Kwong equation of state (Eq. 6.4-1 in Sandler 5th ed)
    """

    description: str = "Soave-Redlich-Kwong Equation of State"

    @property
    def kappa(self):
        """
        kappa parameter for Soave-Redlich-Kwong EOS
        """
        return 0.480 + 1.574 * self.omega - 0.176 * self.omega**2
    
    @property
    def alpha(self):
        """
        alpha parameter for Soave-Redlich-Kwong EOS
        """
        return (1 + self.kappa * (1 - np.sqrt(self.T/self.Tc)))**2

    def _calc_a(self):
        """
        Calculates parameter a for Soave-Redlich-Kwong EOS
        """
        return 0.42748 * self.R**2 * self.Tc**2 / self.Pc * self.alpha
    
    def _calc_b(self):
        """
        Calculates parameter b for Soave-Redlich-Kwong EOS
        """
        return 0.08664 * self.R * self.Tc / self.Pc

    @property
    def cubic_coeff(self):
        """
        cubic coefficients for Soave-Redlich-Kwong EOS
        """
        A, B = self.A, self.B
        return np.array([1.0, -1.0, A - B - B**2, -A * B])

    @property
    def lrfrac(self) -> np.ndarray:
        """
        helper function for ln fugacity, enthalpy, and entropy departures
        """
        z = self.Z
        num_arg = z + self.B
        den_arg = z
        return np.log(num_arg / den_arg)
    
    def _calc_h_departure(self) -> np.ndarray:
        """
        Enthalpy departure at state T and P (from solution to problem 6.36 in Sandler 5th ed)
        """
        z = self.Z
        return self.R * self.T * (z - 1) + (self.T * self.da_dT - self.a) / self.b * self.lrfrac

    def _calc_s_departure(self) -> np.ndarray:
        """
        Entropy departure at state T and P (from solution to problem 6.36 in Sandler 5th ed)
        """
        z = self.Z
        return self.R * np.log(z - self.B) + self.da_dT/self.b * self.lrfrac

    def _calc_logphi(self) -> np.ndarray:
        """
        natural log of fugacity coefficient at state T and P (from solution to problem 6.46 in Sandler 5th ed)
        """
        z = self.Z
        return z - 1 - np.log(z - self.B) - self.a / (self.R * self.T * self.b) * self.lrfrac