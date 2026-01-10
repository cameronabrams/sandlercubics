# Author: Cameron F. Abrams, <cfa22@drexel.edu>

from .eos import CubicEOS, sqrt_2
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger(__name__)

@dataclass
class PengRobinsonEOS(CubicEOS):
    """
    Pure-component Peng-Robinson equation of state
    """

    description: str = "Peng-Robinson Equation of State"

    @property
    def kappa(self):
        """
        kappa parameter for Peng-Robinson EOS
        """
        return 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
    
    @property
    def alpha(self):
        """
        alpha parameter for Peng-Robinson EOS
        """
        return (1 + self.kappa * (1 - np.sqrt(self.T / self.Tc)))**2

    def _calc_a(self):
        """
        Calculates parameter a for Peng-Robinson EOS
        """
        return 0.45724 * self.R**2 * self.Tc**2 / self.Pc * self.alpha
    
    def _calc_b(self):
        """
        Calculates parameter b for Peng-Robinson EOS
        """
        return 0.07780 * self.R * self.Tc / self.Pc

    def _calc_da_dT(self):
        """
        Temperature derivative of parameter a for Peng-Robinson EOS
        """
        return -self.a * self.kappa / np.sqrt(self.alpha * self.T * self.Tc)
    
    @property
    def cubic_coeff(self):
        """
        cubic coefficients for Peng-Robinson EOS
        """
        A, B = self.A, self.B
        return np.array([1.0, -1 + B, A - 3 * B**2 - 2*B, -A * B + B**2 + B**3])

    @property
    def lrfrac(self) -> np.ndarray:
        """
        Helper property for logarithmic fraction term in Peng-Robinson EOS
        """
        z = self.Z
        B = self.B
        num_arg = z + (1 + sqrt_2) * B
        den_arg = z + (1 - sqrt_2) * B
        return np.log(num_arg / den_arg)
        
    def _calc_h_departure(self) -> np.ndarray:
        """
        Enthalpy departure for Peng-Robinson EOS
        """
        z = self.Z
        return self.R * self.T * (z - 1) + (self.T * self.da_dT - self.a) / (2 * sqrt_2 * self.b) * self.lrfrac

    def _calc_s_departure(self) -> np.ndarray:
        """
        Entropy departure for Peng-Robinson EOS
        """
        z = self.Z
        return self.R * np.log(z - self.B) + self.da_dT/(2 * sqrt_2 * self.b) * self.lrfrac
    
    def _calc_logphi(self) -> np.ndarray:
        """
        natural log of fugacity coefficient for Peng-Robinson EOS
        """
        z = self.Z
        return z - 1 - np.log(z - self.B) - self.A / (2 * np.sqrt(2) * self.B) * self.lrfrac

