# Author: Cameron F. Abrams, <cfa22@drexel.edu>
from .eos import CubicEOS
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger(__name__)

@dataclass
class VanDerWaalsEOS(CubicEOS):
    """
    Generalized van der Waals equation of state, in which 'a' and 'b'
    are calculated from critical constants
    """

    description: str = "Van der Waals Equation of State"

    def _calc_a(self):
        """
        Calculates parameter a for Van der Waals EOS
        """
        return (27 / 64) * self.R**2 * self.Tc**2 / self.Pc
    
    def _calc_b(self):
        """
        Calculates parameter b for Van der Waals EOS
        """
        return self.R * self.Tc / (8 * self.Pc)

    @property
    def cubic_coeff(self):
        """
        cubic coefficients for Van der Waals EOS
        """
        A, B = self.A, self.B
        return np.array([1.0, -1.0 - B, A, -A * B])

    def _calc_h_departure(self) -> np.ndarray:
        """
        Enthalpy departure for Van der Waals EOS
        """
        return self.R * self.T * (self.Z - 1 - self.A * np.reciprocal(self.Z))

    def _calc_s_departure(self) -> np.ndarray:
        """
        Entropy departure for Van der Waals EOS
        """
        return self.R * (np.log(self.Z - self.B) - self.A / self.Z)

    def _calc_logphi(self) -> np.ndarray:
        """
        natural log of fugacity coefficient for Van der Waals EOS
        """
        return self.Z - 1 - np.log(self.Z - self.B) - self.A / self.Z