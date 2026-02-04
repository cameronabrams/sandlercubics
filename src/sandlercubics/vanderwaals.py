# Author: Cameron F. Abrams, <cfa22@drexel.edu>
from .eos import CubicEOS
from dataclasses import dataclass
import numpy as np
import logging
from sandlermisc import ureg, R

logger = logging.getLogger(__name__)

@dataclass
class VanDerWaalsEOS(CubicEOS):
    """
    Generalized van der Waals equation of state, in which 'a' and 'b'
    are calculated from critical constants
    """
    
    name: str = "Van der Waals Equation of State"
    description: str = "Van der Waals Equation of State"

    def _calc_a(self):
        """
        Calculates parameter a for Van der Waals EOS
        """
        return (27 / 64) * R**2 * self.Tc**2 / self.Pc
    
    def _calc_da_dT(self) -> float:
        return 0.0

    def _calc_b(self):
        """
        Calculates parameter b for Van der Waals EOS
        """
        return R * self.Tc / (8 * self.Pc)

    def _calc_P(self):
        """
        Calculates pressure from Van der Waals EOS
        """
        v = self.v
        a = self.a
        b = self.b
        T = self.T
        return R * T / (v - b) - a / v**2

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
        return R * self.T * (self.Z - 1 - self.A * np.reciprocal(self.Z))

    def _calc_s_departure(self) -> np.ndarray:
        """
        Entropy departure for Van der Waals EOS
        """
        return R * (np.log(self.Z - self.B) - self.A / self.Z)

    def _calc_logphi(self) -> np.ndarray:
        """
        natural log of fugacity coefficient for Van der Waals EOS
        """
        return self.Z - 1 - np.log(self.Z - self.B) - self.A / self.Z