# Author: Cameron F. Abrams, <cfa22@drexel.edu>
from .eos import CubicEOS
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger(__name__)

@dataclass
class IdealGasEOS(CubicEOS):
    """
    Ideal gas class derived from CubicEOS
    """

    description: str = "Ideal Gas Equation of State"

    def _calc_a(self):
        """ vdW a parameter for ideal gas is zero """
        return 0.0
    
    def _calc_b(self):
        """ vdW b parameter for ideal gas is zero """
        return 0.0
    
    @property
    def cubic_coeff(self) -> np.ndarray:
        """ cubic coefficients for ideal gas """
        return np.array([0.0, 0.0, 1.0, -1.0])

    @property
    def Z(self) -> np.ndarray:
        """ compressibility factor for ideal gas is unity """
        return np.array([1.0])
    
    def _calc_h_departure(self) -> np.ndarray:
        """ Enthalpy departure at state T and P; ideal-gas value """
        return np.array([0.0])

    def _calc_s_departure(self) -> np.ndarray:
        """ Entropy departure at state T and P; ideal-gas value """
        return np.array([0.0])
    
    def _calc_logphi(self) -> np.ndarray:
        """ natural log of fugacity coefficient at state T and P; ideal-gas value """
        return np.array([0.0])

    def solve(self, T: float = None, P: float = None, v: float = None):
        """
        Solve for missing property in ideal gas EOS

        Parameters
        ----------
        T : float, optional
            Temperature
        P : float, optional
            Pressure
        v : float, optional
            Molar volume
        """
        if T is not None:
            self.T = T
        if P is not None:
            self.P = P
        R = self.R_pv # volumetric gas constant
        if self.T is not None and self.P is not None:
            pass
        elif self.T is not None and v is not None:
            self.P = R * self.T / v
        elif self.P is not None and v is not None:
            self.T = self.P * v / R
        else:
            raise ValueError("Insufficient information to solve ideal gas EOS.")
