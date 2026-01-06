# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Pure-component cubic equations of state and their state calculations
"""
from __future__ import annotations
import logging
import numpy as np

from abc import ABC, abstractmethod
from dataclasses import dataclass

from sandlermisc.gas_constant import GasConstant
from sandlermisc.thermals import DeltaH_IG, DeltaS_IG

logger = logging.getLogger(__name__)

sqrt_2 = np.sqrt(2)

@dataclass
class CubicEOS(ABC):
    """
    Abstract class for all Cubic equations of state.
    """

    pressure_unit: str = 'mpa' # MPa
    """ units of pressure, lower-case """
    volume_unit: str = 'm3'
    """ units of volume, lower-case """
    temperature_unit: str = 'K'
    """ units of temperature, not used """

    phase: str = 'unspecified' 
    """ strict phase requirement flag; 'vapor', 'liquid', 'unspecified' (multiple roots) """

    P: float = 0.1
    """ state pressure """
    T: float = 298.15
    """ state temperature """
    Pc: float = 0.1
    """ critical pressure """
    Tc: float = 298.15
    """ critical temperature """
    omega: float = 0.0 
    """ acentricity factor, dimensionless """

    logiter: bool = False
    """ flag for logging iterations in phase calculations """
    maxiter: int = 100
    """ maximum iterations in phase calculations """
    epsilon: float = 1.e-5
    """ fugacity tolerance in phase calculations """
    iter: int = 0
    """ iteration counter in phase calculations """
    err: float = 0.0
    """ current error in phase calculations """

    @property
    def Pv(self):
        """ product of pressure and molar volume """
        return self.P * self.v
    
    @property
    def R_pv(self):
        return GasConstant(self.pressure_unit, self.volume_unit)
    
    @property
    def R(self):
        """ default gas constant value """
        return GasConstant() # J/(mol-K)

    @property
    @abstractmethod
    def a(self):
        """ van der Waals 'a' parameter """
        pass
    
    @property
    def da_dT(self):
        """ derivative of vdW 'a' parameter wrt temperature """
        return 0.0

    @property
    @abstractmethod
    def b(self):
        """ van der Waals 'b' parameter """
        pass
    
    @property
    def A(self):
        """ dimensionless vdW 'a' parameter """
        return self.a * self.P / (self.R * self.T)**2
    
    @property
    def B(self):
        """ dimensionelss vdW 'b' parameter """
        return self.b * self.P / (self.R * self.T)
    
    @property
    def cubic_coeff(self): # default for vdw eos
        """ coefficients of cubic equation of state in Z in descending order """ 
        return np.array([1, -1 - self.B, self.A, -self.A * self.B])
    
    @property
    def Z(self):
        """ Compressibility factor from solution of cubic.  Phase requirement enforced if set. """
        complx_roots = np.roots(self.cubic_coeff)
        real_roots_idx = np.where(complx_roots.imag==0)[0]
        real_roots = complx_roots[real_roots_idx].real
        if len(real_roots) == 1:
            return real_roots[0]
        else:
            if self.phase == 'vapor':
                return np.max(real_roots)
            elif self.phase == 'liquid':
                return np.min(real_roots)
            else: # unspecified/multiple
                return np.array([real_roots[0], real_roots[2]])

    @property
    def v(self):
        """ molar volume """
        return self.Z * self.R_pv * self.T / self. P

    @property
    def h_departure(self):
        """ Enthalpy departure at state T and P; ideal-gas value by default """
        return 0.0

    @property
    def s_departure(self):
        """ Entropy departure at state T and P; ideal-gas value by default """
        return 0.0
    
    @property
    def logphi(self):
        """ natural log of fugacity coefficient at state T and P; ideal-gas value by default """
        return 0.0
    
    @property
    def phi(self): # fugacity coefficient(s)
        logPhi = self.logphi
        return np.exp(logPhi)

    @property
    def f(self): # fugacity/fugacities
        return self.phi * self.P

    @property
    def Pvap(self):
        """ Vapor pressure at current state temperature.  Current state pressure is retained. """
        saveP = self.P
        self.P = self.Pc * (self.T / self.Tc)**8
        keepgoing = True
        self.iter = 0
        while keepgoing:
            self.iter += 1
            try:
                fV, fL = self.f
            except:
                self.P = saveP
                return None
            self.err = np.abs(fL / fV - 1)
            if self.logiter: logger.debug(f'Iter {self.iter}: P {self.P:.6f}, fV {fV:.6f}, fL {fL:.6f}; error {self.err:.4e}')
            self.P *= fL / fV
            if self.err < self.epsilon or self.iter == self.maxiter:
                keepgoing = False
            if self.iter >= self.maxiter:
                logger.warning(f'Reached {self.iter} iterations without convergence; error {np.abs(fL/fV-1):.4e}')
        Pvap = self.P
        self.P = saveP
        return Pvap
    
    @property
    def Tsat(self):
        """ Saturation temperature at state pressure.  Current state temperature is retained. """
        saveT = self.T
        self.T = self.Tc * (self.P / self.Pc)**0.125
        keepgoing = True
        self.iter = 0
        while keepgoing:
            self.iter += 1
            try:
                fV, fL = self.f
            except:
                self.T = saveT
                return None
            self.err = np.abs(fV / fL - 1)
            if self.logiter: logger.debug(f'Iter {self.iter}: T {self.T:.6f}, fV {fV:.6f}, fL {fL:.6f}; error {self.err:.4e}')
            self.T *= (fV / fL)**(1/4)
            if self.err < self.epsilon or self.iter == self.maxiter:
                keepgoing = False
            if self.iter >= self.maxiter:
                logger.warning(f'Reached {self.iter} iterations without convergence; error {np.abs(fL/fV-1):.4e}')
        Tsat = self.T
        self.T = saveT
        return Tsat

    def unit_consistency(self, other: CubicEOS):
        consistent = self.pressure_unit == other.pressure_unit and self.volume_unit == other.volume_unit and self.temperature_unit == other.temperature_unit
        if not consistent:
            raise ValueError('inconsistent units')

    def DeltaH(self, other: CubicEOS, Cp: float | list[float] | dict [str, float]):
        """
        Computes and returns enthalpy change from self to other state
        including ideal gas contribution based on Cp
        """
        self.unit_consistency(other)
        dH_ideal = DeltaH_IG(self.T, other.T, Cp)
        return other.h_departure + dH_ideal - self.h_departure

    def DeltaS(self, other: CubicEOS, Cp: float | list[float] | dict [str, float]):
        """
        Computes and returns entropy change from self to other state
        including ideal gas contribution based on Cp
        """
        self.unit_consistency(other)
        dS_ideal = DeltaS_IG(self.T, self.P, other.T, other.P, Cp, self.R)
        return other.s_departure + dS_ideal - self.s_departure
    
    def DeltaPV(self, other: CubicEOS):
        """
        Returns Delta(PV) in thermal (not PV) units 
        """
        self.unit_consistency(other)
        return (other.Pv - self.Pv) * self.R / self.R_pv
    
    def DeltaU(self, other: CubicEOS, Cp: float | list[float] | dict [str, float]):
        """
        Returns Delta(U) (internal energy)
        """
        return self.DeltaH(other, Cp) - self.DeltaPV(other)
    
@dataclass
class IdealGasEOS(CubicEOS):
    """
    Ideal gas class derived from CubicEOS
    """
    @property
    def a(self):
        """ vdW 'a' parameter for ideal gas is zero """
        return 0.0
    
    @property
    def b(self):
        """ vdW 'b' parameter for ideal gas is zero """
        return 0.0
    
    @property
    def Z(self):
        """ compressibility factor for ideal gas is unity """
        return 1.0
    
    @property
    def f(self):
        """ fugacity for ideal gas is equal to state pressure """
        return self.P

@dataclass
class GeneralizedVDWEOS(CubicEOS):
    """
    Generalized van der Waals equation of state, in which 'a' and 'b'
    are calculated from critical constants
    """

    @property
    def a(self):
        return (27 / 64) * self.R**2 * self.Tc**2 / self.Pc
    
    @property
    def b(self):
        return self.R * self.Tc / (8 * self.Pc)

    @property
    def h_departure(self):
        z = self.Z
        return self.R * self.T * (z - 1 - self.A * np.reciprocal(z))

    @property
    def s_departure(self):
        z = self.Z
        return self.R * (np.log(z - self.B) - np.log(z))
    
    @property
    def logphi(self):
        z = self.Z
        return np.log(z / (z - self.B)) - self.A / z

@dataclass
class PengRobinsonEOS(CubicEOS):
    """
    Pure-component Peng-Robinson equation of state
    """
    @property
    def kappa(self):
        return 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
    
    @property
    def alpha(self):
        return (1 + self.kappa * (1 - np.sqrt(self.T / self.Tc)))**2

    @property
    def a(self):
        return 0.45724 * self.R**2 * self.Tc**2 / self.Pc * self.alpha
    
    @property
    def b(self):
        return 0.07780 * self.R * self.Tc / self.Pc

    @property
    def da_dT(self):
        return -self.a * self.kappa / np.sqrt(self.alpha * self.T * self.Tc)
    
    @property
    def cubic_coeff(self):
        return np.array([1.0, -1 + self.B, self.A - 3 * self.B**2 - 2*self.B, -self.A * self.B + self.B**2 + self.B**3])

    @property
    def lrfrac(self):
        z = self.Z
        num_arg = z + (1 + sqrt_2) * self.B
        den_arg = z + (1 - sqrt_2) * self.B
        return np.log(num_arg / den_arg)
    
    @property
    def h_departure(self):
        z = self.Z
        return self.R * self.T * (z - 1) + (self.T * self.da_dT - self.a) / (2 * sqrt_2 * self.b) * self.lrfrac

    @property
    def s_departure(self):
        z = self.Z
        return self.R * np.log(z - self.B) + self.da_dT/(2 * sqrt_2 * self.b) * self.lrfrac
    
    @property
    def logphi(self):
        z = self.Z
        return z - 1 - np.log(z - self.B) - self.A / (2 * np.sqrt(2) * self.B) * self.lrfrac

# @dataclass
# class SoaveRedlichKwongEOS(CubicEOS):
#     """
#     Pure-component Soave-Redlich-Kwong equation of state; watting for time to derive departures and fugacity coefficients
#     """
#     @property
#     def kappa(self):
#         return 0.480 + 1.574 * self.omega - 0.176 * self.omega**2
    
#     @property
#     def alpha(self):
#         return (1 + self.kappa * (1 - np.sqrt(self.T/self.Tc)))**2

#     @property
#     def a(self):
#         return 0.42748 * self.R**2 * self.Tc**2 / self.Pc * self.alpha
    
#     @property
#     def b(self):
#         return 0.08664 * self.R * self.Tc / self.Pc

#     @property
#     def coeff(self):
#         return np.array([1.0, -1.0, self.A - self.B - self.B**2, -self.A * self.B])