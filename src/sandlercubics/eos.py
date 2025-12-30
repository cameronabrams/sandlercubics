# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Methods for pure-component cubic equations of state
# vdw and PengRobinson
from __future__ import annotations
import numpy as np
from abc import ABC, abstractmethod
from dataclasses import dataclass
from scipy.constants import R # J/mol-K, aka m3-Pa/mol-K
from sandlermisc.gas_constant import GasConstant
from sandlermisc.thermals import DeltaH_IG, DeltaS_IG

sqrt_2 = np.sqrt(2)

@dataclass
class CubicEOS(ABC):
    pressure_unit: str = 'mpa' # MPa
    volume_unit: str = 'm3'
    temperature_unit: str = 'K'

    P: float = 0.1 # MPa
    T: float = 298.15 # K
    Pc: float = 0.1 # critical pressure, MPa
    Tc: float = 298.15 # critical temperature, K
    omega: float = 0.0 # acentricity factor, dimensionless

    showiter: bool = False
    maxiter: int = 100
    epsilon: float = 1.e-5
    iter: int = 0
    err: float = 0.0

    @property
    def Pv(self): # Pv in volumetric energy units
        return self.P * self.v

    @property
    def R_pv(self):
        return GasConstant(self.pressure_unit, self.volume_unit)

    @property
    def R(self):
        return GasConstant("pa", "m3") # aka, J

    @property
    @abstractmethod
    def a(self):
        pass
    
    @property
    def da_dT(self):
        return 0.0

    @property
    @abstractmethod
    def b(self):
        pass
    
    @property
    def A(self): # dimensionless vdw a parameter
        return self.a * self.P / (self.R * self.T)**2
    
    @property
    def B(self): # dimensionelss vdw b parameter
        return self.b * self.P / (self.R * self.T)
    
    @property
    def cubic_coeff(self): # default for vdw eos
        return np.array([1, -1 - self.B, self.A, -self.A * self.B])
    
    @property
    def Z(self):
        complx_roots = np.roots(self.cubic_coeff)
        real_roots_idx = np.where(complx_roots.imag==0)[0]
        real_roots = complx_roots[real_roots_idx].real
        if len(real_roots) == 1:
            return real_roots[0]
        else:
            return np.array([real_roots[0], real_roots[2]])

    @property
    def v(self):
        return self.Z * self.R_pv * self.T / self. P

    @property
    def h_departure(self): # default for vdw eos
        z = self.Z
        return self.R * self.T * (z - 1 - self.A * np.reciprocal(z))

    @property
    def s_departure(self):
        z = self.Z
        return self.R * (np.log(z - self.B) - np.log(z))
    
    @property
    def logphi(self): # natural log of fugacity coefficient(s)
        z = self.Z
        return np.log(z / (z - self.B)) - self.A / z
    
    @property
    def phi(self): # fugacity coefficient(s)
        logPhi = self.logphi
        return np.exp(logPhi)

    @property
    def f(self): # fugacity/fugacities
        return self.phi * self.P

    @property
    def Pvap(self): # vapor pressure at temperature T; P is retained
        saveP = self.P
        self.P = self.Pc * (self.T / self.Tc)**8
        keepgoing = True
        self.iter = 0
        while keepgoing:
            self.iter += 1
            try:
                fV, fL = self.f
            except:
                raise ValueError(f'Error computing pvap at {self.T} K')
            self.err = np.abs(fL / fV - 1)
            if self.showiter: print(f'Iter {self.iter}: P {self.P:.6f}, fV {fV:.6f}, fL {fL:.6f}; error {self.err:.4e}')
            self.P *= fL / fV
            if self.err < self.epsilon or self.iter == self.maxiter:
                keepgoing = False
            if self.iter >= self.maxiter:
                print(f'Reached {self.iter} iterations without convergence; error {np.abs(fL/fV-1):.4e}')
        Pvap = self.P
        self.P = saveP
        return Pvap
    
    def unit_consistency(self, other: CubicEOS):
        consistent = self.pressure_unit == other.pressure_unit and self.volume_unit == other.volume_unit and self.temperature_unit == other.temperature_unit
        if not consistent:
            raise ValueError('inconsistent units')

    def DeltaH(self, other: CubicEOS, Cp: float | list[float] | dict [str, float]):
        self.unit_consistency(other)
        dH_ideal = DeltaH_IG(other.T, self.T, Cp)
        return other.h_departure + dH_ideal - self.h_departure

    def DeltaS(self, other: CubicEOS, Cp: dict | float = None):
        self.unit_consistency(other)
        dS_ideal = DeltaS_IG(other.T, other.P, self.T, self.P, Cp, self.R)
        return other.s_departure + dS_ideal - self.s_departure
    
    def DeltaPV(self, other: CubicEOS):
        """
        Returns Delta(PV) in thermal (not PV) units 
        """
        self.unit_consistency(other)
        return (other.Pv - self.Pv) * self.R / self.R_pv
    
    def DeltaU(self, other: CubicEOS):
        return self.DeltaH(other) - self.DeltaPV(other)
    
@dataclass
class IdealGasEOS(CubicEOS):
    
    @property
    def a(self):
        return 0.0
    
    @property
    def b(self):
        return 0.0
    
    @property
    def Z(self):
        return 1.0
    
    @property
    def f(self):
        return self.P

@dataclass
class GeneralizedVDWEOS(CubicEOS):

    @property
    def a(self):
        return (27 / 64) * self.R**2 * self.Tc**2 / self.Pc
    
    @property
    def b(self):
        return self.R * self.Tc / (8 * self.Pc)

@dataclass
class PengRobinsonEOS(CubicEOS):

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
    def cubic_coeff(self):  # coefficients for cubic form of PR eos
        return np.array([1, -1 + self.B, self.A - 3 * self.B**2 - 2*self.B, -self.A * self.B + self.B**2 + self.B**3])

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
    def logphi(self): # natural log of fugacity coefficient
        z = self.Z
        return z - 1 - np.log(z - self.B) - self.A / (2 * np.sqrt(2) * self.B) * self.lrfrac

