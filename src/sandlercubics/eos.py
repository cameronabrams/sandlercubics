# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Pure-component cubic equations of state and their state calculations
"""
from __future__ import annotations
import logging
import numpy as np

from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import dataclass
from scipy.optimize import root_scalar

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
    """ units of temperature (not used) """

    phase: str = 'unspecified' 
    """ strict phase requirement flag; 'vapor', 'liquid', 'unspecified' (multiple roots) """

    P: float = None
    """ state pressure """
    T: float = None
    """ state temperature """
    x: float = None
    """ vapor fraction in two-phase region; None if single phase """
    Pc: float = None
    """ critical pressure """
    Tc: float = None
    """ critical temperature """
    omega: float = None 
    """ acentricity factor, dimensionless """
    Cp: float | list[float] | dict [str, float] | None = None
    """ heat capacity data for ideal-gas contributions """
    Tref: float = 298.15
    """ reference temperature for 'absolute' internal energy, enthalpy, and entropy calculations """
    Pref_MPa: float = 0.1
    """ reference pressure for 'absolute' entropy calculations, in MPa """

    logiter: bool = False
    """ flag for logging iterations in phase calculations """
    maxiter: int = 100
    """ maximum iterations in phase calculations """
    epsilon: float = 1.e-5
    """ fugacity tolerance in phase calculations (Fig. 7.5-1 in Sandler 5th ed) """
    iter: int = 0
    """ iteration counter in phase calculations """
    err: float = 0.0
    """ current error in phase calculations """

    Tref: float = 298.15
    """ reference temperature for 'absolute' internal energy, enthalpy,
    and entropy calculations """

    def clone(self) -> CubicEOS:
        """ Return a copy of the current EOS object """
        return type(self)(
            pressure_unit=self.pressure_unit,
            volume_unit=self.volume_unit,
            temperature_unit=self.temperature_unit,
            phase=self.phase,
            P=self.P,
            T=self.T,
            Pc=self.Pc,
            Tc=self.Tc,
            omega=self.omega,
            logiter=self.logiter,
            maxiter=self.maxiter,
            epsilon=self.epsilon,
            Tref=self.Tref,
            Pref_MPa=self.Pref_MPa,
            Cp=deepcopy(self.Cp)
        )

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
        """ dimensionless vdW 'a' parameter aP/(R^2 T^2) """
        return self.a * self.P / (self.R * self.T)**2
    
    @property
    def B(self):
        """ dimensionless vdW 'b' parameter bP/(RT) """
        return self.b * self.P / (self.R * self.T)
    
    @property
    @abstractmethod
    def cubic_coeff(self):
        """ coefficients of cubic equation for compressibility factor Z """
        pass

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
        """ Molar volume from compressibility factor Z """
        Z = self.Z
        return Z * self.R_pv * self.T / self.P

    @property
    @abstractmethod
    def h_departure(self):
        """ Enthalpy departure at state T and P """
        pass

    @property
    @abstractmethod
    def s_departure(self):
        """ Entropy departure at state T and P """
        pass
    
    @property
    @abstractmethod
    def logphi(self):
        """ natural log of fugacity coefficient at state T and P """
        pass
    
    @property
    def phi(self): # fugacity coefficient(s)
        logPhi = self.logphi
        return np.exp(logPhi)

    @property
    def f(self): # fugacity/fugacities
        return self.phi * self.P

    @property
    def Pvap(self):
        """ 
        Vapor pressure at current state temperature.  Current state pressure is retained. 
        Implements algorithm in Fig. 7.5-1 of Sandler 5th ed.
        """
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
            if self.logiter: logger.debug(f'Pvap: Iter {self.iter}: P {self.P:.6f}, fV {fV:.6f}, fL {fL:.6f}; error {self.err:.4e}')
            self.P *= fL / fV
            if self.err < self.epsilon or self.iter == self.maxiter:
                keepgoing = False
            if self.iter >= self.maxiter:
                logger.warning(f'Reached {self.iter} iterations without convergence; error {np.abs(fL/fV-1):.4e}')
        Pvap = self.P
        self.P = saveP
        return Pvap
    
    @property
    def Hvap(self):
        """ Heat of vaporization at current state temperature """
        satd_clone = self.clone()
        satd_clone.P = self.Pvap
        return satd_clone.h_departure[0] - satd_clone.h_departure[1]

    @property
    def Svap(self):
        """ Entropy of vaporization at current state temperature """
        satd_clone = self.clone()
        satd_clone.P = self.Pvap
        return satd_clone.s_departure[0] - satd_clone.s_departure[1]

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
            if self.logiter: logger.debug(f'Tsat: Iter {self.iter}: T {self.T:.6f}, fV {fV:.6f}, fL {fL:.6f}; error {self.err:.4e}')
            self.T *= (fV / fL)**0.125
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

    @property
    def h(self):
        """ absolute enthalpy at state T and P """
        if self.Cp is None:
            raise ValueError("Cp data required for absolute enthalpy calculation.")
        dH_ideal = DeltaH_IG(self.Tref, self.T, self.Cp)
        return self.h_departure + dH_ideal

    @property
    def u(self):
        """ absolute internal energy at state T and P """
        if self.Cp is None:
            raise ValueError("Cp data required for absolute internal energy calculation.")
        # u = h - pv
        return self.h - self.Pv * self.R / self.R_pv
    
    @property
    def Pref_local(self):
        """ reference pressure in local units """
        if self.pressure_unit == 'mpa':
            return self.Pref_MPa
        elif self.pressure_unit == 'bar':
            return self.Pref_MPa * 10.0
        elif self.pressure_unit == 'pa':
            return self.Pref_MPa * 1.e6
        elif self.pressure_unit == 'atm':
            return self.Pref_MPa * 9.86923
        else:
            raise ValueError(f"Unsupported pressure unit: {self.pressure_unit}")

    @property
    def s(self):
        """ absolute entropy at state T and P """
        if self.Cp is None:
            raise ValueError("Cp data required for absolute entropy calculation.")
        # make sure Pref_local is in correct units (same as self.P)
        dS_ideal = DeltaS_IG(self.Tref, self.Pref_local, self.T, self.P, self.Cp, self.R)
        return self.s_departure + dS_ideal

    def DeltaS(self, other: CubicEOS, Cp: float | list[float] | dict [str, float]):
        """
        Computes and returns entropy change from self to other state
        including ideal gas contribution based on Cp, no referece needed!
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

    def two_phase_check(self, v: float = None, h: float = None, s: float = None, u: float = None):
        """
        determine vapor fraction in two-phase region based on provided property
        """
        satd_clone = self.clone()
        if satd_clone.T is not None:
            satd_clone.P = self.Pvap
        elif satd_clone.P is not None:
            satd_clone.T = self.Tsat
        if v is not None:
            vL, vV = satd_clone.v
            return (v - vL) / (vV - vL)
        elif h is not None:
            hL, hV = satd_clone.h
            return (h - hL) / (hV - hL)
        elif s is not None:
            sL, sV = satd_clone.s
            return (s - sL) / (sV - sL)
        elif u is not None:
            uL, uV = satd_clone.u
            return (u - uL) / (uV - uL)
        else:
            raise ValueError("One of v, h, s, or u must be provided to determine vapor fraction.")

    def P_solve(self, P_bracket: list[float] = None, v: float = None, h: float = None, s: float = None, u: float = None):
        """
        Solve for pressure given T and one other property (v, h, s, u)
        """
        if v is not None:
            def objective(P):
                self.P = P
                logger.debug(f'Objective eval: P={P}, v={self.v}, target v={v}')
                if hasattr(self.v, '__len__'):
                    # pick the root closest to target v
                    v_diffs = np.abs(self.v - v)
                    v_closest = self.v[np.argmin(v_diffs)]
                    return v - v_closest
                return v - self.v
            res = root_scalar(objective, bracket=P_bracket)
            self.P = res.root
        elif h is not None:
            def objective(P):
                self.P = P
                logger.debug(f'Objective eval: P={P}, h={self.h}, target h={h}')
                if hasattr(self.h, '__len__'):
                    # pick the root closest to target h
                    h_diffs = np.abs(self.h - h)
                    h_closest = self.h[np.argmin(h_diffs)]
                    return h - h_closest
                return h - self.h
            res = root_scalar(objective, bracket=P_bracket)
            self.P = res.root
        elif s is not None:
            raise NotImplementedError("Pressure solve from entropy not implemented yet.")
        elif u is not None:
            raise NotImplementedError("Pressure solve from internal energy not implemented yet.")
        else:
            raise ValueError("One of v, h, s, or u must be provided to solve for pressure.")

    def T_solve(self, T_bracket: list[float] = None, v: float = None, h: float = None, s: float = None, u: float = None):
        """
        Solve for temperature given P and one other property (v, h, s, u)
        """
        if v is not None:
            def objective(T):
                self.T = T
                logger.debug(f'Objective eval: T={T}, v={self.v}, target v={v}')
                if hasattr(self.v, '__len__'):
                    # pick the root closest to target v
                    v_diffs = np.abs(self.v - v)
                    v_closest = self.v[np.argmin(v_diffs)]
                    return v - v_closest
                return v - self.v
            res = root_scalar(objective, bracket=T_bracket)
            self.T = res.root
        elif h is not None:
            def objective(T):
                self.T = T
                logger.debug(f'Objective eval: T={T}, h={self.h}, target h={h}')
                if hasattr(self.h, '__len__'):
                    # pick the root closest to target h
                    h_diffs = np.abs(self.h - h)
                    h_closest = self.h[np.argmin(h_diffs)]
                    return h - h_closest
                return h - self.h
            res = root_scalar(objective, bracket=T_bracket)
            self.T = res.root
        elif s is not None:
            raise NotImplementedError("Temperature solve from entropy not implemented yet.")
        elif u is not None:
            raise NotImplementedError("Temperature solve from internal energy not implemented yet.")
        else:
            raise ValueError("One of v, h, s, or u must be provided to solve for temperature.")

    def solve(self, T: float = None, P: float = None, v: float = None, h: float = None, s: float = None, u: float = None):
        """
        Solve for the missing state variable given two of T, P, v, h, s, u.
        Updates the object's state variables accordingly.
        """
        num_specified = sum(var is not None for var in [T, P, v, h, s, u])
        if num_specified != 2:
            raise ValueError("Exactly two of T, P, v, h, s, or u must be specified.")
        
        self.T = T
        self.P = P

        # Calculate the missing variable using the EOS
        if self.T is not None and self.P is not None:
            pass # do nothing, all props are set
        elif self.P is None: # T and one other property (not P) is set
            if self.T < self.Tc:
                x = self.two_phase_check(v=v, h=h, s=s, u=u)
                print(f'x={x}')
                if 0 < x < 1:
                    self.x = x
                    self.P = self.Pvap
                    # done.  All properties can be retrieved now.ABC
                    return
                elif x < 0:
                    # subcooled liquid
                    self.P = None
                    Plo = self.Pvap
                    Phi = self.Pc
                else:
                    # superheated vapor
                    Phi = self.Pvap
                    Plo = 1.e-5 # guess very low pressure
                self.P_solve([Plo, Phi], v=v, h=h, s=s, u=u)
            else: # given T is above Tc, so supercritical
                self.P_solve([1.e-5, self.Pc*5], v=v, h=h, s=s, u=u)
        elif self.T is None:
            if self.P < self.Pc:
                x = self.two_phase_check(v=v, h=h, s=s, u=u)
                if 0 < x < 1:
                    self.x = x
                    self.T = self.Tsat
                    return
                elif x < 0:
                    # subcooled liquid
                    self.T = None
                    Tlo = 1.0  # guess low temperature
                    Thi = self.Tsat
                else:
                    # superheated vapor
                    Thi = self.Tc
                    Tlo = self.Tsat
                self.T_solve([Tlo, Thi], v=v, h=h, s=s, u=u)
            else:
                self.T_solve([10.0, self.Tc+100], v=v, h=h, s=s, u=u)
        else:
            raise NotImplementedError("Solving when both T and P are not set is not implemented yet.")

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
    def cubic_coeff(self):
        """ cubic coefficients for ideal gas """
        return np.array([0.0, 0.0, 1.0, -1.0])

    @property
    def Z(self):
        """ compressibility factor for ideal gas is unity """
        return 1.0
    
    @property
    def h_departure(self):
        """ Enthalpy departure at state T and P; ideal-gas value """
        return 0.0

    @property
    def s_departure(self):
        """ Entropy departure at state T and P; ideal-gas value """
        return 0.0
    
    @property
    def logphi(self):
        """ natural log of fugacity coefficient at state T and P; ideal-gas value """
        return 0.0

@dataclass
class VanDerWaalsEOS(CubicEOS):
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

@dataclass
class SoaveRedlichKwongEOS(CubicEOS):
    """
    Pure-component Soave-Redlich-Kwong equation of state (Eq. 6.4-1 in Sandler 5th ed)
    """
    @property
    def kappa(self):
        return 0.480 + 1.574 * self.omega - 0.176 * self.omega**2
    
    @property
    def alpha(self):
        return (1 + self.kappa * (1 - np.sqrt(self.T/self.Tc)))**2

    @property
    def a(self):
        return 0.42748 * self.R**2 * self.Tc**2 / self.Pc * self.alpha
    
    @property
    def b(self):
        return 0.08664 * self.R * self.Tc / self.Pc

    @property
    def cubic_coeff(self):
        return np.array([1.0, -1.0, self.A - self.B - self.B**2, -self.A * self.B])

    @property
    def lrfrac(self):
        z = self.Z
        num_arg = z + self.B
        den_arg = z
        return np.log(num_arg / den_arg)
    
    @property
    def h_departure(self):
        """
        Enthalpy departure at state T and P (from solution to problem 6.36 in Sandler 5th ed)
        """
        z = self.Z
        return self.R * self.T * (z - 1) + (self.T * self.da_dT - self.a) / self.b * self.lrfrac

    @property
    def s_departure(self):
        """
        Entropy departure at state T and P (from solution to problem 6.36 in Sandler 5th ed)
        """
        z = self.Z
        return self.R * np.log(z - self.B) + self.da_dT/self.b * self.lrfrac

    @property
    def logphi(self):
        """
        natural log of fugacity coefficient at state T and P (from solution to problem 7.46 in Sandler 5th ed)
        """
        z = self.Z
        return z - 1 - np.log(z - self.B) - self.a / (self.R * self.T * self.b) * self.lrfrac