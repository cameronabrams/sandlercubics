# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Pure-component cubic equations of state and their state calculations
"""
from __future__ import annotations
import logging
import numpy as np

from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import dataclass, field
from scipy.optimize import root_scalar

from sandlermisc.gas_constant import GasConstant
from sandlermisc.thermals import DeltaH_IG, DeltaS_IG
from sandlerprops.compound import Compound
from sandlerprops.properties import get_database

logger = logging.getLogger(__name__)

sqrt_2 = np.sqrt(2)

@dataclass
class CubicEOS(ABC):
    """
    Abstract class for all Cubic equations of state.
    """
    _cache: dict = field(default_factory=dict, init=False, repr=False)
    """ Cache for computed properties """
    _input_state: dict | None = field(default=None, init=False, repr=False)
    """ Snapshot of input field values for cache validation """

    _INPUT_FIELDS = frozenset(['T', 'P', 'x', 'Tc', 'Pc', 'omega', 'Cp', 'Tref', 'Pref_MPa', 'phase', 'pressure_unit', 'volume_unit'])
    """ Fields that affect calculations and cache validity """

    def __setattr__(self, name, value):
        """ Clear cache on input field changes """
        if name in self._INPUT_FIELDS and hasattr(self, '_cache'):
            self._cache.clear()
        super().__setattr__(name, value)

    def _get_current_input_state(self):
        """ Snapshot of current input field values. """
        return {field: getattr(self, field) for field in self._INPUT_FIELDS}

    pressure_unit: str = 'MPa' # MPa
    """ units of pressure, lower-case """
    volume_unit: str = 'm3'
    """ units of volume, lower-case """

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

    def is_cache_stale(self, property_name: str = None):
        """
        Check if cache is stale.
        
        Parameters
        ----------
        property_name: str
            Specific property to check, or None for general staleness
        
        Returns
        -------
        bool or dict: True/False if property_name given, else dict of changed inputs
        """
        if property_name:
            return property_name not in self._cache
        
        # Show which inputs changed
        if self._input_state is None:
            return self._INPUT_FIELDS  # All inputs "changed" (initial state)
        
        current = self._get_current_input_state()
        changed = {k for k in self._INPUT_FIELDS if current[k] != self._input_state.get(k)}
        return changed if changed else False

    def clone(self) -> CubicEOS:
        """ Return a copy of the current EOS object """
        return type(self)(
            pressure_unit=self.pressure_unit,
            volume_unit=self.volume_unit,
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

    @abstractmethod
    def _calc_a(self):
        pass

    @property
    def a(self):
        """ van der Waals a parameter """
        if not 'a' in self._cache: # cache was cleared on an input change
            self._cache['a'] = self._calc_a()
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['a']
    
    def _calc_da_dT(self):
        return 0.0

    @property
    def da_dT(self):
        """ derivative of vdW a parameter wrt temperature """
        if not 'da_dT' in self._cache:
            self._cache['da_dT'] = self._calc_da_dT()
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['da_dT']

    @abstractmethod
    def _calc_b(self):
        pass

    @property
    def b(self):
        """ van der Waals b parameter """
        if not 'b' in self._cache:
            self._cache['b'] = self._calc_b()
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['b']
    
    @property
    def A(self):
        """ dimensionless vdW a parameter aP/(R^2 T^2) """
        if not 'A' in self._cache:
            self._cache['A'] = self.a * self.P / (self.R * self.T)**2
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['A']

    @property
    def B(self):
        """ dimensionless vdW b parameter bP/(RT) """
        if not 'B' in self._cache:
            self._cache['B'] = self.b * self.P / (self.R * self.T)
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['B']

    @property
    @abstractmethod
    def cubic_coeff(self):
        """ coefficients of cubic equation for compressibility factor Z """
        pass

    @property
    def Z(self) -> np.ndarray:
        """
        Computes or retrieves compressibility factor(s) at current state T and P.

        Returns
        -------
        np.ndarray
            Compressibility factor(s) at current state T and P
        """
        if 'Z' not in self._cache:
            self._cache['Z'] = self._solve_for_Z()
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['Z']

    def _solve_for_Z(self) -> np.ndarray:
        """
        Solve the cubic equation for compressibility factor(s) Z if not cached.

        Returns
        -------
        np.ndarray
            Real roots of the cubic equation representing compressibility factors.
        """
        complx_roots = np.roots(self.cubic_coeff)
        real_roots_idx = np.where(complx_roots.imag==0)[0]
        real_roots = complx_roots[real_roots_idx].real
        if len(real_roots) == 1:
            Z = np.array([real_roots[0]])
        else:
            if self.phase == 'vapor':
                Z = np.array([np.max(real_roots)])
            elif self.phase == 'liquid':
                Z = np.array([np.min(real_roots)])
            else: # unspecified/multiple
                Z = np.array([real_roots[0], real_roots[2]]) # sorted by numpy.roots
        return Z

    @property
    def v(self) -> np.ndarray:
        """
        Computes or retrieves molar volume(s) at current state T and P.
        
        Returns
        -------
        np.ndarray
            Molar volume(s) at current state T and P
        """
        if not 'v' in self._cache:
            self._cache['v'] = self.Z * self.R_pv * self.T / self.P
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['v']

    @abstractmethod
    def _calc_h_departure(self) -> np.ndarray:
        """
        Calculate the enthalpy departure at state T and P.
        """
        pass

    @property
    def h_departure(self) -> np.ndarray:
        """
        Computes or retrieves enthalpy departure at current state T and P.

        Returns
        -------
        np.ndarray
            Enthalpy departure at current state T and P
        """
        if not 'h_departure' in self._cache:
            self._cache['h_departure'] = self._calc_h_departure()
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['h_departure']

    @abstractmethod
    def _calc_s_departure(self) -> np.ndarray:
        """
        Calculate the entropy departure at state T and P.
        """
        pass

    @property
    def s_departure(self) -> np.ndarray:
        """
        Computes or retrieves entropy departure at current state T and P.
        
        Returns
        -------
        np.ndarray
            Entropy departure at current state T and P
        """
        if not 's_departure' in self._cache:
            self._cache['s_departure'] = self._calc_s_departure()
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['s_departure']
    
    @abstractmethod
    def _calc_logphi(self) -> np.ndarray:
        """
        Calculate the natural log of fugacity coefficient at state T and P.
        """
        pass

    @property
    def logphi(self) -> np.ndarray:
        """
        Computes or retrieves the natural log of fugacity coefficient at current state T and P.

        Returns
        -------
        np.ndarray
            Natural log of fugacity coefficient(s) at current state T and P
        """
        if not 'logphi' in self._cache:
            self._cache['logphi'] = self._calc_logphi()
            self._input_state = self._get_current_input_state()  # Save state
        return self._cache['logphi']
    
    @property
    def phi(self) -> np.ndarray:
        """
        Computes or retrieves fugacity coefficient(s) at current state T and P.
        
        Returns
        -------
        np.ndarray
            Fugacity coefficient(s) at current state T and P
        """
        if not 'phi' in self._cache:
            self._cache['phi'] = np.exp(self.logphi)
        return self._cache['phi']

    @property
    def f(self) -> np.ndarray:
        """
        Computes or retrieves fugacity/fugacities at current state T and P.

        Returns
        -------
        np.ndarray
            Fugacity/fugacities at current state T and P
        """
        if not 'f' in self._cache:
            self._cache['f'] = self.phi * self.P
        return self._cache['f']

    @property
    def Pvap(self) -> float:
        """ 
        Computes or retrieves vapor pressure at current state temperature.  If not calculable, returns np.nan.
        Implements algorithm in Fig. 7.5-1 of Sandler 5th ed.
        """
        if not 'pvap' in self._cache:
            if self.Pc is None or self.Tc is None:
                self._cache['pvap'] = np.nan
                return self._cache['pvap']
            tmp = self.clone()
            tmp.P = tmp.Pc * (tmp.T / tmp.Tc)**8
            pu = self.R._capitalizations.get(self.pressure_unit, self.pressure_unit)
            logger.debug(f'Initial Pvap guess: {tmp.P:.6f} {pu} at T {tmp.T:.6f} (Tc {tmp.Tc:.6f}): f {tmp.f}')
            keepgoing = True
            i = 0
            while keepgoing:
                i += 1
                try:
                    fV, fL = tmp.f
                except:
                    # raise ValueError(f'Pvap calc error: T {tmp.T:.2f}, P {tmp.P:.2f} {pu}, Tc {tmp.Tc:.2f}, Pc {tmp.Pc:.2f} {pu}')
                    self._cache['pvap'] = np.nan
                    return self._cache['pvap']
                err = np.abs(fL / fV - 1)
                if tmp.logiter: logger.debug(f'Pvap: Iter {i}: P {tmp.P:.6f} {pu}, fV {fV:.6f}, fL {fL:.6f}; error {err:.4e}')
                tmp.P *= fL / fV
                if err < tmp.epsilon or i == tmp.maxiter:
                    keepgoing = False
                if i >= tmp.maxiter:
                    logger.warning(f'Reached {i} iterations without convergence; error {err:.4e}')
            self._cache['pvap'] = tmp.P
        return self._cache['pvap']
    
    @property
    def Hvap(self) -> float:
        """ 
        Computes or retrieves enthalpy of vaporization at current state temperature.

        Returns
        -------
        float
            Enthalpy of vaporization at current state temperature. If not calculable, returns np.nan.
        """
        if not 'hvap' in self._cache:
            satd_clone = self.clone()
            satd_clone.P = self.Pvap
            if not np.isnan(satd_clone.P):
                self._cache['hvap'] = satd_clone.h_departure[0] - satd_clone.h_departure[1]
            else:
                self._cache['hvap'] = np.nan
        return self._cache['hvap']
    
    @property
    def Svap(self) -> float:
        """ 
        Computes or retrieves entropy of vaporization at current state temperature.

        Returns
        -------
        float
            Entropy of vaporization at current state temperature.  If not calculable, returns np.nan.
        """
        if not 'svap' in self._cache:
            satd_clone = self.clone()
            satd_clone.P = self.Pvap
            if not np.isnan(satd_clone.P):
                self._cache['svap'] = satd_clone.s_departure[0] - satd_clone.s_departure[1]
            else:
                self._cache['svap'] = np.nan
        return self._cache['svap']

    @property
    def Tsat(self) -> float:
        """ 
        Computes or retrieves saturation temperature at state pressure.
        
        Returns
        -------
        float
            Saturation temperature at state pressure.  If not calculable, returns np.nan.
        """
        if not 'Tsat' in self._cache:
            tmp = self.clone()
            if self.Pc is None or self.Tc is None:
                self._cache['Tsat'] = np.nan
                return self._cache['Tsat']
            tmp.T = tmp.Tc * (tmp.P / tmp.Pc)**0.125
            pu = self.R._capitalizations.get(self.pressure_unit, self.pressure_unit)
            logger.debug(f'Initial Tsat guess: {tmp.T:.6f} at P {tmp.P:.6f} {pu}')
            keepgoing = True
            i = 0
            while keepgoing:
                i += 1
                try:
                    fV, fL = tmp.f
                except:
                    # raise ValueError(f'Error in fugacity calculation in Tsat at T {tmp.T:.2f}, P {tmp.P:.2f} {pu}')
                    self._cache['Tsat'] = np.nan
                    return self._cache['Tsat']
                err = np.abs(fV / fL - 1)
                if tmp.logiter: logger.debug(f'Tsat: Iter {i}: T {tmp.T:.6f}, fV {fV:.6f} {pu}, fL {fL:.6f} {pu}; error {err:.4e}')
                tmp.T *= (fV / fL)**0.125
                if err < tmp.epsilon or i == tmp.maxiter:
                    keepgoing = False
                if i == tmp.maxiter:
                    logger.warning(f'Reached {i} iterations without convergence; error {err:.4e}')
            self._cache['Tsat'] = tmp.T
        return self._cache['Tsat']

    def unit_consistency(self, other: CubicEOS):
        consistent = self.pressure_unit == other.pressure_unit and self.volume_unit == other.volume_unit
        if not consistent:
            raise ValueError('inconsistent units')

    @property
    def h(self) -> np.ndarray:
        """
        Computes absolute enthalpy at state T and P
        
        Returns
        -------
        np.ndarray
            Absolute enthalpy/enthalpies at state T and P
        """
        if self.Cp is None:
            raise ValueError("Cp data required for absolute enthalpy calculation.")
        dH_ideal = DeltaH_IG(self.Tref, self.T, self.Cp)
        logger.debug(f'dh_ideal: {dH_ideal} for Tref {self.Tref} to T {self.T} with Cp {self.Cp}')
        return self.h_departure + dH_ideal

    @property
    def u(self) -> np.ndarray:
        """
        Computes absolute internal energy at state T and P
        
        Returns
        -------
        np.ndarray
            Absolute internal energy/energies at state T and P
        """
        if self.Cp is None:
            raise ValueError("Cp data required for absolute internal energy calculation.")
        # u = h - pv
        return self.h - self.Pv * self.R / self.R_pv
    
    @property
    def Pref_local(self) -> float:
        """ 
        Returns reference pressure in local units

        Returns
        -------
        float
            Reference pressure in local units
        """
        if self.pressure_unit == 'mpa' or self.pressure_unit == 'MPa':
            return self.Pref_MPa
        elif self.pressure_unit == 'bar':
            return self.Pref_MPa * 10.0
        elif self.pressure_unit == 'pa' or self.pressure_unit == 'Pa':
            return self.Pref_MPa * 1.e6
        elif self.pressure_unit == 'kpa' or self.pressure_unit == 'kPa':
            return self.Pref_MPa * 1.e3
        elif self.pressure_unit == 'atm':
            return self.Pref_MPa * 9.86923
        else:
            raise ValueError(f"Unsupported pressure unit: {self.pressure_unit}")

    @property
    def s(self) -> np.ndarray:
        """
        Computes absolute entropy at state T and P
        
        Returns
        -------
        np.ndarray
            Absolute entropy/entropies at state T and P
        """
        if self.Cp is None:
            raise ValueError("Cp data required for absolute entropy calculation.")
        # make sure Pref_local is in correct units (same as self.P)
        dS_ideal = DeltaS_IG(self.Tref, self.Pref_local, self.T, self.P, self.Cp, self.R)
        return self.s_departure + dS_ideal

    def delta_h(self, other: CubicEOS) -> np.ndarray:
        """
        Computes and returns enthalpy change from self to other state
        """
        self.unit_consistency(other)
        return other.h - self.h

    def delta_s(self, other: CubicEOS) -> np.ndarray:
        """
        Computes and returns entropy change from self to other state
        """
        self.unit_consistency(other)
        return other.s - self.s
    
    def delta_pv(self, other: CubicEOS) -> np.ndarray:
        """
        Returns Delta(PV) in thermal (not PV) units 
        """
        self.unit_consistency(other)
        return (other.Pv - self.Pv) * self.R / self.R_pv
    
    def delta_u(self, other: CubicEOS) -> np.ndarray:
        """
        Returns Delta(U) (internal energy)
        """
        return self.delta_h(other) - self.delta_pv(other)
        
    def two_phase_check(self, v: float = None, h: float = None, s: float = None, u: float = None):
        """
        Returns the vapor fraction given one property in the two-phase region.

        Parameters
        ----------
        v: float
            Molar volume
        h: float
            Enthalpy
        s: float
            Entropy
        u: float
            Internal energy

        Returns
        -------
        float
            Vapor fraction, where 0 to 1 signfifies two-phase region, <0 is subcooled liquid, >1 is superheated vapor
        """
        satd_clone = self.clone()
        if satd_clone.T is not None:
            satd_clone.P = self.Pvap
        elif satd_clone.P is not None:
            satd_clone.T = self.Tsat
        if v is not None:
            vV, vL = satd_clone.v
            return (v - vL) / (vV - vL)
        elif h is not None:
            hV, hL = satd_clone.h
            return (h - hL) / (hV - hL)
        elif s is not None:
            sV, sL = satd_clone.s
            return (s - sL) / (sV - sL)
        elif u is not None:
            uV, uL = satd_clone.u
            return (u - uL) / (uV - uL)
        else:
            raise ValueError("One of v, h, s, or u must be provided to determine vapor fraction.")

    def P_solve(self, P_bracket: list[float] = None, v: float = None, h: float = None, s: float = None, u: float = None):
        """
        Solve for pressure given T and one other property (v, h, s, u). Stores result in self.P.

        Parameters
        ----------
        P_bracket: list[float]
            Bracket [P_low, P_high] for root finding
        v: float
            Molar volume
        h: float
            Enthalpy
        s: float
            Entropy
        u: float
            Internal energy
        """
        if v is not None:
            def objective(P):
                self.P = P
                logger.debug(f'Objective eval: P={P}, v={self.v}, target v={v}')
                # pick the root closest to target v
                v_diffs = np.abs(self.v - v)
                v_closest = self.v[np.argmin(v_diffs)]
                return v - v_closest
        elif h is not None:
            def objective(P):
                self.P = P
                logger.debug(f'Objective eval: P={P}, h={self.h}, target h={h}')
                    # pick the root closest to target h
                h_diffs = np.abs(self.h - h)
                h_closest = self.h[np.argmin(h_diffs)]
                return h - h_closest
        elif s is not None:
            def objective(P):
                self.P = P
                logger.debug(f'Objective eval: P={P}, s={self.s}, target s={s}')
                # pick the root closest to target s
                s_diffs = np.abs(self.s - s)
                s_closest = self.s[np.argmin(s_diffs)]
                return s - s_closest
        elif u is not None:
            def objective(P):
                self.P = P
                logger.debug(f'Objective eval: P={P}, u={self.u}, target u={u}')
                # pick the root closest to target u
                u_diffs = np.abs(self.u - u)
                u_closest = self.u[np.argmin(u_diffs)]
                return u - u_closest
        else:
            raise ValueError("One of v, h, s, or u must be provided to solve for pressure.")
        res = root_scalar(objective, bracket=P_bracket)
        self.P = res.root

    def T_solve(self, T_bracket: list[float] = None, v: float = None, h: float = None, s: float = None, u: float = None):
        """
        Solve for temperature given P and one other property (v, h, s, u). Stores result in self.T.

        Parameters
        ----------
        T_bracket: list[float]
            Bracket [T_low, T_high] for root finding
        v: float
            Molar volume
        h: float
            Enthalpy
        s: float
            Entropy
        u: float
            Internal energy
        """
        if v is not None:
            def objective(T):
                self.T = T
                logger.debug(f'Objective eval: T={T}, v={self.v}, target v={v}')
                # pick the root closest to target v
                v_diffs = np.abs(self.v - v)
                v_closest = self.v[np.argmin(v_diffs)]
                return v - v_closest
        elif h is not None:
            def objective(T):
                self.T = T
                logger.debug(f'Objective eval: T={T}, h={self.h}, target h={h}')
                # pick the root closest to target h
                h_diffs = np.abs(self.h - h)
                h_closest = self.h[np.argmin(h_diffs)]
                return h - h_closest
        elif s is not None:
            def objective(T):
                self.T = T
                logger.debug(f'Objective eval: T={T}, s={self.s}, target s={s}')
                # pick the root closest to target s
                s_diffs = np.abs(self.s - s)
                s_closest = self.s[np.argmin(s_diffs)]
                return s - s_closest
        elif u is not None:
            def objective(T):
                self.T = T
                logger.debug(f'Objective eval: T={T}, u={self.u}, target u={u}')
                # pick the root closest to target u
                u_diffs = np.abs(self.u - u)
                u_closest = self.u[np.argmin(u_diffs)]
                return u - u_closest
        else:
            raise ValueError("One of v, h, s, or u must be provided to solve for temperature.")

        res = root_scalar(objective, bracket=T_bracket)
        self.T = res.root

    def solve(self, T: float = None, P: float = None, v: float = None, h: float = None, s: float = None, u: float = None):
        """
        Solve for the missing state variable given two of T, P, v, h, s, u.
        Updates the object's state variables accordingly.

        Parameters
        ----------
        T: float
            Temperature
        P: float
            Pressure
        v: float
            Molar volume
        h: float
            Enthalpy
        s: float
            Entropy
        u: float
            Internal energy
        """
        num_specified = sum(var is not None for var in [T, P, v, h, s, u])
        if num_specified != 2:
            raise ValueError("Exactly two of T, P, v, h, s, or u must be specified.")
        
        self.T = T
        self.P = P

        # Calculate the missing variable using the EOS
        if self.T is not None and self.P is not None:
            pass # all properties can be retrieved now
        elif self.P is None: # T and one other property (not P) is set
            if self.T < self.Tc:
                x = self.two_phase_check(v=v, h=h, s=s, u=u)
                if 0 < x < 1:
                    self.x = x
                    self.P = self.Pvap
                    # done.  All properties can be retrieved now.
                    return
                elif x < 0:
                    # subcooled liquid
                    Plo = self.Pvap
                    Phi = self.Pc
                    phase = 'liquid'
                else:
                    # superheated vapor
                    Phi = self.Pvap
                    Plo = 1.e-5 # guess very low pressure
                    phase = 'vapor'
                self.phase = phase
                self.P_solve([Plo, Phi], v=v, h=h, s=s, u=u)
                self.phase = 'unspecified'
            else: # given T is above Tc, so supercritical
                self.P_solve([1.e-5, self.Pc*5], v=v, h=h, s=s, u=u)
        elif self.T is None:
            if self.P < self.Pc:
                x = self.two_phase_check(v=v, h=h, s=s, u=u)
                logger.debug(f'P {self.P} subcrit: x = {x}')
                if 0 < x < 1:
                    self.x = x
                    self.T = self.Tsat
                    return
                elif x < 0:
                    # subcooled liquid
                    Tlo = 1.0  # guess low temperature
                    Thi = self.Tsat
                    phase = 'liquid'
                else:
                    # superheated vapor
                    Thi = self.Tc*3
                    Tlo = self.Tsat
                    phase = 'vapor'
                self.phase = phase
                self.T_solve([Tlo, Thi], v=v, h=h, s=s, u=u)
                self.phase = 'unspecified'
            else:
                self.T_solve([10.0, self.Tc*3], v=v, h=h, s=s, u=u)
        else:
            raise NotImplementedError("Solving when neither T nor P are set is not implemented yet.")

    def set_compound(self, compound: str | Compound):
        """
        Set critical properties and Cp data from a compound name.

        Parameters
        ----------
        compound: str | Compound
            Name of the compound to retrieve properties for
        """
        db = get_database()
        compound = db.get_compound(compound) if isinstance(compound, str) else compound
        compound_name = compound.Name if compound is not None else str(compound)
        if compound is None:
            raise ValueError(f"Compound '{compound_name}' not found in database.")
        return self.transfer_crits_from_compound(compound)

    def transfer_crits_from_compound(self, compound: Compound = None):
        """
        Set critical properties and Cp data from a Compound object.

        Parameters
        ----------
        compound: Compound
            Compound object containing critical properties and Cp data
        """
        if compound is not None:
            self.Tc = compound.Tc
            Pc_bar = compound.Pc
            if self.pressure_unit == 'MPa' or self.pressure_unit == 'mpa':
                self.Pc = Pc_bar / 10.0
            elif self.pressure_unit == 'bar':
                self.Pc = Pc_bar
            elif self.pressure_unit == 'kPa' or self.pressure_unit == 'kpa':
                self.Pc = Pc_bar * 100.0
            elif self.pressure_unit == 'Pa' or self.pressure_unit == 'pa':
                self.Pc = Pc_bar * 1.e5
            elif self.pressure_unit == 'atm':
                self.Pc = Pc_bar / 1.01325
            else:
                raise ValueError(f"Unsupported pressure unit: {self.pressure_unit}")
            self.omega = compound.Omega
            self.Cp = deepcopy(compound.Cp)
        return self