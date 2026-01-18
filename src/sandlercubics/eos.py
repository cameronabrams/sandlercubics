# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Pure-component cubic equations of state and their state calculations
"""
from __future__ import annotations
import logging
import numpy as np
import pint
from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import dataclass, field
from scipy.optimize import root_scalar

from sandlermisc.thermals import DeltaH_IG, DeltaS_IG, R
from sandlermisc.thermodynamicstate import ThermodynamicState, _ureg
from sandlerprops.compound import Compound
from sandlerprops.properties import get_database

logger = logging.getLogger(__name__)

sqrt_2 = np.sqrt(2)

@dataclass  
class CubicEOS(ThermodynamicState):
    """
    Abstract class for all Cubic equations of state.
    """
    phase: str = 'unspecified' 
    """ strict phase requirement flag; 'vapor', 'liquid', 'unspecified' (multiple roots) """

    Pc: pint.Quantity = None
    """ critical pressure """
    Tc: pint.Quantity = None
    """ critical temperature """
    omega: float = None
    """ acentricity factor, dimensionless """

    Cp: float | list[float] | dict [str, float] | None = None
    """ heat capacity data for ideal-gas contributions """

    _PARAMETER_FIELDS = frozenset({'Tc', 'Pc', 'omega', 'Cp'})

    Z: float = None
    """ compressibility factor """
    
    Tref: pint.Quantity = 298.15 * _ureg.kelvin
    """ reference temperature for 'absolute' internal energy, enthalpy, and entropy calculations """
    Pref: pint.Quantity = 0.1 * _ureg.megapascal
    """ reference pressure for 'absolute' entropy calculations, in MPa """

    logiter: bool = False
    """ flag for logging iterations in phase calculations """
    maxiter: int = 100
    """ maximum iterations in phase calculations """
    epsilon: float = 1.e-5
    """ fugacity tolerance in phase calculations (Fig. 7.5-1 in Sandler 5th ed) """

    _default_unit_map = {
        'P': 'MPa',
        'T': 'K',
        'v': 'm**3 / mol',
        'u': 'J / mol',
        'h': 'J / mol',
        's': 'J / (mol * K)',
        'Pv': 'J / mol',
    }

    _pint_unit_map = {
        P: _ureg.megapascal,
        T: _ureg.kelvin,
        v: _ureg.meter**3 / _ureg.mol,
        u: _ureg.joule / _ureg.mol,
        h: _ureg.joule / _ureg.mol,
        s: _ureg.joule / (_ureg.mol * _ureg.kelvin),
        Pv: _ureg.joule / _ureg.mol,
    }

    def __setattr__(self, name, value):
        """
        Override __setattr__ to clear cached properties when inputs change.
        """
        super().__setattr__(name, value)
        if name in self._PARAMETER_FIELDS:
            self._cache['_is_complete'] = False

    def _resolve(self):
        """
        Resolve all state variables from the two input variables.   This method
        is called automatically when any property assignment results
        in a fully specified state.  There should never be a situation in which
        this method is called explicitly; it is only called by __setattr__ or __post_init__.
        """
        states_speced = list(self._cache['_input_values'].keys())
        logger.debug(f'_resolve: State {self.name}: Attempting to resolve state with inputs: {states_speced}')
        hasT, hasP, hasx = 'T' in states_speced, 'P' in states_speced, 'x' in states_speced
        if hasT and hasP:
            assert not hasx, f'State {self.name}: Cannot specify T, P, and x simultaneously.'
            resolved = self._resolve_TP()
        elif hasT:
            if hasx: # explicitly saturated state
                resolved = self._resolve_saturated_Tx()
            else: # may be saturated or unsaturated
                resolved = self._resolve_Top()
        elif hasP:
            if hasx: # explicitly saturated state
                resolved = self._resolve_saturated_Px()
            else: # may be saturated or unsaturated
                resolved = self._resolve_Pop()
        elif hasx:
            resolved = self._resolve_saturated_xop()
        else: 
            resolved = self._resolve_op_op()
        if not resolved:
            raise ValueError(f'State {self.name}: Unable to resolve state with inputs: {states_speced}')
        self._scalarize()
        logger.debug(f'_resolve: State {self.name}: Successfully resolved state with inputs: {states_speced}')
        self._cache['_is_complete'] = True
        logger.debug(f'_resolve: State {self.name}: State resolution complete {self._cache["_is_complete"]}')

    def _resolve_TP(self) -> bool:
        """
        Resolve state given T and P.
        """
        Zlist = self._solve_for_Z()
        if len(Zlist) == 1:
            Z = Zlist[0]
            setattr(self, 'Z', Z)
        else:
            if self.T < self.Tc:
                Psat = self.Psat
                if self.P < Psat:
                    setattr(self, 'Z', Zlist[0])  # vapor root
                elif self.P >= Psat:
                    setattr(self, 'Z', Zlist[1])   # liquid root
            else:  # P < Pc
                Tsat = self.Tsat
                if self.T < Tsat:
                    setattr(self, 'Z', Zlist[1])   # liquid root
                elif self.T >= Tsat:
                    setattr(self, 'Z', Zlist[0])  # vapor root
        return self._calculate_vhus()

    def _calculate_vhus(self) -> bool:
        """
        Resolve v, h, u, s from Z, T, P
        """
        vunit = self.get_default_unit('v')
        punit = self.get_default_unit('P')
        tunit = self.get_default_unit('T')
        eunit = self.get_default_unit('h')
        sunit = self.get_default_unit('s')
        v = self.Z * R.to(punit*vunit/tunit) * self.T / self.P
        if self.Cp is None:
            raise ValueError("Cp data required for absolute enthalpy/entropy calculation.")
        dH_ideal = DeltaH_IG(self.Tref.to('K'), self.T.to('K'), self.Cp).to(eunit)
        h = self.h_departure + dH_ideal
        dS_ideal = DeltaS_IG(self.Tref.to('K'), self.Pref.to(punit), 
                                self.T.to('K'), self.P.to(punit), 
                                self.Cp, R.to(eunit/tunit)).to(sunit)
        s = self.s_departure + dS_ideal
        Pv = (self.P * v).to(eunit)
        u = h - Pv

        setattr(self, 'v', v)
        setattr(self, 'h', h)
        setattr(self, 'u', u)
        setattr(self, 's', s)
        setattr(self, 'Pv', Pv)

        return True

    def _resolve_saturated_Tx(self) -> bool:
        specs = list(self._cache['_input_values'].keys())
        setattr(self, 'P', self.Psat)
        Z = self._solve_for_Z()
        if 0 < self.x < 1:
            self.Liquid = CubicEOS(x=0.0, T=self.T, 
                name=f'{self.name}_L' if self.name else 'Saturated Liquid',
                Tc=self.Tc, Pc=self.Pc, omega=self.omega, Cp=self.Cp)
            self.Vapor = CubicEOS(x=1.0, T=self.T, 
                name=f'{self.name}_V' if self.name else 'Saturated Vapor',
                Tc=self.Tc, Pc=self.Pc, omega=self.omega, Cp=self.Cp)
            for op in self._STATE_VAR_FIELDS - {'T', 'P', 'x'}:
                setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
            setattr(self, 'Z', None)
        elif self.x == 1.0:
            setattr(self, 'Z', Z[0])
            self._calculate_vhus()
        elif self.x == 0.0:
            setattr(self, 'Z', Z[1])
            self._calculate_vhus()

    def _resolve_saturated_Px(self) -> bool:
        specs = list(self._cache['_input_values'].keys())
        setattr(self, 'T', self.Tsat)
        Z = self._solve_for_Z()
        if 0 < self.x < 1:
            self.Liquid = CubicEOS(x=0.0, P=self.P, name=f'{self.name}_L' if self.name else 'Saturated Liquid', Tc=self.Tc, Pc=self.Pc, omega=self.omega, Cp=self.Cp)
            self.Vapor = CubicEOS(x=1.0, P=self.P, name=f'{self.name}_V' if self.name else 'Saturated Vapor', Tc=self.Tc, Pc=self.Pc, omega=self.omega, Cp=self.Cp)
            for op in self._STATE_VAR_FIELDS - {'T', 'P', 'x'}:
                setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
        elif self.x == 1.0:
            setattr(self, 'Z', Z[0])
            self._calculate_vhus()
        elif self.x == 0.0:
            setattr(self, 'Z', Z[1])
            self._calculate_vhus()

    def _resolve_Top(self) -> bool:
        """
        Resolve state given T and one other property.  Input state is set, so any state variable assignments not to the inputs will be pass-through settings.
        """
        punit = self.get_default_unit('P')
        specs = list(self._cache['_input_values'].keys())
        op = specs[0] if specs[1] == 'T' else specs[1]
        op_value = getattr(self, op)
        Psat = self.Psat
        constants = {k: getattr(self, k) for k in self._PARAMETER_FIELDS }
        satd_liquid_state = CubicEOS(x=0.0, T=self.T, name=f'{self.name}_satL', **constants)
        satd_vapor_state = CubicEOS(x=1.0, T=self.T, name=f'{self.name}_satV', **constants)
        op_value_sat_vapor = getattr(satd_vapor_state, op)
        op_value_sat_liquid = getattr(satd_liquid_state, op)
        companion_state = CubicEOS(T=self.T, P=Psat, name=f'{self.name}_companion', **constants)
        def objective(P):
            setattr(companion_state, 'P', P)
            return getattr(companion_state, op) - op_value
        if op_value > op_value_sat_vapor: # superheated vapor
            # find P that gives desired property value at T

            # P always > Psat, so all trials will be in superheated region
            p_bracket = [Psat, self.Pc.m_as(punit)*10.0]
        else:
            p_bracket = [0.01 * self.Pc.m_as(punit), Psat]

        sol = root_scalar(objective, bracket=p_bracket, method='bisect', xtol=self.epsilon, maxiter=self.maxiter)
        if not sol.converged:
            raise ValueError(f'State {self.name}: Unable to converge to solve for P in unsaturated region given T and {op}.')
        setattr(self, 'P', sol.root * punit)
        for var in self._STATE_VAR_FIELDS - {'T', 'P', 'x', op}:
            setattr(self, var, getattr(companion_state, var))
        return True

    def _resolve_Pop(self) -> bool:
        tunit = self.get_default_unit('T')
        specs = list(self._cache['_input_values'].keys())
        op = specs[0] if specs[1] == 'P' else specs[1]
        op_value = getattr(self, op)
        Tsat = self.Tsat
        constants = {k: getattr(self, k) for k in self._PARAMETER_FIELDS }
        satd_liquid_state = CubicEOS(x=0.0, P=self.P, name=f'{self.name}_satL', **constants)
        satd_vapor_state = CubicEOS(x=1.0, P=self.P, name=f'{self.name}_satV', **constants)
        op_value_sat_vapor = getattr(satd_vapor_state, op)
        op_value_sat_liquid = getattr(satd_liquid_state, op)
        companion_state = CubicEOS(P=self.P, name=f'{self.name}_companion', **constants)
        def objective(T):
            setattr(companion_state, 'T', T)
            return getattr(companion_state, op) - op_value
        if op_value > op_value_sat_vapor: # superheated vapor
            # find T that gives desired property value at P
            # T always > Tsat, so all trials will be in superheated region
            t_bracket = [Tsat.m_as(tunit), self.Tc.m_as(tunit)*5.0]
        else:
            t_bracket = [0.5 * self.Tc.m_as(tunit), Tsat.m_as(tunit)]

        sol = root_scalar(objective, bracket=t_bracket, method='bisect', xtol=self.epsilon, maxiter=self.maxiter)
        if not sol.converged:
            raise ValueError(f'State {self.name}: Unable to converge to solve for T in unsaturated region given P and {op}.')
        setattr(self, 'T', sol.root * tunit)
        for var in self._STATE_VAR_FIELDS - {'T', 'P', 'x', op}:
            setattr(self, var, getattr(companion_state, var))

    def _resolve_saturated_xop(self) -> bool:
        specs = list(self._cache['_input_values'].keys())
        op = specs[0] if specs[1] == 'x' else specs[1]
        op_value = getattr(self, op)
        constants = {k: getattr(self, k) for k in self._PARAMETER_FIELDS }
        companion_state = CubicEOS(x=self.x, name=f'{self.name}_companion', **constants)
        def objective(T):
            setattr(companion_state, 'T', T)
            return getattr(companion_state, op) - op_value
        T_bracket = [0.5 * self.Tc.m_as('K'), self.Tc.m_as('K')]
        sol = root_scalar(objective, bracket=T_bracket, method='bisect', xtol=self.epsilon, maxiter=self.maxiter)
        if not sol.converged:
            raise ValueError(f'State {self.name}: Unable to converge to solve for T in saturated region given x and {op}.')
        setattr(self, 'T', sol.root * _ureg.kelvin)
        self.Liquid = CubicEOS(x=0.0, T=self.T, name=f'{self.name}_L' if self.name else 'Saturated Liquid', **constants)
        self.Vapor = CubicEOS(x=1.0, T=self.T, name=f'{self.name}_V' if self.name else 'Saturated Vapor', **constants)
        for var in self._STATE_VAR_FIELDS - {'T', 'P', 'x', op}:
            setattr(self, var, self.x * getattr(self.Vapor, var) + (1 - self.x) * getattr(self.Liquid, var))
        return True

    def _resolve_op_op(self) -> bool:
        specs = list(self._cache['_input_values'].keys())
        op1, op2 = specs[0], specs[1]
        op1_value = getattr(self, op1)
        op2_value = getattr(self, op2)
        constants = {k: getattr(self, k) for k in self._PARAMETER_FIELDS}
        companion_state = CubicEOS(name=f'{self.name}_companion', T=T, **constants)
        def objective(T):
            setattr(companion_state, op1, op1_value)
            companion_op2_value = getattr(companion_state, op2)
            return companion_op2_value - op2_value
        T_bracket = [0.5 * self.Tc.m_as('K'), self.Tc.m_as('K')*5.0]
        sol = root_scalar(objective, bracket=T_bracket, method='bisect', xtol=self.epsilon, maxiter=self.maxiter)
        if not sol.converged:
            raise ValueError(f'State {self.name}: Unable to converge to solve for T in unsaturated region given {op1} and {op2}.')
        setattr(self, 'T', sol.root * _ureg.kelvin)
        for var in self._STATE_VAR_FIELDS - {'T', op1, op2}:
            companion_var_value = getattr(companion_state, var)
            setattr(self, var, companion_var_value)
        setattr(self, 'Z', companion_state.Z)

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

    # @property
    # def v(self) -> np.ndarray:
    #     """
    #     Computes or retrieves molar volume(s) at current state T and P.
        
    #     Returns
    #     -------
    #     np.ndarray
    #         Molar volume(s) at current state T and P
    #     """
    #     if not 'v' in self._cache:
    #         self._cache['v'] = self.Z * self.R_pv * self.T / self.P
    #         self._input_state = self._get_current_input_state()  # Save state
    #     return self._cache['v']

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

    # def unit_consistency(self, other: CubicEOS):
    #     consistent = self.pressure_unit == other.pressure_unit and self.volume_unit == other.volume_unit
    #     if not consistent:
    #         raise ValueError('inconsistent units')

    # @property
    # def h(self) -> np.ndarray:
    #     """
    #     Computes absolute enthalpy at state T and P
        
    #     Returns
    #     -------
    #     np.ndarray
    #         Absolute enthalpy/enthalpies at state T and P
    #     """
    #     if self.Cp is None:
    #         raise ValueError("Cp data required for absolute enthalpy calculation.")
    #     dH_ideal = DeltaH_IG(self.Tref, self.T, self.Cp)
    #     logger.debug(f'dh_ideal: {dH_ideal} for Tref {self.Tref} to T {self.T} with Cp {self.Cp}')
    #     return self.h_departure + dH_ideal

    # @property
    # def u(self) -> np.ndarray:
    #     """
    #     Computes absolute internal energy at state T and P
        
    #     Returns
    #     -------
    #     np.ndarray
    #         Absolute internal energy/energies at state T and P
    #     """
    #     if self.Cp is None:
    #         raise ValueError("Cp data required for absolute internal energy calculation.")
    #     # u = h - pv
    #     return self.h - self.Pv

    # @property
    # def s(self) -> np.ndarray:
    #     """
    #     Computes absolute entropy at state T and P
        
    #     Returns
    #     -------
    #     np.ndarray
    #         Absolute entropy/entropies at state T and P
    #     """
    #     if self.Cp is None:
    #         raise ValueError("Cp data required for absolute entropy calculation.")
    #     # make sure Pref_local is in correct units (same as self.P)
    #     dS_ideal = DeltaS_IG(self.Tref, self.Pref_local, self.T, self.P, self.Cp, self.R)
    #     return self.s_departure + dS_ideal

    # def delta_h(self, other: CubicEOS) -> np.ndarray:
    #     """
    #     Computes and returns enthalpy change from self to other state
    #     """
    #     self.unit_consistency(other)
    #     return other.h - self.h

    # def delta_s(self, other: CubicEOS) -> np.ndarray:
    #     """
    #     Computes and returns entropy change from self to other state
    #     """
    #     self.unit_consistency(other)
    #     return other.s - self.s
    
    # def delta_pv(self, other: CubicEOS) -> np.ndarray:
    #     """
    #     Returns Delta(PV) in thermal (not PV) units 
    #     """
    #     self.unit_consistency(other)
    #     return (other.Pv - self.Pv) * self.R / self.R_pv
    
    # def delta_u(self, other: CubicEOS) -> np.ndarray:
    #     """
    #     Returns Delta(U) (internal energy)
    #     """
    #     return self.delta_h(other) - self.delta_pv(other)
        
    # def _check_saturation(self, specs: dict[str, float | pint.Quantity]) -> bool:
    #     """
    #     """

    #     if 'x' in specs:
    #         """ A vapor fraction was specified as input """
    #         return True

    #     p, op = None, None
    #     hasT, hasP = 'T' in specs, 'P' in specs
    #     has_only_T_or_P = hasT ^ hasP
    #     if hasT and has_only_T_or_P:
    #         p, v = 'T', self.T
    #         if v > self.Tc:
    #             return False
    #         op = specs[0] if specs[1] == v else specs[1]
    #         pcomp, pcomp_satd = 'P', self.Pvap
    #     elif hasP and has_only_T_or_P:
    #         p, v = 'P', self.P
    #         if v > self.Pc:
    #             return False
    #         pcomp, pcomp_satd = 'T', self.Tsat
    #         op = specs[0] if specs[1] == v else specs[1]

    #     if p is not None and op is not None:
    #         # Either T or P, and one other property from h, u, v, s is specified



    #         if p is None or np.isnan(p):
    #             return False
    #         # check if specified variable equals saturation value within tolerance
    #         tol = 1e-6 * p.m_as(op.units)
    #         return np.abs(op.m_as(op.units) - p.m_as(op.units)) < tol


    #     satd_clone = self.clone()
    #     if satd_clone.T is not None:
    #         satd_clone.P = self.Pvap
    #     elif satd_clone.P is not None:
    #         satd_clone.T = self.Tsat
    #     if v is not None:
    #         vV, vL = satd_clone.v
    #         return (v - vL) / (vV - vL)
    #     elif h is not None:
    #         hV, hL = satd_clone.h
    #         return (h - hL) / (hV - hL)
    #     elif s is not None:
    #         sV, sL = satd_clone.s
    #         return (s - sL) / (sV - sL)
    #     elif u is not None:
    #         uV, uL = satd_clone.u
    #         return (u - uL) / (uV - uL)
    #     else:
    #         raise ValueError("One of v, h, s, or u must be provided to determine vapor fraction.")

    # def P_solve(self, P_bracket: list[float] = None, v: float = None, h: float = None, s: float = None, u: float = None):
    #     """
    #     Solve for pressure given T and one other property (v, h, s, u). Stores result in self.P.

    #     Parameters
    #     ----------
    #     P_bracket: list[float]
    #         Bracket [P_low, P_high] for root finding
    #     v: float
    #         Molar volume
    #     h: float
    #         Enthalpy
    #     s: float
    #         Entropy
    #     u: float
    #         Internal energy
    #     """
    #     if v is not None:
    #         def objective(P):
    #             self.P = P
    #             logger.debug(f'Objective eval: P={P}, v={self.v}, target v={v}')
    #             # pick the root closest to target v
    #             v_diffs = np.abs(self.v - v)
    #             v_closest = self.v[np.argmin(v_diffs)]
    #             return v - v_closest
    #     elif h is not None:
    #         def objective(P):
    #             self.P = P
    #             logger.debug(f'Objective eval: P={P}, h={self.h}, target h={h}')
    #                 # pick the root closest to target h
    #             h_diffs = np.abs(self.h - h)
    #             h_closest = self.h[np.argmin(h_diffs)]
    #             return h - h_closest
    #     elif s is not None:
    #         def objective(P):
    #             self.P = P
    #             logger.debug(f'Objective eval: P={P}, s={self.s}, target s={s}')
    #             # pick the root closest to target s
    #             s_diffs = np.abs(self.s - s)
    #             s_closest = self.s[np.argmin(s_diffs)]
    #             return s - s_closest
    #     elif u is not None:
    #         def objective(P):
    #             self.P = P
    #             logger.debug(f'Objective eval: P={P}, u={self.u}, target u={u}')
    #             # pick the root closest to target u
    #             u_diffs = np.abs(self.u - u)
    #             u_closest = self.u[np.argmin(u_diffs)]
    #             return u - u_closest
    #     else:
    #         raise ValueError("One of v, h, s, or u must be provided to solve for pressure.")
    #     res = root_scalar(objective, bracket=P_bracket)
    #     self.P = res.root

    # def T_solve(self, T_bracket: list[float] = None, v: float = None, h: float = None, s: float = None, u: float = None):
    #     """
    #     Solve for temperature given P and one other property (v, h, s, u). Stores result in self.T.

    #     Parameters
    #     ----------
    #     T_bracket: list[float]
    #         Bracket [T_low, T_high] for root finding
    #     v: float
    #         Molar volume
    #     h: float
    #         Enthalpy
    #     s: float
    #         Entropy
    #     u: float
    #         Internal energy
    #     """
    #     if v is not None:
    #         def objective(T):
    #             self.T = T
    #             logger.debug(f'Objective eval: T={T}, v={self.v}, target v={v}')
    #             # pick the root closest to target v
    #             v_diffs = np.abs(self.v - v)
    #             v_closest = self.v[np.argmin(v_diffs)]
    #             return v - v_closest
    #     elif h is not None:
    #         def objective(T):
    #             self.T = T
    #             logger.debug(f'Objective eval: T={T}, h={self.h}, target h={h}')
    #             # pick the root closest to target h
    #             h_diffs = np.abs(self.h - h)
    #             h_closest = self.h[np.argmin(h_diffs)]
    #             return h - h_closest
    #     elif s is not None:
    #         def objective(T):
    #             self.T = T
    #             logger.debug(f'Objective eval: T={T}, s={self.s}, target s={s}')
    #             # pick the root closest to target s
    #             s_diffs = np.abs(self.s - s)
    #             s_closest = self.s[np.argmin(s_diffs)]
    #             return s - s_closest
    #     elif u is not None:
    #         def objective(T):
    #             self.T = T
    #             logger.debug(f'Objective eval: T={T}, u={self.u}, target u={u}')
    #             # pick the root closest to target u
    #             u_diffs = np.abs(self.u - u)
    #             u_closest = self.u[np.argmin(u_diffs)]
    #             return u - u_closest
    #     else:
    #         raise ValueError("One of v, h, s, or u must be provided to solve for temperature.")

    #     res = root_scalar(objective, bracket=T_bracket)
    #     self.T = res.root

    # def solve(self, T: float = None, P: float = None, v: float = None, h: float = None, s: float = None, u: float = None):
    #     """
    #     Solve for the missing state variable given two of T, P, v, h, s, u.
    #     Updates the object's state variables accordingly.

    #     Parameters
    #     ----------
    #     T: float
    #         Temperature
    #     P: float
    #         Pressure
    #     v: float
    #         Molar volume
    #     h: float
    #         Enthalpy
    #     s: float
    #         Entropy
    #     u: float
    #         Internal energy
    #     """
    #     num_specified = sum(var is not None for var in [T, P, v, h, s, u])
    #     if num_specified != 2:
    #         raise ValueError("Exactly two of T, P, v, h, s, or u must be specified.")
        
    #     self.T = T
    #     self.P = P

    #     # Calculate the missing variable using the EOS
    #     if self.T is not None and self.P is not None:
    #         pass # all properties can be retrieved now
    #     elif self.P is None: # T and one other property (not P) is set
    #         if self.T < self.Tc:
    #             x = self.two_phase_check(v=v, h=h, s=s, u=u)
    #             if 0 < x < 1:
    #                 self.x = x
    #                 self.P = self.Pvap
    #                 # done.  All properties can be retrieved now.
    #                 return
    #             elif x < 0:
    #                 # subcooled liquid
    #                 Plo = self.Pvap
    #                 Phi = self.Pc
    #                 phase = 'liquid'
    #             else:
    #                 # superheated vapor
    #                 Phi = self.Pvap
    #                 Plo = 1.e-5 # guess very low pressure
    #                 phase = 'vapor'
    #             self.phase = phase
    #             self.P_solve([Plo, Phi], v=v, h=h, s=s, u=u)
    #             self.phase = 'unspecified'
    #         else: # given T is above Tc, so supercritical
    #             self.P_solve([1.e-5, self.Pc*5], v=v, h=h, s=s, u=u)
    #     elif self.T is None:
    #         if self.P < self.Pc:
    #             x = self.two_phase_check(v=v, h=h, s=s, u=u)
    #             logger.debug(f'P {self.P} subcrit: x = {x}')
    #             if 0 < x < 1:
    #                 self.x = x
    #                 self.T = self.Tsat
    #                 return
    #             elif x < 0:
    #                 # subcooled liquid
    #                 Tlo = 1.0  # guess low temperature
    #                 Thi = self.Tsat
    #                 phase = 'liquid'
    #             else:
    #                 # superheated vapor
    #                 Thi = self.Tc*3
    #                 Tlo = self.Tsat
    #                 phase = 'vapor'
    #             self.phase = phase
    #             self.T_solve([Tlo, Thi], v=v, h=h, s=s, u=u)
    #             self.phase = 'unspecified'
    #         else:
    #             self.T_solve([10.0, self.Tc*3], v=v, h=h, s=s, u=u)
    #     else:
    #         raise NotImplementedError("Solving when neither T nor P are set is not implemented yet.")

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
            self.Pc = compound.Pc.to('MPa')
            self.Cp = deepcopy(compound.Cp)
        return self