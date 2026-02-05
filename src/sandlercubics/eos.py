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
from scipy.optimize import root_scalar, brentq

from sandlermisc.thermals import DeltaH_IG, DeltaS_IG
from sandlermisc import ThermodynamicState, ureg, R, cached_property
from sandlerprops.compound import Compound
from sandlerprops.properties import get_database

logger = logging.getLogger(__name__)

sqrt_2 = np.sqrt(2)

def vapor_pressure_ambrose_walton(T, Tc, Pc, omega):
    Tr = T / Tc
    tau = 1 - Tr
    if tau <= 0:
        raise ValueError(f'vapor_pressure_ambrose_walton: T={T} is above critical temperature Tc={Tc}.')
    f0 = (-5.97616*tau + 1.29874*tau**1.5 - 0.60394*tau**2.5 
          - 1.06841*tau**5) / Tr
    
    f1 = (-5.03365*tau + 1.11505*tau**1.5 - 5.41217*tau**2.5 
          - 7.46628*tau**5) / Tr
    
    ln_Pr = f0 + omega * f1
    logger.debug(f'vapor_pressure_ambrose_walton: T={T}, Tc={Tc}, Pc={Pc}, omega={omega} => Tr {Tr}, tau {tau}, f0 {f0}, f1 {f1}, ln_Pr={ln_Pr}')
    return Pc * np.exp(ln_Pr)

@dataclass
class CubicEOS(ThermodynamicState):
    """
    Abstract class for all Cubic equations of state.
    """
    phase: str = 'unspecified' 
    """ strict phase requirement flag; 'vapor', 'liquid', 'unspecified' (multiple roots) """

    name: str = 'cubic'

    Pc: pint.Quantity = None
    """ critical pressure """
    Tc: pint.Quantity = None
    """ critical temperature """
    omega: float = None
    """ acentricity factor, dimensionless """

    Cp: float | list[float] | dict [str, float] | None = None
    """ heat capacity data for ideal-gas contributions """

    _PARAMETER_ORDERED_FIELDS = ['Tc', 'Pc', 'omega', 'Cp']
    _PARAMETER_FIELDS = frozenset(_PARAMETER_ORDERED_FIELDS)
    """ Fields that define parameters for the EOS; to be defined in subclasses """
    
    Z: float = None
    """ compressibility factor """

    Tref: pint.Quantity = 298.15 * ureg.kelvin
    """ reference temperature for 'absolute' internal energy, enthalpy, and entropy calculations """
    Pref: pint.Quantity = 0.1 * ureg.megapascal
    """ reference pressure for 'absolute' entropy calculations, in MPa """

    logiter: bool = False
    """ flag for logging iterations in phase calculations """
    maxiter: int = 100
    """ maximum iterations in phase calculations """
    epsilon: float = 1.e-5
    """ fugacity tolerance in phase calculations (Fig. 7.5-1 in Sandler 5th ed) """

    @abstractmethod
    def _calc_P(self):
        """
        Calculate pressure from the equation of state.  This is
        equation-of-state specific and must be implemented in subclasses.
        """
        pass

    def _spawn_helper(self) -> CubicEOS:
        """
        Create a helper CubicEOS instance for internal calculations
        """
        # logger.debug(f'_spawn_helper: Creating helper state for {self.name}')
        # logger.debug(f'_spawn_helper: Creating helper state for {self.name} with T={self.T:.3f}, P={self.P:.3f}, Tc={self.Tc:.3f}, Pc={self.Pc:.3f}, omega={self.omega:.3f}, Cp={self.Cp}')
        helper_state = self.__class__.simple(
            name=f'{self.name}_helper' if self.name else 'CubicEOS_helper',
            T=self.T,
            P=self.P,
            Tc=self.Tc,
            Pc=self.Pc,
            omega=self.omega,
            Cp=self.Cp,
            logiter=self.logiter,
            maxiter=self.maxiter,
            epsilon=self.epsilon
        )
        # logger.debug(f'_spawn_helper: Created helper state {helper_state.name} with T={helper_state.T:.3f}, P={helper_state.P:.3f}, Tc={helper_state.Tc:.3f}, Pc={helper_state.Pc:.3f}, omega={helper_state.omega:.3f}, Cp={helper_state.Cp}')
        return helper_state

    def get_default_unit(self, field_name: str) -> pint.Unit:
        """
        Get the default unit for a given field
        """
        _default_unit_map = {
            'P': ureg.megapascal,
            'T': ureg.kelvin,
            'v': ureg.meter**3 / ureg.mol,
            'u': ureg.joule / ureg.mol,
            'h': ureg.joule / ureg.mol,
            's': ureg.joule / (ureg.mol * ureg.kelvin),
        }
        return _default_unit_map.get(field_name, ureg.dimensionless)

    def resolve(self) -> bool:
        """
        Resolve all state variables from the two input variables.   This method
        is called automatically when any property assignment results
        in a fully specified state.  There should never be a situation in which
        this method is called explicitly; it is only called by __setattr__ or __post_init__.
        """
        hasT, hasP, hasx = self._is_specified_input_var('T'), self._is_specified_input_var('P'), self._is_specified_input_var('x')
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
            resolved = self._resolve_saturated_opx()
        else: 
            resolved = self._resolve_op_op()
        if not resolved:
            raise ValueError(f'State {self.name}: Unable to resolve state with inputs: {statevars_speced}')
        return resolved

    def _resolve_TP(self) -> bool:
        """
        Resolve state given T and P.
        """
        logger.debug(f'{self.__class__.__name__}._resolve_TP: Resolving state for T={self.T:.3f}, P={self.P:.3f}')
        Zlist = self._solve_for_Z()
        # logger.debug(f' -> Z roots found: {Zlist}')
        if len(Zlist) > 1:
            if self.T < self.Tc:
                logger.debug(f'****** Pvap calculation in _resolve_TP triggered T={self.T:.3f} < Tc={self.Tc:.3f} ******')
                Pvap = self.Pvap
                logger.debug(f'****** Pvap calculation in _resolve_TP completed {Pvap:.3f} ******')
                if self.P < Pvap:
                    setattr(self, 'Z', Zlist[0])  # vapor root
                elif self.P >= Pvap:
                    setattr(self, 'Z', Zlist[1])   # liquid root
            elif self.P < self.Pc:
                logger.debug(f'****** Tsat calculation in _resolve_TP triggered P={self.P:.3f} < Pc={self.Pc:.3f} ******')
                Tsat = self.Tsat
                logger.debug(f'****** Tsat calculation in _resolve_TP completed {Tsat:.3f} ******')
                if self.T < Tsat:
                    setattr(self, 'Z', Zlist[1])   # liquid root
                elif self.T >= Tsat:
                    setattr(self, 'Z', Zlist[0])  # vapor root
        else:
            setattr(self, 'Z', Zlist[0])
        logger.debug(f'{self.__class__.__name__}._resolve_TP: T={self.T:.3f}, P={self.P:.3f} -> Z={self.Z:.6f}')
        return self._calculate_vhus()


    def _calculate_op(self, op: str) -> bool:
        """
        Calculate other properties (v, h, u, s) from Z, T, P
        """
        vunit = self.get_default_unit('v')
        punit = self.get_default_unit('P')
        tunit = self.get_default_unit('T')
        eunit = self.get_default_unit('h')
        sunit = self.get_default_unit('s')
        if 'v' in op:
            v = self.Z * R.to(punit * vunit / tunit) * self.T / self.P
            setattr(self, 'v', v)
        if 'h' in op or 'u' in op or 's' in op:
            if self.Cp is None:
                raise ValueError("Cp data required for absolute enthalpy/entropy calculation.")
        if 'h' in op:
            dH_ideal = DeltaH_IG(self.Tref.to('K'), self.T.to('K'), self.Cp).to(eunit)
            h = self.h_departure + dH_ideal
            setattr(self, 'h', h)
            if 'u' in op:
                Pv = (self.P * self.v).to(eunit)
                u = h - Pv
                setattr(self, 'u', u)
                setattr(self, 'Pv', Pv)
        if 'u' in op and not 'h' in op:
            dH_ideal = DeltaH_IG(self.Tref.to('K'), self.T.to('K'), self.Cp).to(eunit)
            Pv = (self.P * self.v).to(eunit)
            u = h - Pv
            setattr(self, 'u', u)
            setattr(self, 'Pv', Pv)
        if 's' in op:
            dS_ideal = DeltaS_IG(self.Tref.to('K'), self.Pref.to(punit), 
                                self.T.to('K'), self.P.to(punit), 
                                self.Cp, R.to(eunit/tunit)).to(sunit)
            s = self.s_departure + dS_ideal
            setattr(self, 's', s)
        return True

    def _calculate_hus(self) -> bool:
        return self._calculate_op('hus')

    def _calculate_vhus(self) -> bool: # T, P
        return self._calculate_op('vhus')

    def _calculate_Thus(self) -> bool: # P, v
        punit = self.get_default_unit('P')
        vunit = self.get_default_unit('v')
        tunit = self.get_default_unit('T')
        T = self.P * self.v / R.to(punit * vunit / tunit) / self.Z
        setattr(self, 'T', T)
        return self._calculate_hus()

    def _calculate_Phus(self) -> bool: # T, v
        punit = self.get_default_unit('P')
        vunit = self.get_default_unit('v')
        tunit = self.get_default_unit('T')
        P = self.Z * R.to(punit * vunit / tunit) * self.T / self.v
        setattr(self, 'P', P)
        return self._calculate_hus()

    def _resolve_saturated_Tx(self) -> bool:
        if self.T > self.Tc:
            raise ValueError(f'State {self.name}: Cannot have saturated state at temperature T={self.T} above critical temperature Tc={self.Tc}.')
        self.P = self.Pvap
        if 0 < self.x < 1:
            self.Liquid = self.__class__(x=0.0, T=self.T, phase='liquid',
                name=f'{self.name}_L' if self.name else 'Saturated Liquid',
                Tc=self.Tc, Pc=self.Pc, omega=self.omega, Cp=self.Cp)
            self.Vapor = self.__class__(x=1.0, T=self.T, phase='vapor',
                name=f'{self.name}_V' if self.name else 'Saturated Vapor',
                Tc=self.Tc, Pc=self.Pc, omega=self.omega, Cp=self.Cp)
            for op in self._STATE_VAR_FIELDS - {'T', 'P', 'x'}:
                setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
            return True
        else:
            return self._resolve_TP()

    def _resolve_saturated_Px(self) -> bool:
        if self.P > self.Pc:
            raise ValueError(f'State {self.name}: Cannot have saturated state at pressure P={self.P} above critical pressure Pc={self.Pc}.')
        self.T = self.Tsat
        if 0 < self.x < 1:
            self.Liquid = self.__class__(x=0.0, P=self.P, phase='liquid', name=f'{self.name}_L' if self.name else 'Saturated Liquid', Tc=self.Tc, Pc=self.Pc, omega=self.omega, Cp=self.Cp)
            self.Vapor = self.__class__(x=1.0, P=self.P, phase='vapor', name=f'{self.name}_V' if self.name else 'Saturated Vapor', Tc=self.Tc, Pc=self.Pc, omega=self.omega, Cp=self.Cp)
            for op in self._STATE_VAR_FIELDS - {'T', 'P', 'x'}:
                setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
            return True
        else:
            return self._resolve_TP()

    def _resolve_Top(self) -> bool:
        """
        Resolve state given T and one other property.
        """
        punit = self.get_default_unit('P')
        specs = self.get_input_varnames()
        op = specs[0] if specs[1] == 'T' else specs[1]
        if op == 'v':
            self.P = self._calc_P()
            self.Z = self.P * self.v / (R * self.T)
            self._calculate_hus()
            return True
        bounds = (0.01 * self.Pc.m_as(punit), self.Pc.m_as(punit)*5.0)
        op_value_target = getattr(self, op)
        self.swap_input_vars(op, 'P')
        def residual(P_guess):
            self.P = P_guess
            self._resolve_TP()
            op_value_current = getattr(self, op)
            return op_value_current.m - op_value_target.m
        r1 = residual(bounds[0])
        r2 = residual(bounds[1])
        logger.debug(f'{self.__class__.__name__}._resolve_Top: residual at bound1 {bounds[0]:.3f}: {r1:.3f}, residual at bound2 {bounds[1]:.3f}: {r2:.3f}')
        niter = 0
        bail = False
        while r1 * r2 > 0:
            bounds = (bounds[0] + 0.1, bounds[1])
            r1 = residual(bounds[0])
            r2 = residual(bounds[1])
            logger.debug(f'{self.__class__.__name__}._resolve_Top: iter {niter}, residual at bound1 {bounds[0]:.3f}: {r1:.3f}, residual at bound2 {bounds[1]:.3f}: {r2:.3f}')
            niter += 1
            if niter > 20:
                logger.debug(f'State {self.name}: Unable to bracket root for {op}={op_value_target:.3f} after {niter} iterations.')
                bail = True
                break
        if not bail:
            P_solution = brentq(residual, bounds[0], bounds[1], xtol=1e-6)
            self.P = P_solution
            self._resolve_TP()
        self.swap_input_vars('P', op)
        if bail:
            return False
        return True

    def _resolve_Pop(self) -> bool:
        tunit = self.get_default_unit('T')
        specs = self.get_input_varnames()
        op = specs[0] if specs[1] == 'P' else specs[1]
        bounds = (0.1 * self.Tc.m_as(tunit), self.Tc.m_as(tunit)*3.0)
        op_value_target = getattr(self, op)
        self.swap_input_vars(op, 'T')
        def residual(T_guess):
            self.T = T_guess
            self._resolve_TP()
            op_value_current = getattr(self, op)
            return op_value_current.m - op_value_target.m
        T_solution = brentq(residual, bounds[0], bounds[1], xtol=1e-6)
        self.T = T_solution
        self._resolve_TP()
        self.swap_input_vars('T', op)
        return True

    def _resolve_saturated_opx(self) -> bool:
        tunit = self.get_default_unit('T')
        specs = self.get_input_varnames()
        op = specs[0] if specs[1] == 'x' else specs[1]
        op_value = getattr(self, op)
        self.swap_input_vars(op, 'T')
        bounds = (self.Tc.m_as(tunit)*0.3, self.Tc.m_as(tunit)*0.95)
        def residual(T_guess):
            self.T = T_guess
            self._resolve_saturated_Tx()
            op_value_current = getattr(self, op)
            logger.debug(f'_resolve_saturated_opx: T_guess={T_guess:.3f}, op={op}, op_value_current={op_value_current:.3f}, op_value_target={op_value:.3f}')
            return op_value_current.m - op_value.m
        if op == 's':
            r1 = residual(bounds[0])
            r2 = residual(bounds[1])
            niter = 0
            while r1 * r2 > 0:
                bounds = (bounds[0] + 10.0, bounds[1])
                r1 = residual(bounds[0])
                r2 = residual(bounds[1])
                niter += 1
                if niter > 20:
                    raise ValueError(f'State {self.name}: Unable to bracket root for saturated {op}={op_value:.3f} after {niter} iterations.')
            logger.debug(f'_resolve_saturated_opx: Needed {niter} expansions of bounds to find sign change.')
        T_solution = brentq(residual, bounds[0], bounds[1], xtol=1e-6)
        self.T = T_solution
        self._resolve_saturated_Tx()
        self.swap_input_vars('T', op)
        return True

    def _resolve_op_op(self) -> bool:
        raise NotImplementedError('_resolve_op_op not implemented yet')
        specs = self.get_input_varnames()
        op1, op2 = specs[0], specs[1]
        op1_value = getattr(self, op1)
        op2_value = getattr(self, op2)
        logger.debug(f'{self.__class__.__name__}._resolve_op_op: {op1} = {op1_value:.3f}, {op2} = {op2_value:.3f}')

        self.swap_input_vars(op1, 'T')
        bounds = (0.25 * self.Tc.m_as(self.get_default_unit('T')), self.Tc.m_as(self.get_default_unit('T'))*2.0)
        def residual(T_guess):
            self.T = T_guess
            niter = 0
            while not self._resolve_Top():
                logger.debug(f'{self.__class__.__name__}._resolve_op_op: _resolve_Top failed at T={self.T:.3f}, trying next T guess.')
                T_guess += 5.0
                self.T = T_guess
                niter += 1
                if niter > 50:
                    raise ValueError(f'State {self.name}: Unable to resolve state for {op1}={op1_value:.3f} after {niter} attempts adjusting T.')
            op2_value_current = getattr(self, op2)
            return op2_value_current.m - op2_value.m
        r1 = residual(bounds[0])
        r2 = residual(bounds[1])
        niter = 0
        while r1 * r2 > 0:
            bounds = (bounds[0] + 10.0, bounds[1])
            r1 = residual(bounds[0])
            r2 = residual(bounds[1])
            niter += 1
            if niter > 20:
                raise ValueError(f'State {self.name}: Unable to bracket root for saturated {op}={op_value:.3f} after {niter} iterations.')
        logger.debug(f'_resolve_saturated_opx: Needed {niter} expansions of bounds to find sign change.')
        T_solution = brentq(residual, bounds[0], bounds[1], xtol=1e-6)
        self.T = T_solution
        self._resolve_Top()
        self.swap_input_vars('T', op1)
        logger.debug(f'{self.__class__.__name__}._resolve_op_op: Resolved T = {self.T:.3f}, P = {self.P:.3f}')
        return True

    @abstractmethod
    def _calc_a(self):
        """
        Calculate the van der Waals a parameter.  This is equation-of-state specific
        and must be implemented in subclasses.
        """
        pass

    @cached_property
    def a(self):
        """
        van der Waals a parameter, and a cached property.  This method should not need to be further
        overridden in subclasses.
        """
        pass
    
    @abstractmethod
    def _calc_da_dT(self):
        """
        Returns the derivative of vdW a parameter wrt temperature.  This is
        equation-of-state specific and must be implemented in subclasses.
        """
        pass

    @cached_property
    def da_dT(self):
        pass

    @abstractmethod
    def _calc_b(self):
        pass

    @cached_property
    def b(self):
        pass

    def _calc_A(self):
        return self.a * self.P / (R * self.T)**2

    @cached_property
    def A(self):
        pass

    def _calc_B(self):
        return self.b * self.P / (R * self.T)

    @cached_property
    def B(self):
        pass

    @property
    @abstractmethod
    def cubic_coeff(self):
        """
        Coefficients of cubic equation for compressibility factor Z. These are 
        equation-of-state specific and must be implemented in subclasses.
        """
        pass

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
                # logger.debug(f'{self.__class__.__name__}._solve_for_Z: Multiple roots found: {real_roots}, sorting and returning vapor and liquid roots if present')
                positive_roots = real_roots[real_roots > 0]
                if len(positive_roots) < 2:
                    Z = positive_roots
                else:
                    # sorted_roots = np.sort(positive_roots)
                    Z = np.array([real_roots[0], real_roots[-1]])
        return Z

    @abstractmethod
    def _calc_h_departure(self) -> float:
        """
        Calculate the enthalpy departure at state T and P.  This is equation-of-state specific 
        and must be implemented in subclasses.
        """
        pass

    @cached_property
    def h_departure(self) -> float:
        pass

    @abstractmethod
    def _calc_s_departure(self) -> float:
        """
        Calculate the entropy departure at state T and P.  This is equation-of-state specific 
        and must be implemented in subclasses.
        """
        pass

    @cached_property
    def s_departure(self) -> float:
        pass
    
    @abstractmethod
    def _calc_logphi(self) -> float:
        """
        Calculate the natural log of fugacity coefficient at state T and P.
        This is equation-of-state specific and must be implemented in subclasses.
        """
        pass

    @cached_property
    def logphi(self) -> float:
        pass
    
    def _calc_phi(self) -> float:
        """
        Calculate fugacity coefficient(s) at current state T and P.
        """
        return np.exp(self.logphi)

    @cached_property
    def phi(self) -> float:
        pass

    def _calc_f(self) -> float:
        """
        Calculate fugacity/fugacities at current state T and P.
        """
        return self.phi * self.P

    @cached_property
    def f(self) -> float:
        pass

    @cached_property
    def Pvap(self) -> float:
        pass

    def z_helpers(self, hV: CubicEOS, hL: CubicEOS):
        Zlist = hV._solve_for_Z()
        # logger.debug(f'z_helpers: Z roots found: {Zlist}')
        ZV, ZL = Zlist[0], Zlist[1]
        hV.Z = ZV
        hL.Z = ZL
        # logger.debug(f'z_helpers: ZV {ZV:.6f}, ZL {ZL:.6f}')
        # fV = hV.f
        # logger.debug(f'z_helpers: fV {fV:.6f}')
        # fL = hL.f
        # logger.debug(f'z_helpers: fL {fL:.6f}')

    def _calc_Pvap(self) -> float:
        if self.T >= self.Tc:
            raise ValueError(f'State {self.name}: Cannot calculate Pvap at T={self.T} above critical temperature Tc={self.Tc}.')
        hV = self._spawn_helper()
        hL = self._spawn_helper()
        hV.P = (vapor_pressure_ambrose_walton(self.T.m_as('K'), self.Tc.m_as('K'), self.Pc.m_as('Pa'), self.omega) * ureg.pascal).to(self.get_default_unit('P'))
        logger.debug(f'{self.__class__.__name__}._calc_Pvap: Initial Pvap guess: {hV.P:.3f} at T {hV.T:.3f} (Pc {hV.Pc:.3f}, Tc {hV.Tc:.3f})')
        hV.T = self.T
        Zlist = hV._solve_for_Z()
        niter = 0
        while len(Zlist) < 2:
            z = Zlist[0]
            if z > 0.27:
                logger.debug(f'{self.__class__.__name__}._calc_Pvap: Pvap at {self.T:.3f} (Tc={self.Tc:.3f}) initial guess too low; z = {z} increasing P guess from {hV.P:.3f} to {hV.P*1.05:.3f}')
                hV.P *= 1.05
            else:
                logger.debug(f'{self.__class__.__name__}._calc_Pvap: {self.__class__.__name__}._calc_Pvap: Pvap at {self.T:.3f} (Tc={self.Tc:.3f}) initial guess too high; z = {z} decreasing P guess from {hV.P:.3f} to {hV.P*0.95:.3f}')
                hV.P *= 0.95
            Zlist = hV._solve_for_Z()
            niter += 1
            if niter > 50:
                raise ValueError(f'State {self.name}: Unable to find initial Pvap guess after {niter} adjustments.')
        hL.P = hV.P
        hL.T = hV.T
        logger.debug(f'{self.__class__.__name__}._calc_Pvap: Initial Pvap guess: {hV.P:.3f} at T {hV.T:.3f} (Pc {hV.Pc:.3f}, Tc {hV.Tc:.3f})')
        keepgoing = True
        i = 0
        while keepgoing:
            self.z_helpers(hV, hL)
            i += 1
            try:
                fV, fL = hV.f, hL.f
            except:
                return np.nan
            err = np.abs(fL / fV - 1)
            if hV.logiter: logger.debug(f'Pvap: Iter {i}: P {hV.P:.3f}, fV {fV:.3f}, fL {fL:.3f}; error {err:.4e}')
            hV.P *= fL / fV
            hL.P = hV.P
            if err < hV.epsilon or i == hV.maxiter:
                keepgoing = False
            if i >= hV.maxiter:
                logger.warning(f'{self.__class__.__name__}._calc_Pvap: Reached {i} iterations without convergence; error {err:.4e}')
        return hV.P
    
    def _calc_Hvap(self) -> float:
        hV = self._spawn_helper()
        hV.T = self.T
        hV.P = self.Pvap
        hL = self._spawn_helper()
        hL.T = self.T
        hL.P = hV.P
        self.z_helpers(hV, hL)
        return hV.h_departure - hL.h_departure

    @cached_property
    def Hvap(self) -> float:
        pass
    
    def _calc_Svap(self) -> float:
        hV = self._spawn_helper()
        hV.T = self.T
        hV.P = self.Pvap
        hL = self._spawn_helper()
        hL.T = self.T
        hL.P = hV.P
        self.z_helpers(hV, hL)
        return hV.s_departure - hL.s_departure

    @cached_property
    def Svap(self) -> float:
        pass

    def _calc_Tsat(self) -> float:
        if self.P >= self.Pc:
            raise ValueError(f'State {self.name}: Cannot calculate Tsat at P={self.P} above critical pressure Pc={self.Pc}.')
        hV = self._spawn_helper()
        hL = self._spawn_helper()
        hV.P = self.P
        hV.T = hV.Tc * (hV.P / hV.Pc)**0.125
        Zlist = hV._solve_for_Z()
        niter = 0
        while len(Zlist) < 2:
            z = Zlist[0]
            if z > 0.27:
                logger.debug(f'{self.__class__.__name__}._calc_Tsat: Tsat at {self.P:.3f} (Pc={self.Pc:.3f}) initial guess too low; z = {z} increasing T guess from {hV.T:.3f} to {hV.T*1.05:.3f}')
                hV.T *= 1.05
            else:
                logger.debug(f'{self.__class__.__name__}._calc_Tsat: Tsat at {self.P:.3f} (Pc={self.Pc:.3f}) initial guess too high; z = {z} decreasing T guess from {hV.T:.3f} to {hV.T*0.95:.3f}')
                hV.T *= 0.95
            Zlist = hV._solve_for_Z()
            niter += 1
            if niter > 50:
                raise ValueError(f'{self.__class__.__name__}._calc_Tsat: Unable to find initial Tsat guess after {niter} adjustments.')
        hL.P = hV.P
        hL.T = hV.T
        logger.debug(f'{self.__class__.__name__}._calc_Tsat: Initial Tsat guess: {hL.T:.3f} at {hL.P:.3f} (Pc {hL.Pc:.3f}, Tc {hL.Tc:.3f})')
        keepgoing = True
        i = 0
        while keepgoing:
            self.z_helpers(hV, hL)
            i += 1
            try:
                fV, fL = hV.f, hL.f
            except:
                # raise ValueError(f'Error in fugacity calculation in Tsat at T {hL.T:.2f}, P {hL.P:.2f} {pu}')
                return np.nan
            err = np.abs(fV / fL - 1)
            if hL.logiter: logger.debug(f'{self.__class__.__name__}._calc_Tsat: Iter {i}: T {hL.T:.6f}, fV {fV:.6f}, fL {fL:.6f}; error {err:.4e}')
            hV.T *= (fV / fL)**0.125
            hL.T = hV.T
            if err < hL.epsilon or i == hL.maxiter:
                keepgoing = False
            if i == hL.maxiter:
                logger.warning(f'{self.__class__.__name__}._calc_Tsat: Reached {i} iterations without convergence; error {err:.4e}')
        return hV.T

    @cached_property
    def Tsat(self) -> float:
        pass

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
            self.omega = compound.Omega
            self.Cp = deepcopy(compound.Cp)
        return self
