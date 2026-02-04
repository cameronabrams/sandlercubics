# Author: Cameron F. Abrams, <cfa22@drexel.edu>
from .eos import CubicEOS
from sandlermisc import ureg, R
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger(__name__)

@dataclass
class IdealGasEOS(CubicEOS):
    """
    Ideal gas class derived from CubicEOS
    """
    name: str = "Ideal Gas Equation of State"
    description: str = "Ideal Gas Equation of State"

    _PARAMETER_FIELDS = frozenset({'Cp'})

    def _calc_a(self):
        """ vdW a parameter for ideal gas is zero """
        return 0.0
    
    def _calc_da_dT(self):
        """ Temperature derivative of vdW a parameter for ideal gas is zero """
        return 0.0

    def _calc_b(self):
        """ vdW b parameter for ideal gas is zero """
        return 0.0
    
    def _calc_P(self):
        """ Calculates pressure from ideal gas EOS """
        v = self.v
        T = self.T
        return R * T / v

    @property
    def cubic_coeff(self) -> np.ndarray:
        """ cubic coefficients for ideal gas """
        return np.array([1.0, 0.0, 0.0, -1.0])

    def _solve_for_Z(self) -> np.ndarray:
        """ Solve for compressibility factor Z; ideal-gas value """
        return np.array([1.0])

    def _calc_h_departure(self) -> float:
        """ Enthalpy departure at state T and P; ideal-gas value """
        return 0.0

    def _calc_s_departure(self) -> float:
        """ Entropy departure at state T and P; ideal-gas value """
        return 0.0
    
    def _calc_logphi(self) -> float:
        """ natural log of fugacity coefficient at state T and P; ideal-gas value """
        return 0.0
    
    def _resolve_Top(self) -> bool:
        """ Resolve state given T and one other state property that is not P """
        punit = self.get_default_unit('P')
        tunit = self.get_default_unit('T')
        vunit = self.get_default_unit('v')
        specs = self._cache['_input_vars_specified']
        op = specs[0] if specs[1] == 'T' else specs[1]
        op_value = getattr(self, op)
        self.Z = 1.0
        if op in 'hu':
            return False
        if op == 'v': # T, v given
            return self._calculate_Phus()
        elif op == 's': # T, s
            temp_State = IdealGasEOS(T=self.T, name=f'{self.name}_temp_resolve', Cp=self.Cp)
            # determine P such that s value from temp_state matches
            from scipy.optimize import bisect
            def objective(P_guess):
                temp_State.P = P_guess
                return getattr(temp_State, op).m - op_value.m
            P_solution = bisect(objective, 1.0 * punit.m, 5000.0 * punit.m)
            self.P = P_solution * punit
            return self._calculate_hus()
        else:
            raise NotImplementedError(f"Cannot resolve IdealGasEOS with T and {op} yet.")
        return True
    
    def _resolve_Pop(self) -> bool:
        """ Resolve state given P and one other state property that is not T """
        tunit = self.get_default_unit('T')
        punit = self.get_default_unit('P')
        vunit = self.get_default_unit('v')
        specs = self._cache['_input_vars_specified']
        op = specs[0] if specs[1] == 'P' else specs[1]
        op_value = getattr(self, op)
        self.Z = 1.0
        if op == 'v': # P, v given
            return self._calculate_Thus()
        else:
            temp_State = IdealGasEOS(P=self.P, name=f'{self.name}_temp_resolve', Cp=self.Cp)
            # determine T such that op value from temp_state matches
            from scipy.optimize import bisect
            def objective(T_guess):
                temp_State.T = T_guess
                return getattr(temp_State, op).m - op_value.m
            T_solution = bisect(objective, 1.0 * tunit.m, 5000.0 * tunit.m)
            self.T = T_solution * tunit
            return self._calculate_hus()
        return True

    def _resolve_saturated_Tx(self) -> bool:
        return False

    def _resolve_saturated_Px(self) -> bool:
        return False

    def _resolve_saturated_xop(self) -> bool:
        return False

    def _resolve_op_op(self) -> bool: # vu vh vs uh sh us hs
        specs = self._cache['_input_vars']
        op1, op2 = specs[0], specs[1]
        # the only combo not solvable for ideal gas is h and u together
        if set([op1, op2]) == set(['h', 'u']):
            return False
        op1_value = getattr(self, op1)
        op2_value = getattr(self, op2)
        constants = {k: getattr(self, k) for k in self._PARAMETER_FIELDS}
        self.Z = 1.0
        if op1 in set(['h', 'u']) or op2 in set(['h', 'u']): # vu vh sh us hs
            if op1 == 'h':
                opx = op2
                def objective(T):
                    dH_ideal = DeltaH_IG(self.Tref, T, self.Cp)
                    return dH_ideal.m - op1_value.m
            elif op2 == 'h':
                opx = op1
                def objective(T):
                    dH_ideal = DeltaH_IG(self.Tref, T, self.Cp)
                    return dH_ideal.m - op2_value.m
            elif op1 == 'u':
                opx = op2
                def objective(T):
                    dH_ideal = DeltaH_IG(self.Tref, T, self.Cp)
                    dU_ideal = dH_ideal - R.to(ureg.J / ureg.mol / ureg.K) * T + R.to(ureg.J / ureg.mol / ureg.K) * self.Tref
                    return dU_ideal.m - op1_value.m
            elif op2 == 'u':
                opx = op1
                def objective(T):
                    dH_ideal = DeltaH_IG(self.Tref, T, self.Cp)
                    dU_ideal = dH_ideal - R.to(ureg.J / ureg.mol / ureg.K) * T + R.to(ureg.J / ureg.mol / ureg.K) * self.Tref
                    return dU_ideal.m - op2_value.m
            T_sol = bisect(objective, 1.0 * ureg.kelvin.m, 5000.0 * ureg.kelvin.m)
            setattr(self, 'T', T_sol * ureg.kelvin)
            self._cache['_input_vars'] = ['T', opx]
            return self._resolve_Top()
        assert set([op1, op2]) == set(['v', 's'])  # vs
        # find T and P such that v and s match
        temp_State = IdealGasEOS(name=f'{self.name}_temp_resolve', Cp=self.Cp)
        def objective(vars):
            T_guess, P_guess = vars
            temp_State.T = T_guess
            temp_State.P = P_guess
            err1 = getattr(temp_State, 'v').m - op1_value.m
            err2 = getattr(temp_State, 's').m - op2_value.m
            return [err1, err2]
        from scipy.optimize import fsolve
        T_sol, P_sol = fsolve(objective, [300.0 * ureg.kelvin.m, 101325.0 * ureg.pascal.m])
        self.T = T_sol * ureg.kelvin
        self.P = P_sol * ureg.pascal
        self._cache['_input_vars'] = ['T', 'P']
        return self._calculate_vhus()

    @property
    def Pvap(self) -> ureg.Quantity:
        """ Vapor pressure is undefined for ideal gas """
        return ureg.Quantity(np.nan, self.get_default_unit('P'))

    @property
    def Hvap(self) -> ureg.Quantity:
        """ Enthalpy of vaporization is undefined for ideal gas """
        return ureg.Quantity(np.nan, self.get_default_unit('h'))

    @property
    def Svap(self) -> ureg.Quantity:
        """ Entropy of vaporization is undefined for ideal gas """
        return ureg.Quantity(np.nan, self.get_default_unit('s'))

    @property
    def Tsat(self) -> ureg.Quantity:
        """ Saturation temperature is undefined for ideal gas """
        return ureg.Quantity(np.nan, self.get_default_unit('T'))