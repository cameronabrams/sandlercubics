from unittest import TestCase
from sandlercubics import IdealGasEOS
from sandlerprops.properties import Compound, get_database
from sandlermisc import R, ureg
import numpy as np
import logging
import pint

logger = logging.getLogger(__name__)

class TestIdealGasEOS(TestCase):

    def test_idealgas_instance(self):
        eos = IdealGasEOS(name='idealgas_test')
        self.assertIsInstance(eos, IdealGasEOS)
        self.assertIsNone(eos.P)
        eos.P = 1.0
        # default units are MPa
        self.assertEqual(eos.P.dimensionality, ureg('MPa').dimensionality)

    def test_idealgas_check_spec(self):
        eos = IdealGasEOS(name='idealgas_test')
        self.assertFalse(eos._cache['_is_specified'])
        self.assertFalse(eos._cache['_is_complete'])
        eos.T = 300
        self.assertEqual(eos.T.dimensionality, ureg('K').dimensionality)
        self.assertFalse(eos._cache['_is_specified'])
        self.assertFalse(eos._cache['_is_complete'])
        eos.P = 5.0
        self.assertTrue(eos._cache['_is_specified'])
        logger.debug(f'Test: EOS {eos.name} is now fully specified, checking completeness.')
        logger.debug(f'cache: {eos._cache}')
        self.assertFalse(eos._cache['_is_complete'])
        eos.Cp = 50.0  # J/mol-K
        self.assertTrue(eos._cache['_is_specified'])
        self.assertTrue(eos._cache['_is_complete'])
        self.assertIsNotNone(eos.h)
        self.assertIsNotNone(eos.u)
        self.assertIsNotNone(eos.v)
        self.assertIsNotNone(eos.s)
        self.assertIsNone(eos.x)

    def test_idealgas_v(self):
        eos = IdealGasEOS(name='idealgas_test')
        eos.T = 300.0 * ureg.kelvin
        eos.P = 101325.0 * ureg.pascal # this is converted to MPa internally
        eos.Cp = 50.0  # J/mol-K
        v_calculated = eos.v.m
        v_expected = (R.to(ureg.megapascal * ureg.meter**3 / ureg.mol / ureg.kelvin) * eos.T / eos.P).m
        self.assertAlmostEqual(v_calculated, v_expected, places=3)

    def test_idealgas_solve_Tv(self):
        P_as_input = 101325.0 * ureg.pascal
        mock_eos = IdealGasEOS(T=300.0 * ureg.kelvin, P=P_as_input, Cp=50.0).set_compound('methane')
        solve_eos = IdealGasEOS(T=300.0 * ureg.kelvin, v=mock_eos.v, Cp=50.0).set_compound('methane')
        P_as_dependent = solve_eos.P
        self.assertAlmostEqual(P_as_input, P_as_dependent, places=3)

    def test_idealgas_solve_Pv(self):
        T_as_input = 300.0  * ureg.kelvin
        mock_eos = IdealGasEOS(T=T_as_input, P=101325.0*ureg.pascal, Cp=50.0).set_compound('methane')
        solve_eos = IdealGasEOS(P=mock_eos.P, v=mock_eos.v, Cp=50.0).set_compound('methane')
        T_as_dependent = solve_eos.T
        self.assertAlmostEqual(T_as_input, T_as_dependent, places=3)

    def test_idealgas_properties(self):
        eos = IdealGasEOS()
        eos.T = 300.0  * ureg.kelvin
        eos.P = 101325.0 * ureg.pascal
        eos.Cp = 50.0  # J/mol-K
        # Ideal gas compressibility factor should be 1
        self.assertAlmostEqual(eos.Z, 1.0, places=6)

        # Ideal gas molar volume should be RT/P
        V_ideal = R * eos.T / eos.P
        self.assertAlmostEqual(eos.v, V_ideal, places=4)

        # Fugacity coefficient should be 1 for ideal gas
        self.assertAlmostEqual(eos.phi, 1.0, places=6)

        self.assertEqual(eos.h_departure, 0.0)
        self.assertEqual(eos.s_departure, 0.0)
        self.assertNotEqual(eos.Pvap, eos.Pvap)  # Pvap is NaN for ideal gas
        self.assertNotEqual(eos.Tsat, eos.Tsat)  # Tsat is NaN for ideal gas
        self.assertNotEqual(eos.Hvap, eos.Hvap)  # Hvap is NaN for ideal gas
        self.assertNotEqual(eos.Svap, eos.Svap)  # Svap is NaN for ideal gas

