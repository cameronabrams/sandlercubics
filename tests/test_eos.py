from unittest import TestCase
from sandlercubics import IdealGasEOS, PengRobinsonEOS, VanDerWaalsEOS, SoaveRedlichKwongEOS
from sandlerprops.properties import get_database
import numpy as np
import logging

logger = logging.getLogger(__name__)

class TestCubicEOS(TestCase):

    def setUp(self):
        # Load a sample component from the database for testing
        self.methane = get_database().get_compound('methane')
        self.ethane = get_database().get_compound('ethane')
        self.water = get_database().get_compound('water')

    def test_ideal_gas_solve(self):
        eos = IdealGasEOS(pressure_unit="pa", T=300.0, P=101325.0)
        v_calculated = eos.v[0]
        R = 8.314  # J/(mol·K)
        v_expected = R * eos.T / eos.P
        self.assertAlmostEqual(v_calculated, v_expected, places=3)

    def test_ideal_gas_solve_Tv(self):
        P_as_input = 101325.0  # Pa
        eos = IdealGasEOS(pressure_unit="pa", T=300.0, P=P_as_input).set_compound(self.methane)
        v_start = eos.v[0]
        eos.P = None
        eos.solve(T=eos.T, v=v_start)
        P_as_dependent = eos.P
        self.assertAlmostEqual(P_as_input, P_as_dependent, places=3)

    def test_ideal_gas_solve_Pv(self):
        T_as_input = 300.0  # K
        eos = IdealGasEOS(pressure_unit="pa", T=T_as_input, P=101325.0).set_compound(self.methane)
        v_start = eos.v[0]
        eos.T = None
        eos.solve(P=eos.P, v=v_start)
        T_as_dependent = eos.T
        self.assertAlmostEqual(T_as_input, T_as_dependent, places=3)

    def test_ideal_gas_properties(self):
        eos = IdealGasEOS(pressure_unit="pa")
        eos.T = 300.0  # K
        eos.P = 101325.0  # Pa

        # Ideal gas compressibility factor should be 1
        Z = eos.Z[0]
        self.assertAlmostEqual(Z, 1.0, places=6)

        # Ideal gas molar volume should be RT/P
        R = 8.314  # J/(mol·K)
        V_ideal = R * eos.T / eos.P
        v = eos.v[0]
        self.assertAlmostEqual(v, V_ideal, places=4)

        # Fugacity coefficient should be 1 for ideal gas
        phi = eos.phi[0]
        self.assertAlmostEqual(phi, 1.0, places=6)

        self.assertEqual(eos.h_departure[0], 0.0)
        self.assertEqual(eos.s_departure[0], 0.0)
        self.assertNotEqual(eos.Pvap, eos.Pvap)  # Pvap is NaN for ideal gas
        self.assertNotEqual(eos.Tsat, eos.Tsat)  # Tsat is NaN for ideal gas
        self.assertNotEqual(eos.Hvap, eos.Hvap)  # Hvap is NaN for ideal gas
        self.assertNotEqual(eos.Svap, eos.Svap)  # Svap is NaN for ideal gas

    def test_peng_robinson_solve(self):
        eos = PengRobinsonEOS(pressure_unit="bar", T=300.0, P=1.01325).set_compound(self.ethane)
        z = eos.Z[0]
        self.assertAlmostEqual(z, 0.99168, places=3)

        eos.P = 50.0  # bar
        z = eos.Z[0]
        self.assertAlmostEqual(z, 0.196798, places=3)

    def test_peng_robinson_water_saturation(self):
        eos = PengRobinsonEOS(pressure_unit="bar", logiter=True, T=373.15).set_compound(self.water)
        Pvap = eos.Pvap
        self.assertAlmostEqual(Pvap, 0.963, places=2)  # Expect ~1 atm in bar, but PR is not so accurate for water

        eos.P = 2.0
        Tsat = eos.Tsat
        self.assertAlmostEqual(Tsat, 394.413, places=2)

    def test_soave_redlich_kwong_solve(self):
        eos = SoaveRedlichKwongEOS(pressure_unit="bar", T=300.0, P=1.01325).set_compound(self.methane)

        z = eos.Z[0]
        self.assertAlmostEqual(z, 0.9985, places=3)
    
        eos.P = 50.0  # bar
        z = eos.Z[0]
        self.assertAlmostEqual(z, 0.9241, places=3)

    def test_soave_redlich_kwong_water_saturation(self):
        Tc = self.water.Tc
        Pc = self.water.Pc
        omega = self.water.Omega
        eos = SoaveRedlichKwongEOS(Tc=Tc, Pc=Pc, omega=omega, pressure_unit="bar", logiter=True)
        eos.T = 373.15
        
        Pvap = eos.Pvap
        assert abs(Pvap - 0.967) < 0.05  # Expect ~1 atm in bar, but SRK is not so accurate for water

        eos.P = 2.0
        Tsat = eos.Tsat
        assert abs(Tsat - 394.413) < 1.0

    def test_peng_robinson_heat_of_vaporization(self):
        eos = PengRobinsonEOS(pressure_unit="bar", logiter=True, T=373.15).set_compound(self.water)
        
        Hvap = eos.Hvap
        self.assertAlmostEqual(Hvap, 42073, delta=500)  # J/mol, PR is approximate

    def test_peng_robinson_entropy_of_vaporization(self):
        eos = PengRobinsonEOS(pressure_unit="bar", logiter=True, T=373.15).set_compound(self.water)
        
        Svap = eos.Svap
        self.assertAlmostEqual(Svap, 112.8, delta=5)  # J/mol-K, PR is approximate

        # check clausius-clapeyron relation
        Svap_CC = eos.Hvap / eos.T
        self.assertAlmostEqual(Svap, Svap_CC, delta=5)

    def test_peng_robinson_solve_Tv(self):
        P_as_input = 10.0
        eos = PengRobinsonEOS(pressure_unit="bar", T=300.0, P=P_as_input).set_compound(self.methane)
        v_start = eos.v[0]
        eos.P = None
        eos.solve(T=eos.T, v=v_start)
        P_as_dependent = eos.P
        self.assertAlmostEqual(P_as_input, P_as_dependent, places=3)

    def test_peng_robinson_solve_Pv(self):
        T_as_input = 350.0
        eos = PengRobinsonEOS(pressure_unit="bar", T=T_as_input, P=15.0).set_compound(self.ethane)
        v_start = eos.v[0]
        eos.T = None
        eos.solve(P=eos.P, v=v_start)
        T_as_dependent = eos.T
        self.assertAlmostEqual(T_as_input, T_as_dependent, places=3)

    def test_peng_robinson_solve_Th(self):
        P_as_input = 11.0
        eos = PengRobinsonEOS(pressure_unit='bar', T=397, P=P_as_input).set_compound(self.ethane)
        h_sav = eos.h[0]
        eos.P = None
        eos.solve(T=eos.T, h=h_sav)
        P_as_dependent = eos.P
        self.assertAlmostEqual(P_as_input, P_as_dependent, places=3)

    def test_peng_robinson_solve_Ts(self):
        P_as_input = 9.0
        eos = PengRobinsonEOS(pressure_unit='bar', T=350, P=P_as_input).set_compound(self.ethane)
        s_sav = eos.s[0]
        eos.P = None
        eos.solve(T=eos.T, s=s_sav)
        P_as_dependent = eos.P
        self.assertAlmostEqual(P_as_input, P_as_dependent, places=3)

    def test_peng_robinson_solve_Ps(self):
        T_as_input = 350.0
        eos = PengRobinsonEOS(pressure_unit='bar', T=T_as_input, P=15.0).set_compound(self.ethane)
        s_sav = eos.s[0]
        eos.T = None
        eos.solve(P=eos.P, s=s_sav)
        T_as_dependent = eos.T
        self.assertAlmostEqual(T_as_input, T_as_dependent, places=3)

    def test_peng_robinson_solve_Tu(self):
        P_as_input = 9.0
        eos = PengRobinsonEOS(pressure_unit='bar', T=350, P=P_as_input).set_compound(self.ethane)
        u_sav = eos.u[0]
        eos.P = None
        eos.solve(T=eos.T, u=u_sav)
        P_as_dependent = eos.P
        self.assertAlmostEqual(P_as_input, P_as_dependent, places=3)

    def test_peng_robinson_solve_Pu(self):
        T_as_input = 350.0
        eos = PengRobinsonEOS(pressure_unit='bar', T=T_as_input, P=15.0).set_compound(self.ethane)
        u_sav = eos.u[0]
        eos.T = None
        eos.solve(P=eos.P, u=u_sav)
        T_as_dependent = eos.T
        self.assertAlmostEqual(T_as_input, T_as_dependent, places=3)

    def test_vdw_solve(self):
        eos = VanDerWaalsEOS(pressure_unit='bar', T=400, P=8.0).set_compound(self.ethane)
        self.assertAlmostEqual(eos.v[0], 0.00405312, places=6)

    def test_peng_robinson_ethane_pvap(self):
        eos = PengRobinsonEOS(pressure_unit='bar', logiter=True, T=270).set_compound(self.ethane)
        Pvap = eos.Pvap
        self.assertAlmostEqual(Pvap, 22.21845, places=2)  # Expect ~1 atm in bar