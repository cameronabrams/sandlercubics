from unittest import TestCase
from sandlercubics import PengRobinsonEOS
from sandlerprops.properties import Compound, get_database
from sandlermisc import R, ureg
import numpy as np
import logging
import pint

logger = logging.getLogger(__name__)

logging.getLogger('sandlermisc.thermodynamicstate').setLevel(logging.ERROR)

class TestPengRobinsonEOS(TestCase):
    
    def test_pengrobinson_instance(self):
        eos = PengRobinsonEOS(name='pengrobinson_test')
        self.assertIsInstance(eos, PengRobinsonEOS)
        self.assertIsNone(eos.P)
        eos.P = 1.0
        # default units are MPa
        self.assertEqual(eos.P.dimensionality, ureg('MPa').dimensionality)

    def test_pengrobinson_simple(self):
        eos = PengRobinsonEOS(T=300.0).set_compound('methane')
        self.assertTrue(hasattr(eos, '_calc_a'))
        self.assertTrue(hasattr(eos, '_cache'))
        eos_simple = eos._spawn_helper()
        self.assertIsInstance(eos_simple, PengRobinsonEOS)
        self.assertTrue(hasattr(eos_simple, '_calc_a'))
        self.assertEqual(eos_simple._cache, {})

    def test_pengrobinson_Pvap_at_T(self):
        eos = PengRobinsonEOS(logiter=True, T=373.15).set_compound('water')
        self.assertIsNone(eos.P) # we have not set P yet
        Pvap = eos.Pvap # triggers Pvap calculation
        self.assertAlmostEqual(Pvap, 0.0963*ureg.megapascal, places=2)
        self.assertIsNone(eos.P) # we have not set P yet

        eos.set_compound('ethane') # this should invalidate the whole thing
        eos.T = 200.0 # so should this
        Pvap = eos.Pvap # triggers Pvap calculation
        self.assertAlmostEqual(Pvap, 0.21772279*ureg.megapascal, places=2)

    def test_pengrobinson_resolve_TP(self):
        eos = PengRobinsonEOS(T=300.0, P=1.01325*ureg.bar).set_compound('ethane')
        z = eos.Z
        self.assertAlmostEqual(z, 0.99168, places=3)
        self.assertAlmostEqual(eos.v.magnitude, 0.024412399675717563, places=4)  # m3/mol

    def test_pengrobinson_resolve_post_set(self):
        eos = PengRobinsonEOS().set_compound('ethane')
        eos.T = 300.0
        eos.P = 1.01325*ureg.bar
        z = eos.Z
        self.assertAlmostEqual(z, 0.99168, places=3)
        self.assertAlmostEqual(eos.v.magnitude, 0.024412399675717563, places=4)  # m3/mol

    def test_pengrobinson_water_saturation(self):
        eos = PengRobinsonEOS(logiter=True, T=373.15).set_compound('water')
        Pvap = eos.Pvap
        self.assertAlmostEqual(Pvap, 0.0963*ureg.megapascal, places=2)  # Expect ~1 atm in bar, but PR is not so accurate for water
        eos.P = 0.2
        Tsat = eos.Tsat
        self.assertAlmostEqual(Tsat, 394.413*ureg.kelvin, places=2)

    def test_pengrobinson_heat_of_vaporization(self):
        eos = PengRobinsonEOS(logiter=True, T=373.15).set_compound('water')
        Hvap = eos.Hvap
        self.assertAlmostEqual(Hvap, 42073*ureg.joule/ureg.mol, delta=500*ureg.joule/ureg.mol)  # J/mol, PR is approximate

    def test_pengrobinson_entropy_of_vaporization(self):
        eos = PengRobinsonEOS(logiter=True, T=373.15).set_compound('water')
        Svap = eos.Svap
        self.assertAlmostEqual(Svap, 112.8*ureg.joule/(ureg.mol*ureg.kelvin), delta=5*ureg.joule/(ureg.mol*ureg.kelvin))  # J/mol-K, PR is approximate
        # check clausius-clapeyron relation
        Svap_CC = eos.Hvap / eos.T
        self.assertAlmostEqual(Svap, Svap_CC, delta=5*ureg.joule/(ureg.mol*ureg.kelvin))
        
    def test_pengrobinson_resolve_Tv(self):
        P_as_input = 10.0*ureg.bar
        eos = PengRobinsonEOS(T=300.0, P=P_as_input).set_compound('methane')
        v_start = eos.v
        h_start = eos.h
        s_start = eos.s
        u_start = eos.u
        eos = PengRobinsonEOS(T=300.0, v=v_start).set_compound('methane')
        P_as_dependent = eos.P
        h_as_dependent = eos.h
        s_as_dependent = eos.s
        u_as_dependent = eos.u
        self.assertAlmostEqual(P_as_input, P_as_dependent, places=3)
        self.assertAlmostEqual(h_start.magnitude, h_as_dependent.magnitude, places=3)
        self.assertAlmostEqual(s_start.magnitude, s_as_dependent.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, u_as_dependent.magnitude, places=3)

    def test_pengrobinson_resolve_Pv(self):
        T_as_input = 350.0
        eos = PengRobinsonEOS(T=T_as_input, P=1.50).set_compound('ethane')
        v_start = eos.v
        h_start = eos.h
        s_start = eos.s
        u_start = eos.u
        eos = PengRobinsonEOS(P=eos.P, v=v_start).set_compound('ethane')
        P_as_dependent = eos.P
        h_as_dependent = eos.h
        s_as_dependent = eos.s
        u_as_dependent = eos.u
        T_as_dependent = eos.T
        self.assertAlmostEqual(T_as_input*ureg.kelvin, T_as_dependent, places=3)
        self.assertAlmostEqual(h_start.magnitude, h_as_dependent.magnitude, places=3)
        self.assertAlmostEqual(s_start.magnitude, s_as_dependent.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, u_as_dependent.magnitude, places=3)

    def test_pengrobinson_resolve_Th(self):
        P_as_input = 1.1
        eos = PengRobinsonEOS(T=397, P=P_as_input).set_compound('ethane')
        h_start = eos.h
        v_start = eos.v
        u_start = eos.u
        s_start = eos.s
        eos = PengRobinsonEOS(T=eos.T, h=h_start).set_compound('ethane')
        P_as_dependent = eos.P
        self.assertAlmostEqual(P_as_input*ureg.megapascal, P_as_dependent, places=3)
        self.assertAlmostEqual(v_start.magnitude, eos.v.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, eos.u.magnitude, places=3)
        self.assertAlmostEqual(s_start.magnitude, eos.s.magnitude, places=3)

    def test_pengrobinson_resolve_Ph(self):
        T_as_input = 397.0
        eos = PengRobinsonEOS(T=T_as_input, P=1.1).set_compound('ethane')
        h_start = eos.h
        v_start = eos.v
        u_start = eos.u
        s_start = eos.s
        eos = PengRobinsonEOS(P=eos.P, h=h_start).set_compound('ethane')
        T_as_dependent = eos.T
        self.assertAlmostEqual(T_as_input*ureg.kelvin, T_as_dependent, places=3)
        self.assertAlmostEqual(v_start.magnitude, eos.v.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, eos.u.magnitude, places=3)
        self.assertAlmostEqual(s_start.magnitude, eos.s.magnitude, places=3)

    def test_pengrobinson_resolve_Ts(self):
        P_as_input = 9.0
        eos = PengRobinsonEOS(T=350, P=P_as_input).set_compound('ethane')
        h_start = eos.h
        v_start = eos.v
        u_start = eos.u
        s_start = eos.s
        eos = PengRobinsonEOS(T=eos.T, s=s_start).set_compound('ethane')
        P_as_dependent = eos.P
        self.assertAlmostEqual(P_as_input*ureg.megapascal, P_as_dependent, places=3)
        self.assertAlmostEqual(v_start.magnitude, eos.v.magnitude, places=3)
        self.assertAlmostEqual(h_start.magnitude, eos.h.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, eos.u.magnitude, places=3)

    def test_pengrobinson_resolve_Ps(self):
        T_as_input = 350.0
        eos = PengRobinsonEOS(T=T_as_input, P=1.5).set_compound('ethane')
        h_start = eos.h
        v_start = eos.v
        u_start = eos.u
        s_start = eos.s
        eos = PengRobinsonEOS(P=eos.P, s=s_start).set_compound('ethane')
        T_as_dependent = eos.T
        self.assertAlmostEqual(T_as_input*ureg.kelvin, T_as_dependent, places=3)
        self.assertAlmostEqual(v_start.magnitude, eos.v.magnitude, places=3)
        self.assertAlmostEqual(h_start.magnitude, eos.h.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, eos.u.magnitude, places=3)

    def test_pengrobinson_resolve_Tu(self):
        P_as_input = 9.0
        eos = PengRobinsonEOS(T=350, P=P_as_input).set_compound('ethane')
        h_start = eos.h
        v_start = eos.v
        u_start = eos.u
        s_start = eos.s
        eos = PengRobinsonEOS(T=350, u=u_start).set_compound('ethane')
        P_as_dependent = eos.P
        self.assertAlmostEqual(P_as_input*ureg.megapascal, P_as_dependent, places=3)

    def test_pengrobinson_resolve_Pu(self):
        T_as_input = 350.0
        eos = PengRobinsonEOS(T=T_as_input, P=1.5).set_compound('ethane')
        h_start = eos.h
        v_start = eos.v
        u_start = eos.u
        s_start = eos.s
        eos = PengRobinsonEOS(P=eos.P, u=u_start).set_compound('ethane')
        T_as_dependent = eos.T
        self.assertAlmostEqual(T_as_input*ureg.kelvin, T_as_dependent, places=3)

    def test_pengrobinson_invalid_saturation(self):
        eos = PengRobinsonEOS().set_compound('water')
        eos.T = 650.0  # above critical temperature
        with self.assertRaises(ValueError):
            Pvap = eos.Pvap
        eos.P = 30.0  # above critical pressure
        self.assertTrue(eos.P > eos.Pc)
        with self.assertRaises(ValueError):
            Tsat = eos.Tsat

    def test_pengrobinson_resolve_Tx(self):
        test_eos = PengRobinsonEOS(logiter=True, T=200.0).set_compound('ethane')
        Pvap = test_eos.Pvap
        eos = PengRobinsonEOS(T=200.0, x=0.5).set_compound('ethane')
        self.assertAlmostEqual(eos.P.magnitude, Pvap.magnitude, places=3)
        self.assertAlmostEqual(eos.T.magnitude, 200.0, places=3)
        self.assertAlmostEqual(eos.x, 0.5, places=3)
        self.assertTrue(hasattr(eos, 'Liquid'))
        self.assertTrue(hasattr(eos, 'Vapor'))
        self.assertAlmostEqual(eos.v, eos.x * eos.Vapor.v + (1 - eos.x) * eos.Liquid.v, places=3)
    
    def test_pengrobinson_resolve_Px(self):
        test_eos = PengRobinsonEOS(logiter=True, P=2.5).set_compound('ethane')
        Tsat = test_eos.Tsat
        eos = PengRobinsonEOS(P=2.5, x=0.35).set_compound('ethane')
        self.assertAlmostEqual(eos.T.magnitude, Tsat.magnitude, places=3)
        self.assertAlmostEqual(eos.x, 0.35, places=3)
        self.assertTrue(hasattr(eos, 'Liquid'))
        self.assertTrue(hasattr(eos, 'Vapor'))
        self.assertAlmostEqual(eos.v, eos.x * eos.Vapor.v + (1 - eos.x) * eos.Liquid.v, places=3)
    
    def test_pengrobinson_resolve_vx(self):
        test_eos = PengRobinsonEOS(logiter=True, T=200.0, x=0.5).set_compound('ethane')
        Pvap = test_eos.Pvap
        h_start = test_eos.h
        v_start = test_eos.v
        u_start = test_eos.u
        s_start = test_eos.s
        eos = PengRobinsonEOS(v=0.5 * test_eos.Vapor.v + 0.5 * test_eos.Liquid.v, x=0.5).set_compound('ethane')
        self.assertAlmostEqual(eos.P.magnitude, Pvap.magnitude, places=3)
        self.assertAlmostEqual(h_start.magnitude, eos.h.magnitude, places=3)
        self.assertAlmostEqual(v_start.magnitude, eos.v.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, eos.u.magnitude, places=3)
        self.assertAlmostEqual(s_start.magnitude, eos.s.magnitude, places=3)

    def test_pengrobinson_resolve_hx(self):
        test_eos = PengRobinsonEOS(logiter=True, T=200.0, x=0.5).set_compound('ethane')
        Pvap = test_eos.Pvap
        h_start = test_eos.h
        v_start = test_eos.v
        u_start = test_eos.u
        s_start = test_eos.s
        eos = PengRobinsonEOS(h=h_start, x=0.5).set_compound('ethane')
        self.assertAlmostEqual(eos.P.magnitude, Pvap.magnitude, places=3)
        self.assertAlmostEqual(h_start.magnitude, eos.h.magnitude, places=3)
        self.assertAlmostEqual(v_start.magnitude, eos.v.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, eos.u.magnitude, places=3)
        self.assertAlmostEqual(s_start.magnitude, eos.s.magnitude, places=3)

    def test_pengrobinson_resolve_sx(self):
        test_eos = PengRobinsonEOS(logiter=True, T=200.0, x=0.5).set_compound('ethane')
        Pvap = test_eos.Pvap
        h_start = test_eos.h
        v_start = test_eos.v
        u_start = test_eos.u
        s_start = test_eos.s
        eos = PengRobinsonEOS(s=s_start, x=0.5).set_compound('ethane')
        self.assertAlmostEqual(eos.P.magnitude, Pvap.magnitude, places=3)
        self.assertAlmostEqual(h_start.magnitude, eos.h.magnitude, places=3)
        self.assertAlmostEqual(v_start.magnitude, eos.v.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, eos.u.magnitude, places=3)
        self.assertAlmostEqual(s_start.magnitude, eos.s.magnitude, places=3)

    def test_pengrobinson_resolve_ux(self):
        test_eos = PengRobinsonEOS(logiter=True, T=200.0, x=0.5).set_compound('ethane')
        Pvap = test_eos.Pvap
        h_start = test_eos.h
        v_start = test_eos.v
        u_start = test_eos.u
        s_start = test_eos.s
        eos = PengRobinsonEOS(u=u_start, x=0.5).set_compound('ethane')
        self.assertAlmostEqual(eos.P.magnitude, Pvap.magnitude, places=3)
        self.assertAlmostEqual(h_start.magnitude, eos.h.magnitude, places=3)
        self.assertAlmostEqual(v_start.magnitude, eos.v.magnitude, places=3)
        self.assertAlmostEqual(u_start.magnitude, eos.u.magnitude, places=3)
        self.assertAlmostEqual(s_start.magnitude, eos.s.magnitude, places=3)
