from unittest import TestCase
from sandlercubics.eos import IdealGasEOS, PengRobinsonEOS, GeneralizedVDWEOS, SoaveRedlichKwongEOS
from sandlerprops.properties import get_database
class TestCubicEOS(TestCase):

    def setUp(self):
        # Load a sample component from the database for testing
        self.methane = get_database().get_compound('methane')
        self.water = get_database().get_compound('water')

    def test_ideal_gas_properties(self):
        eos = IdealGasEOS(pressure_unit="pa")
        eos.T = 300.0  # K
        eos.P = 101325.0  # Pa

        # Ideal gas compressibility factor should be 1
        self.assertAlmostEqual(eos.Z, 1.0, places=6)

        # Ideal gas molar volume should be RT/P
        R = 8.314  # J/(mol·K)
        V_ideal = R * eos.T / eos.P
        self.assertAlmostEqual(eos.v, V_ideal, places=4)

        # Fugacity coefficient should be 1 for ideal gas
        self.assertAlmostEqual(eos.phi, 1.0, places=6)

    def test_peng_robinson_gas_properties(self):
        Tc = self.methane.Tc
        Pc = self.methane.Pc
        omega = self.methane.Omega
        eos = PengRobinsonEOS(Tc=Tc, Pc=Pc, omega=omega, pressure_unit="bar")
        eos.T = 300.0  # K
        eos.P = 1.01325  # bar

        z = eos.Z
        self.assertAlmostEqual(z, 0.9978, places=3)

        eos.P = 50.0  # bar
        z = eos.Z
        self.assertAlmostEqual(z, 0.9021, places=3)

    def test_peng_robinson_water_saturation(self):
        Tc = self.water.Tc
        Pc = self.water.Pc
        omega = self.water.Omega
        eos = PengRobinsonEOS(Tc=Tc, Pc=Pc, omega=omega, pressure_unit="bar", logiter=True)
        eos.T = 373.15
        
        Pvap = eos.Pvap
        self.assertAlmostEqual(Pvap, 0.963, places=2)  # Expect ~1 atm in bar, but PR is not so accurate for water

        eos.P = 2.0
        Tsat = eos.Tsat
        self.assertAlmostEqual(Tsat, 394.413, places=2)

    def test_soave_redlich_kwong_gas_properties(self):
        Tc = self.methane.Tc
        Pc = self.methane.Pc
        omega = self.methane.Omega
        eos = SoaveRedlichKwongEOS(Tc=Tc, Pc=Pc, omega=omega, pressure_unit="bar")
        eos.T = 300.0  # K
        eos.P = 1.01325  # bar

        z = eos.Z
        self.assertAlmostEqual(z, 0.9985, places=3)

        eos.P = 50.0  # bar
        z = eos.Z
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
        Tc = self.water.Tc
        Pc = self.water.Pc
        omega = self.water.Omega
        eos = PengRobinsonEOS(Tc=Tc, Pc=Pc, omega=omega, pressure_unit="bar", logiter=True)
        eos.T = 373.15
        
        Hvap = eos.Hvap
        self.assertAlmostEqual(Hvap, 42073, delta=500)  # J/mol, PR is approximate

    def test_peng_robinson_entropy_of_vaporization(self):
        Tc = self.water.Tc
        Pc = self.water.Pc
        omega = self.water.Omega
        eos = PengRobinsonEOS(Tc=Tc, Pc=Pc, omega=omega, pressure_unit="bar", logiter=True)
        eos.T = 373.15
        
        Svap = eos.Svap
        self.assertAlmostEqual(Svap, 112.8, delta=5)  # J/mol-K, PR is approximate

        # check clausius-clapeyron relation
        Svap_CC = eos.Hvap / eos.T
        self.assertAlmostEqual(Svap, Svap_CC, delta=5)

