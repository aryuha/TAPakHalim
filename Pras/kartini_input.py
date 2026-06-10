import math

from pyrk.utilities.ur import units
from pyrk import th_component as th
from pyrk.timer import Timer
from pyrk.materials.material import Material
from pyrk.materials.liquid_material import LiquidMaterial
from pyrk.density_model import DensityModel
from pyrk.reactivity_insertion import ReactivityInsertion
from pyrk.materials.graphite import Graphite

# -----------------------------------------------------------------------------
# USER PARAMETERS
# -----------------------------------------------------------------------------

TOTAL_WORTH_DOLLAR = 1.389  
BETA_EFF           = 0.0069
T_RAMP             = 60    # Waktu penarikan batang kendali

N_FUEL_ELEMENTS    = 71      # jumlah bahan bakar 
N_DUMMY_ELEMENTS   = 15      # jumlah dummy fuel

rho_total_dk = TOTAL_WORTH_DOLLAR * BETA_EFF  
T_START      = 0.0     # delay before rod starts moving [seconds]

# -----------------------------------------------------------------------------
# SIMULATION TIME
# -----------------------------------------------------------------------------

t0 = 0.0  * units.seconds
tf = 80.0 * units.seconds   
dt = 0.005 * units.seconds

ti = Timer(t0=t0, tf=tf, dt=dt)

n_pg        = 6       
n_dg        = 0       
kappa       = 0.0     
n_ref       = 0      
fission_iso = "u235"  
spectrum    = "thermal"
feedback    = False  
nsteps      = 16000

class SineSquaredReactivityInsertion(ReactivityInsertion):
    
    def __init__(self, timer, rho_total_dk, t_ramp_s, t_start_s=5.0):
        self._rho_total = rho_total_dk
        self._t_ramp    = t_ramp_s
        self._t_start   = t_start_s   # delay before rod moves [s]
        self._dt        = timer.dt.magnitude  # [s] - must be set before super()
        ReactivityInsertion.__init__(self, timer=timer)

    def f(self, t_idx):
        t = t_idx * self._dt

        if t <= self._t_start:
            return 0.0 * units.delta_k
        t_rod = t - self._t_start

        if t_rod >= self._t_ramp:
            return self._rho_total * units.delta_k

        x   = t_rod / self._t_ramp
        val = self._rho_total * (x - math.sin(2.0 * math.pi * x)
                                    / (2.0 * math.pi))
        return val * units.delta_k


rho_ext = SineSquaredReactivityInsertion(
    timer        = ti,
    rho_total_dk = rho_total_dk,
    t_ramp_s     = T_RAMP,
    t_start_s    = T_START
)

class LightWater(LiquidMaterial):

    def __init__(self, name="light_water"):
        LiquidMaterial.__init__(
            self,
            name = name,
            k    = self.thermal_conductivity(),
            cp   = self.specific_heat_capacity(),
            dm   = self.density()
        )

    def thermal_conductivity(self):
        return 0.609 * units.watt / (units.meter * units.kelvin)

    def specific_heat_capacity(self):
        return 4182.0 * units.joule / (units.kg * units.kelvin)

    def density(self):
        return DensityModel(a=997.0 * units.kg / (units.meter**3),
                            model="constant")


class UZrHFuel(Material):

    def __init__(self, name="uzrh_fuel"):
        Material.__init__(
            self,
            name = name,
            k    = self.thermal_conductivity(),
            cp   = self.specific_heat_capacity(),
            dm   = self.density()
        )

    def thermal_conductivity(self):
        return 18.0 * units.watt / (units.meter * units.kelvin)

    def specific_heat_capacity(self):
        return 473.0 * units.joule / (units.kg * units.kelvin)

    def density(self):
        return DensityModel(a=5900.0 * units.kg / (units.meter**3),
                            model="constant")

h_core = 0.58  * units.meter    
r_fuel = 0.018 * units.meter    


vol_one_fuel  = math.pi * r_fuel**2 * h_core       
vol_one_dummy = math.pi * r_fuel**2 * h_core       


vol_fuel  = vol_one_fuel  * N_FUEL_ELEMENTS         
vol_dummy = vol_one_dummy * N_DUMMY_ELEMENTS        
vol_cool  = 5.0e-5 * units.meter**2 * h_core       


a_fuel  = 2.0 * math.pi * r_fuel * h_core * N_FUEL_ELEMENTS
a_dummy = 2.0 * math.pi * r_fuel * h_core * N_DUMMY_ELEMENTS

v_cool = 0.10  * units.meter / units.second
h_conv = 3.0e4 * units.watt / (units.meter**2 * units.kelvin)

TOTAL_POWER    = 1.0 * units.watt

FUEL_FRAC      = 0.97    
DUMMY_FRAC     = 0.03    

power_fuel  = TOTAL_POWER * FUEL_FRAC   
power_dummy = TOTAL_POWER * DUMMY_FRAC  

t_fuel  = units.Quantity(30.0, units.degC).to(units.kelvin)
t_dummy = units.Quantity(30.0, units.degC).to(units.kelvin)
t_cool  = units.Quantity(30.0, units.degC).to(units.kelvin)
t_inlet = units.Quantity(27.0, units.degC).to(units.kelvin)


alpha_fuel  = 0.0 * units.delta_k / units.kelvin   
alpha_dummy  = 0.0 * units.delta_k / units.kelvin  
alpha_cool  = 0.0 * units.delta_k / units.kelvin   

fuel = th.THComponent(
    name      = "fuel",
    mat       = UZrHFuel(name="uzrh_fuel"),
    vol       = vol_fuel,
    T0        = t_fuel,
    alpha_temp= alpha_fuel,
    timer     = ti,
    heatgen   = True,
    power_tot = power_fuel
)

dummy = th.THComponent(
    name      = "dummy",
    mat       = Graphite(name="graphite_dummy"),
    vol       = vol_dummy,
    T0        = t_dummy,
    alpha_temp= alpha_dummy,
    timer     = ti,
    heatgen   = True,        
    power_tot = power_dummy
)

cool = th.THComponent(
    name      = "cool",
    mat       = LightWater(name="cool_water"),
    vol       = vol_cool,
    T0        = t_cool,
    alpha_temp= alpha_cool,
    timer     = ti
)

inlet = th.THComponent(
    name      = "inlet",
    mat       = LightWater(name="inlet_water"),
    vol       = vol_cool,
    T0        = t_inlet,
    alpha_temp= 0.0 * units.delta_k / units.kelvin,
    timer     = ti
)

# Fuel <-> Coolant
fuel.add_convection('cool',  h=h_conv, area=a_fuel)
cool.add_convection('fuel',  h=h_conv, area=a_fuel)

# Dummy <-> Coolant
dummy.add_convection('cool', h=h_conv, area=a_dummy)
cool.add_convection('dummy', h=h_conv, area=a_dummy)

# Coolant natural circulation
cool.add_mass_trans('inlet', H=h_core, u=v_cool)

components = [fuel, dummy, cool, inlet]