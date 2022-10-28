
import numpy as np


class viper():
    def __init__(self):

        # --- miscellaneaous parameters ---
        self.engine_ID = None
        self.cruise = False
        self.A_t = None
        self.A_n = None
        self.pres_ratio = 5.5          # [-]
        self.gross_thrust = 15167   # [N]
        self.m_dot_air = 23.81      # [kg/s]
        self.m_dot_f = 0.4267       # [kg/s]
        self.T_0 = 288              # [K]
        self.p_0 = 100000           # [Pa]
        self.p_ref = 100000         # [Pa]
        self.T_ref = 288            # [K]
        self.tur_eff = 0.8          # [-]
        self.comp_eff = 0.886530850        # [-]0.886530850
        self.comb_eff = 1           # [-]
        self.nozz_eff = 1           # [-]
        self.cp_a = 1000            # [J/kg*K]
        self.cp_g = 1150            # [J/kg*K]
        self.k_a = 1.4              # [-]
        self.k_g = 1.33             # [-]
        self.LHV = 43*10**6         # [J]
        self.mech_eff = 0.99        # [-]
        self.R = 287                # [J/kg*K]
        self.nozzle_area = None     # [m^2]
        self.gross_thrust = None    # [N]

        # --- station parameters ---
        self.Tt_1 = None
        self.pt_1 = None
        self.Tt_2 = None
        self.pt_2 = None
        self.Tt_3 = None
        self.pt_3 = None

        # --- after combustion chamber ---
        self.Tt_4 = None
        self.pt_4 = None
        self.m_dot_4 = None
        self.Tt_5 = None
        self.pt_5 = None
        self.Tt_8 = None
        self.pt_8 = None
        self.rho_8 = None
        self.choked = False
        self.v_8 = None

    def pressure_turbine(self):
        '''
        calculate the exit presssure across a turbine stage as a function of the temperature ratio and the initial temperature
        '''

        tempr_ratio = self.Tt_5 / self.Tt_4
        self.pt_5 = self.pt_4 * (1 - (1 / self.tur_eff) * (1 - tempr_ratio)) ** (self.k_g / (self.k_g - 1))

    def compressor(self):
        '''
        calculate the exit temperature across a compressor stage as a function of the initial temperature and the pressure ratio
        '''

        self.pt_3 = self.pt_2 * self.pres_ratio
        self.Tt_3 = self.Tt_2 * (1 + (1 / self.comp_eff) * ((self.pres_ratio) ** ((self.k_a - 1) / self.k_a) - 1))

    def combustion(self):
        '''
        calculate the exit temperature of the combustion chamber as a function of CC inlet temp
        assuming fuel massflow is known
        '''

        self.pt_4 = self.pt_3
        self.Tt_4 = (self.m_dot_f*self.LHV*self.comb_eff)/(self.m_dot_air*self.cp_g) + self.Tt_3
        self.m_dot_4 = self.m_dot_f + self.m_dot_air

    def temperature_turbine(self):
        '''
        calculate exit temperature of turbine stage based on power balance with compressor
        '''

        self.Tt_5 = self.Tt_4 - (self.m_dot_air*self.cp_a*(self.Tt_3-self.Tt_2))/(self.m_dot_4*self.cp_g*self.mech_eff)

    def m_fuel_rate(self):
        '''
        calculate fuel massflow based on temperature difference across combustion chamber
        '''

        self.m_dot_f = (self.m_dot_air*self.cp_g*(self.Tt_4 - self.Tt_3))/(self.comb_eff*self.LHV)
        self.m_dot_4 = self.m_dot_air + self.m_dot_f
        self.pt_4 = self.pt_3

    def nozzle(self):
        '''
        caclulate flow condition at nozzle throat for either a choked or unchoked nozzle
        '''

        p_critical_ratio = (1-(1/self.nozz_eff)*((self.k_g-1)/(self.k_g+1)))**(-self.k_g/(self.k_g-1))

        if self.pt_5/self.p_0 > p_critical_ratio:
            self.choked = True

        if self.choked:
            self.pt_8 = self.pt_5 / p_critical_ratio                         # Calculate nozzle pressure using expansion up to M=1
            self.Tt_8 = self.Tt_5 * (2 / (self.k_g + 1))                     # calculate nozzle temperature using expansion up to M=1
            self.v_8 = np.sqrt(self.k_g * self.R * self.Tt_8)                # calculate speed of sound in nozzle, equal to jet velocity
            self.rho_8 = self.pt_8 / (self.R * self.Tt_8)                      # calculate density in the nozzle
            if not self.cruise:
                self.nozzle_area = self.m_dot_4 / (self.v_8 * self.rho_8)        # calculate nozzle area using continuity
            else:
                self.m_dot_4 = self.nozzle_area * self.v_8 * self.rho_8
            self.gross_thrust = self.m_dot_4 * self.v_8 + self.nozzle_area * (self.pt_8 - self.p_0)  # calculate gross thrust

        if not self.choked:
            self.pt_8 = self.p_0                                             # equate nozzle pressure to ambient pressure
            self.Tt_8 = self.Tt_5 - self.Tt_5 * self.nozz_eff * (1 - (self.pt_8 / self.pt_5) ** ((self.k_g-1) / self.k_g))    # calculate nozzle temperature based on expansion to P_ambient
            self.v_8 = np.sqrt(2 * self.cp_g * (self.Tt_5 - self.Tt_8))      # calculate jet velocity based on conservation of total enthalpy
            self.rho_8 = self.pt_8 / (self.R * self.Tt_8)                    # calculate density in the nozzle using ideal gas law
            if not self.cruise:
                self.nozzle_area = self.m_dot_4 / (self.v_8 * self.rho_8)        # calculate nozzle area using continuity
            else:
                self.m_dot_4 = self.nozzle_area * self.v_8 * self.rho_8
            self.gross_thrust = self.m_dot_4 * self.v_8                      # calculate gross thrust

        self.A_t = self.nozzle_area * (self.pt_5 / self.pt_4) ** ((2 * self.k_g - (self.k_g - 1)) / (2 * self.k_g))


    def printing(self, names, values, units):
        print('\n--- results for engine', self.engine_ID, '---')
        for value, name, unit in zip(values, names, units):
            print(name,' = ',  '%.6f' % value, ' [', unit, ']')

        if self.choked:
            print('this engine is choked')
        else:
            print('this engine is not choked')


#--------------------- exercise 1 ---------------------------------
engine1 = viper()
engine1.engine_ID = 1
engine1.Tt_2 = engine1.T_0
engine1.pt_2 = engine1.p_0

# calculate conditions after compressor
engine1.compressor()

# calculate conditions after combustion chamber
engine1.combustion()

# calculate conditions after the turbine
engine1.temperature_turbine()
engine1.pressure_turbine()

# calculate conditions in nozzle
engine1.nozzle()

# -------------------- exercise 2 --------------------
engine2 = viper()
engine2.engine_ID = 2
engine2.Tt_2 = engine2.T_0
engine2.pt_2 = engine2.p_0

# calculate properties after compressor stage
engine2.compressor()

# calculate properties after combustion stage
engine2.Tt_4 = 900 # [K]
engine2.m_fuel_rate()

engine2.temperature_turbine()
engine2.pressure_turbine()

# calculate conditions in nozzle
engine2.nozzle()

# -------------------- exercise 3 --------------------
engine3 = viper()
engine3.engine_ID = 3
engine3.cruise = True

engine3.p_0 = 0.2454 *10**5 # [Pa]
engine3.T_0 = 220 # [K]
engine3.pt_2 = engine3.p_0 * (1 + ((engine3.k_a - 1)/2)*0.78**2)**(engine3.k_a/(engine3.k_a-1))
engine3.Tt_2 = engine3.T_0 * (1+ ((engine3.k_a - 1)/ 2)*0.78**2)

engine3.Tt_4 = 1150  # [K]
engine3.m_dot_4 = engine3.m_dot_air
engine3.nozzle_area = engine2.nozzle_area

# calculate properties after compressor stage
for iteration in range(10):
    print(engine3.m_dot_air)
    engine3.compressor()

    # calculate properties after combustion stage
    if iteration != 0:
        engine3.m_fuel_rate()
    else:
        engine3.pt_4 = engine3.pt_3

    engine3.temperature_turbine()
    engine3.pressure_turbine()

    # calculate conditions in nozzle
    engine3.nozzle()
    if iteration == 0:
        engine3.m_dot_air = engine3.m_dot_4
        engine3.cruise = False

#--------------- printing results ------------------
names1 = ['turbine to nozzle area ratio', 'gross thrust', 'TIT']
values1 = [engine1.A_t/engine1.nozzle_area, engine1.gross_thrust, engine1.Tt_4]
units1 = ['-', 'N', 'K']
engine1.printing(names1, values1, units1)

engine2.printing(names1, values1, units1)

names3 = ['exit massflow', 'gross thrust', 'fuel massflow', 'inlet total pressure', 'inlet total temperature', 'compressor work', 'temp ratio']
values3 = [engine3.m_dot_4, engine3.gross_thrust, engine3.m_dot_f, engine3.pt_2, engine3.Tt_2, (engine3.Tt_3-engine3.Tt_2)*engine3.cp_a*engine3.m_dot_air, engine3.Tt_3/engine3.Tt_2]
units3 = ['kg/s', 'N', 'kg/s', 'Pa', 'K', 'w', '-']


engine3.printing(names3, values3, units3)
# ------------ Questions --------------------

# - How can the are ratio of turbine to nozzle be smaller than 1 if the turbine is already choked, wouldn't that mean that flow would go supersonic??


