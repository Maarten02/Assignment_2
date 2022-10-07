
import numpy as np
import pandas as pd

class viper():
    def __init__(self):
        self.Pi_comp = 5.5          # [-]
        self.gross_thrust = 15167   # [N]
        self.m_dot_air = 23.81      # [kg/s]
        self.m_dot_f = 0.4267       # [kg/s]
        self.T_0 = 288              # [K]
        self.p_0 = 1                # [bar]
        self.tur_eff = 0.9          # [-]
        self.comp_eff = 0.92        # [-]
        self.comb_eff = 1           # [-]
        self.nozz_eff = 1           # [-]
        self.cp_a = 1000            # [J/kg*K]
        self.cp_g = 1150            # [J/kg*K]
        self.k_a = 1.4              # [-]
        self.k_g = 1.33             # [-]



    def scenario_1(self):
        '''
        Scenario 1:
        test-bed assuming choked turbine and propulsive nozzle
        calculate:

        is the nozzle choked?
        what is the turbine to propulsive nozzle area ratio
        '''



        return None


    def pressure_turbine(self, efficiency, initial_pressure, tempr_ratio, gas=True):

        exit_pressure = initial_pressure * (1 + efficiency)
        return exit_pressure


    def temperature_compressor(self, efficiency, initial_temperature, pres_ratio, gas=False):
        if gas == False:
            T_2 = initial_temperature*(1 + (1/efficiency)*((pres_ratio)**((self.k_a-1)/self.k_a)-1))

        else:
            T_2 = initial_temperature * (1 + (1 / efficiency) * ((pres_ratio) ** ((self.k_g - 1) / self.k_g) - 1))

        return T_2










