
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
        self.LHV = 43*10**6         # [J]
        self.mech_eff = 0.99        # [-]



    def scenario_1(self):
        '''
        Scenario 1:
        test-bed assuming choked turbine and propulsive nozzle
        calculate:

        is the nozzle choked?
        what is the turbine to propulsive nozzle area ratio
        '''



        return None


    def pressure_turbine(self, initial_pressure, tempr_ratio, gas=True):
        if gas == True:
            exit_pressure = initial_pressure * (1-(1/self.tur_eff)*(1-tempr_ratio))**(self.k_g/(self.k_g-1))
        else:
            exit_pressure = initial_pressure * (1 - (1 / self.tur_eff) * (1 - tempr_ratio)) ** (self.k_a / (self.k_a - 1))
        return exit_pressure


    def temperature_compressor(self, initial_temperature, pres_ratio, gas=False):
        if gas == False:
            exit_temperature = initial_temperature*(1 + (1/self.comp_eff)*((pres_ratio)**((self.k_a-1)/self.k_a)-1))

        else:
            exit_temperature = initial_temperature * (1 + (1 / self.comp_eff) * ((pres_ratio) ** ((self.k_g - 1) / self.k_g) - 1))

        return exit_temperature

    def combustion(self, initial_temperature):
        return (self.m_dot_f*self.LHV*self.comb_eff)/(self.m_dot_air*self.cp_a) + initial_temperature

    def temperature_turbine(self, T_2, T_3, T_4):
        return T_4 - (self.m_dot_air*self.cp_a*(T_3-T_2))/(self.m_dot_f*self.cp_g*self.mech_eff)

#p0 = pt0 = pt1 = pt2
#T0 = Tt0 = Tt1= Tt2












