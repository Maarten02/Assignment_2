
import numpy as np
import pandas as pd

class viper():
    def __init__(self):
        self.pres_ratio = 5.5          # [-]
        self.gross_thrust = 15167   # [N]
        self.m_dot_air = 23.81      # [kg/s]
        self.m_dot_f = 0.4267       # [kg/s]
        self.T_0 = 288              # [K]
        self.p_0 = 100000            # [Pa]
        self.tur_eff = 0.8          # [-]
        self.comp_eff = 0.645        # [-]
        self.comb_eff = 1           # [-]
        self.nozz_eff = 1           # [-]
        self.cp_a = 1000            # [J/kg*K]
        self.cp_g = 1150            # [J/kg*K]
        self.k_a = 1.4              # [-]
        self.k_g = 1.33             # [-]
        self.LHV = 43*10**6         # [J]
        self.mech_eff = 0.99        # [-]
        self.R = 287                # [J/kg*K]



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


    def temperature_compressor(self, initial_temperature, gas=False):
        if gas == False:
            exit_temperature = initial_temperature*(1 + (1/self.comp_eff)*((self.pres_ratio)**((self.k_a-1)/self.k_a)-1))

        else:
            exit_temperature = initial_temperature * (1 + (1 / self.comp_eff) * ((self.pres_ratio) ** ((self.k_g - 1) / self.k_g) - 1))

        return exit_temperature

    def combustion(self, initial_temperature):
        return (self.m_dot_f*self.LHV*self.comb_eff)/(self.m_dot_air*self.cp_a) + initial_temperature

    def temperature_turbine(self, T_2, T_3, T_4, m_dot_5):
        return T_4 - (self.m_dot_air*self.cp_a*(T_3-T_2))/(m_dot_5*self.cp_g*self.mech_eff)

    def nozzle(self, p5, t5, m5, v_inf = 0):
        choked = False
        p_critical_ratio = (1- (1/self.nozz_eff)*((self.k_g-1)/(self.k_g+1)))**(-self.k_g/(self.k_g-1))
        if p5/self.p_0 > p_critical_ratio:
            choked = True
        if choked:
            nozzle_pressure = p5 / p_critical_ratio
            nozzle_temperature = t5 * (2 / (self.k_g + 1))
            v8 = np.sqrt(self.k_g*self.R*nozzle_temperature)
            rho8 = nozzle_pressure/(self.R*nozzle_temperature)
            nozzle_area = m5/(v8*rho8)
            T_gross = m5* v8 + nozzle_area*(nozzle_pressure-self.p_0)

        if not choked:
            nozzle_pressure = self.p_0
            nozzle_temperature = None


        return nozzle_pressure, nozzle_temperature, nozzle_area, T_gross

#p0 = pt0 = pt1 = pt2
#T0 = Tt0 = Tt1= Tt2


engine1 = viper()

ambient_pressure = engine1.p_0

# calculate conditions after compressor

Pt3 = engine1.pres_ratio * engine1.p_0
Tt3 = engine1.temperature_compressor(engine1.T_0)

# print(Pt3, '\t', Tt3)

# calculate conditions after combustion chamber

Pt4 = Pt3
Tt4 = engine1.combustion(Tt3)
m_dot4 = engine1.m_dot_f + engine1.m_dot_air

# print(Pt4, '\t', Tt4)

# calculate conditions after the turbine
m_dot5 = m_dot4
Tt5 = engine1.temperature_turbine(engine1.T_0, Tt3, Tt4, m_dot5)
Pt5 = engine1.pressure_turbine(Pt4, Tt5/Tt4)


# print(Pt5, '\t', Tt5)

# calculate conditions in nozzle
nozzle_calculation = engine1.nozzle(Pt5, Tt5, m_dot5)
print(nozzle_calculation)



