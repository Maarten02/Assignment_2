
import numpy as np
import pandas as pd

class viper():
    def __init__(self):
        self.pres_ratio = 5.5          # [-]
        self.gross_thrust = 15167   # [N]
        self.m_dot_air = 23.81      # [kg/s]
        self.m_dot_f = 0.4267       # [kg/s]
        self.T_0 = 288              # [K]
        self.p_0 = 100000           # [Pa]
        self.tur_eff = 0.8          # [-]
        self.comp_eff = 0.645       # [-]
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
        return (self.m_dot_f*self.LHV*self.comb_eff)/(self.m_dot_air*self.cp_g) + initial_temperature

    def temperature_turbine(self, T_2, T_3, T_4, m_dot_5):
        return T_4 - (self.m_dot_air*self.cp_a*(T_3-T_2))/(m_dot_5*self.cp_g*self.mech_eff)

    def m_fuel_rate(self, T_3, T_4):
        return (self.m_dot_air*self.cp_g*(T_4-T_3))/(self.comb_eff*self.LHV)

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
            T_gross = m5 * v8 + nozzle_area*(nozzle_pressure-self.p_0)

        if not choked:
            nozzle_pressure = self.p_0
            nozzle_temperature = t5 - t5*self.nozz_eff*(1-(nozzle_pressure/p5)**((self.k_g-1)/self.k_g))
            v8 = np.sqrt(2*self.cp_g*(t5 - nozzle_temperature))



        return nozzle_pressure, nozzle_temperature, nozzle_area, T_gross, choked

#p0 = pt0 = pt1 = pt2
#T0 = Tt0 = Tt1= Tt2


engine1 = viper()

ambient_pressure = engine1.p_0

# calculate conditions after compressor

Pt3_1 = engine1.pres_ratio * engine1.p_0
Tt3_1 = engine1.temperature_compressor(engine1.T_0)

# print(Pt3, '\t', Tt3)

# calculate conditions after combustion chamber

Pt4_1 = Pt3_1
Tt4_1 = engine1.combustion(Tt3_1)
m_dot4_1 = engine1.m_dot_f + engine1.m_dot_air

# print(Pt4, '\t', Tt4)

# calculate conditions after the turbine
m_dot5_1 = m_dot4_1
Tt5_1 = engine1.temperature_turbine(engine1.T_0, Tt3_1, Tt4_1, m_dot5_1)
Pt5_1 = engine1.pressure_turbine(Pt4_1, Tt5_1/Tt4_1)


# print(Pt5, '\t', Tt5)

# calculate conditions in nozzle
nozzle_calculation = engine1.nozzle(Pt5_1, Tt5_1, m_dot5_1)
A_n = nozzle_calculation[2]
A_t = A_n*(Pt5_1/Pt4_1)**((2*engine1.k_g-(engine1.k_g-1))/(2*engine1.k_g))
choked = nozzle_calculation[4]
if choked:
    print('the nozzle is choked')
if not choked:
    print('the nozzle is not unchoked')

# -------------------- exercise 2 --------------------
engine2 = viper()


# calculate properties after compressor stage
Pt3_2 = engine2.pres_ratio * engine2.p_0
Tt3_2 = engine2.temperature_compressor(engine2.T_0)

# calculate properties after combustion stage
Tt4_2 = 900 # [K]
engine2.m_dot_f = engine2.m_fuel_rate(Tt3_2, Tt4_2)

Tt5_2 = engine2.temperature_turbine(engine2.p_0, Tt3_2, Tt4_2, engine2.m_dot_f+engine2.m_dot_air)
Pt_5 = engine2.pressure_turbine(Pt3_2, (Tt5_2/Tt4_2))


print(engine2.m_dot_f, '\t', Tt5_2, '\t', Pt_5)

