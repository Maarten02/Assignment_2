import numpy as np
import math
from matplotlib import pyplot as plt
import pyromat as pm

class compressor():
    def __init__(self):
        self.k = 1.4 # [-]
        self.R = 287 # [J/kg*K]
        self.inlet_pressure = 36678.87 # [pa]
        self.inlet_temperature = 246.77 # [K]
        self.inlet_static_pressure = 0.2454 * 10 ** 5 # [pa]
        self.inlet_static_temperature = 220 [K]
        self.inlet_density = self.inlet_pressure / (self.R * self.inlet_temperature)

        self.n_stages = 5
        self.phi_stage = 0.7
        self.psi_stage = 0.4
        self.deg_of_reaction_stage = 0.5

        self.stage_areas = []
        self.stage_pressures = []
        self.stage_pressure_ratios = []
        self.stage_pressures_static = []
        self.stage_temperatures = []
        self.stage_temperatures_static = []
        self.stage_densities = []
        self.velocities = []

        self.hub_radii = []
        self.tip_radii = []
        self.mean_radius = None
        self.stage_efficiency = None

        self.compressor_power = 1973832.940054       # [W], from main.py, it1: 2512166.60 it2:1928808.35
        self.Omega = 10000 * 2 * math.pi / 60   # [rad/s]
        self.work_per_stage = self.compressor_power / self.n_stages


        self.U_meanline = None


        self.m_dot_air = 11.299775384849212 # [kg/s] it1:12.6534 it2: 11.20978 (0.92 efficiency poly)
        self.spec_work_stage = self.work_per_stage / self.m_dot_air

        self.alpha_1 = math.radians(23.2)
        self.alpha_2 = math.radians(45)
        self.beta_1 = math.radians(-45)
        self.beta_2 = math.radians(-23.2)

        '''
        ----------------------Iterate--------------------------
        1. guess number of stages
        2. divide specific work by number of stages
        3. compute the U_mean by using specific work and work coeff (U_mean = sqrt(specific_work/work_coeff))
        4. compute R_mean using U_mean and Omega (R_mean = U_mean/Omega)
        5. (velocity traingles)
        6a. Compute Tt_i+1/Tt_i and after that the pressure ratio per stage, using the polytropic efficiency = 0.92 from the plot, formula from reader
        6b. Use the over all p_ratio and T_ratio to find isentropic efficiency
        6c. update the total compressor power
        7. set up continuity equation for the first stage
        8. determine the A using the continuity equation
        9. Determine R_tip (the blade height) per stage (R_tip differs per stage due to varying density)
        10. CHECK: U_tip < 450-600 m/s, using R and Omega
        -------------------------------------
        11. Determine the blade axial chords based on a suitable value of the aspect
            ratio (h/b) or by selecting an appropriate value of blade chord (e.g. in the range of 0.03-0.05 m)
            (see literature: Lewis R.I., “Turbomachinery performance analysis”)
        12. meridional gas path and hs diagram
        13. Overleaf!!!!!!!!!!!!!!!!!!!!!!!!!!!
        '''

    def U_mean(self):
        self.U_meanline = np.sqrt(self.spec_work_stage/self.psi_stage)

    def Radius_mean(self):
        self.mean_radius = self.U_meanline / self.Omega

    def v_abs_exit(self, U):
        return U / (np.sin(self.alpha_2)-np.cos(self.alpha_2)*np.tan(self.beta_2))

    def Mach(self, T, abs_velocity):
        return abs_velocity/np.sqrt(T*self.k*self.R)

    def isen_eff(self, T_ratio, p_ratio):
         return (p_ratio ** ((self.k - 1) / self.k) - 1) / (T_ratio - 1)

    def T_ratio(self, t_last):
        delta_t = self.spec_work_stage / 1000
        t_new = t_last + delta_t
        return t_new / t_last

    def p_ratio(self, T_ratio):
        return (T_ratio) ** ((self.k * self.stage_efficiency) / (self.k - 1))

    def density(self, pressure, temperature):
        return pressure / (self.R * temperature)

    def blade_length(self, density, abs_velocity):
        ''' assumption: velocity at mean radius is the mean velocity along the entire blade '''
        area = self.m_dot_air/(density*abs_velocity)
        return area/(2*np.pi*self.mean_radius)

    def tip_radius(self, blade_length):
        return self.mean_radius + blade_length*0.5

    def axial_outlet_velocity(self, v_exit):
        return v_exit*np.cos(self.alpha_2)

    def total_temp_to_static(self, T_t, v):
        return T_t - (v**2/(2*1000))

    def total_pressure_to_static(self, T_s, T_t, p_t):
        return p_t*((T_t/T_s)**((self.k-1)/self.k))

    def stage_loop(self):

        self.stage_temperatures.append(self.inlet_temperature)
        self.stage_pressures.append(self.inlet_pressure)
        self.stage_densities.append(self.inlet_density)
        self.stage_pressures_static.append(self.inlet_static_pressure)
        self.stage_temperatures_static.append(self.inlet_static_temperature)
        self.stage_efficiency = 0.91

        for i in range(self.n_stages):

            # here the meanline velocity, radius, abs and axial velocity are calculated
            self.U_mean()
            self.Radius_mean()
            exit_v = self.v_abs_exit(self.U_meanline)       # exit_v = v2
            axial_v = self.axial_outlet_velocity(exit_v)
            #axial v van 1 en 3 bepalen
            # total temp across a stage
            total_temp_1 = self.stage_temperatures[-1]
            total_temp_3 = total_temp_1 * self.T_ratio(total_temp_1)
            self.stage_temperatures.append(total_temp_3)

            # total pres across a stage
            total_pres_1 = self.stage_pressures[-1]
            total_pres_3 = total_pres_1 * self.p_ratio(total_temp_3/total_temp_1)
            self.stage_pressures.append(total_pres_3)

            # static temp across a stage
            static_temp_1 = self.stage_temperatures_static[-1]
            static_temp_2 = self.total_temp_to_static(total_temp_1, exit_v)
            static_temp_3 = static_temp_1 + (total_temp_3 - total_temp_1)
            self.stage_temperatures_static.extend([static_temp_2, static_temp_3])

            # static pres across a stage
            static_pres_1 = self.stage_pressures_static[-1]
            static_pres_2 = self.total_pressure_to_static(static_temp_2, total_temp_3, total_pres_3)
            static_pres_3 = static_temp_1 + (total_temp_3 - total_temp_1)
            self.stage_pressures_static.extend([static_pres_2, static_pres_3])

            # below is not correct
            density_1 = static_pres_1 / (static_temp_1 * self.R)
            density_2 = static_pres_2 / (static_temp_2 * self.R)
            density_3 = static_pres_3 / (static_temp_3 * self.R)

            # self.stage_densities.extend(density_2, density_3)

            blade_length_1 = self.blade_length(density_1, axial_v_1)
            blade_length_2 = self.blade_length(density_2, axial_v_2)
            blade_length_3 = self.blade_length(density_3, axial_v_3)
            tip_radius_1 = self.tip_radius(blade_length_1)
            tip_radius_2 = self.tip_radius(blade_length_2)
            tip_radius_3 = self.tip_radius(blade_length_3)
            self.tip_radii.extend([tip_radius_1, tip_radius_2, tip_radius_3])

            tip_speed_1 = self.Omega * tip_radius_1
            tip_speed_2 = self.Omega * tip_radius_2
            tip_speed_3 = self.Omega * tip_radius_3


            self.stage_pressure_ratios.append(next_pres/last_pres)
            print('tip speed at stage', i, '=', '%.2f' % tip_speed, '[m/s]')
            print('pressure ratio for stage', i, '=', '%.2f' % (next_pres/last_pres))
            print('blade length = %.4f' % blade_length, '[m]')
            print('tip radius = %.4f' % tip_radius, '[m]')
            print('mean radius = %.4f' % self.mean_radius, '[m]')
            if tip_speed > 450:
                print("the speed limit is reached in stage ", i)
            print('---------------------------------')
        print('opr =', (self.stage_pressures[-1]/self.stage_pressures[0]))
        print('otr =', (self.stage_temperatures[-1]/self.stage_temperatures[0]))

    def power_loop(self):
        # determine isentropic efficiency
        T_ratio = self.stage_temperatures[-1] / self.stage_temperatures[0]
        p_ratio = self.stage_pressures[-1] / self.stage_pressures[0]
        isen_eff = self.isen_eff(T_ratio, p_ratio)
        print('isentropic efficiency = %.9f' % isen_eff)
        # update compressor power using new found isentropic efficiency





    # --- QUESTIONS ---

    # - (can we assume) does the axial velocity component remain constant throughout the entire compressor?
    # yes - for constant mean radius and repeated stage it is a valid assumption
    # - should we match compressor with turbine and nozzle?
    #  no really
    # - what is stage total to total efficiency?
    # use
    # - pressure ratio from temp ratio




viper_compressor = compressor()
viper_compressor.stage_loop()
viper_compressor.power_loop()

print(viper_compressor.stage_pressure_ratios)



n_stages = list(range(1,viper_compressor.n_stages+1))
print(n_stages)
plt.figure(1)
plt.plot(n_stages,viper_compressor.stage_pressure_ratios, marker = '.')
plt.title("Flow properties along Gas Path")
plt.xlabel('Stages')
plt.ylabel(r'Pressure Ratio, $\beta_{stage}$')
plt.grid()

air = pm.get('ig.air')
enthalpy_stage = air.h(T = viper_compressor.stage_temperatures)
entropy_stage = air.s(T=viper_compressor.stage_temperatures, p= viper_compressor.stage_pressures)

plt.figure(2)
plt.plot(entropy_stage,enthalpy_stage, marker = '.')
for i in range(viper_compressor.n_stages):
    plt.plot(air.s(T=viper_compressor.stage_temperatures, p= viper_compressor.stage_pressures[i]), air.h(T = viper_compressor.stage_temperatures))
plt.title("h-s diagram")
plt.xlabel(r'Entropy, $s$')
plt.ylabel(r'Enthalpy, $h$')
plt.grid()
plt.show()




