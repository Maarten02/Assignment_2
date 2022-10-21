import numpy as np
import math


class compressor():
    def __init__(self):
        self.k = 1.4 # [-]
        self.R = 287 # [J/kg*K]
        self.stages = []
        self.overall_pressure_ratio = 5.5
        self.inlet_pressure = 36678.87
        self.inlet_temperature = 246.77
        self.inlet_density = self.inlet_pressure / (self.R * self.inlet_temperature)

        self.n_stages = 7
        self.phi_stage = 0.7
        self.psi_stage = 0.4
        self.deg_of_reaction_stage = 0.5

        self.stage_areas = []
        self.stage_pressures = []
        self.stage_temperatures = []
        self.stage_densities = []
        self.velocities = []

        self.hub_radii = []
        self.tip_radii = []
        self.mean_radius = None

        self.compressor_power = 2512166.60      # [W], from main.py
        self.Omega = 10000 * 2 * math.pi / 60   # [rad/s]
        self.work_per_stage = self.compressor_power / self.n_stages


        self.U_meanline = None


        self.m_dot_air = 12.6534 # [kg/s]
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
        6. Compute Tt_i+1/Tt_i and after that the pressure ratio per stage, using the isentropic total to total efficiency
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

    def stage_eff(self, v):
        return 1/(1+ (v**2/(2*self.spec_work_stage)))

    def T_ratio(self, M):
        return 1 + M**2*(self.k-1)*self.psi_stage

    def p_ratio(self, T_ratio, stage_eff):
        return (1 + stage_eff*(T_ratio-1))**(self.k / (self.k-1))

    def density(self, pressure, temperature):
        return pressure / (self.R* temperature)

    def blade_length(self, density, abs_velocity):
        '''assumption: velocity at mean radius is the mean velocity along the entire blade'''
        area = self.m_dot_air/(density*abs_velocity)
        return area/(2*np.pi*self.mean_radius)

    def tip_radius(self, blade_length):
        return self.mean_radius + blade_length*0.5

    def axial_outlet_velocity(self, v_exit):
        return v_exit*np.cos(self.alpha_2)

    def stage_loop(self):

        self.stage_temperatures.append(self.inlet_temperature)
        self.stage_pressures.append(self.inlet_pressure)
        self.stage_densities.append(self.inlet_density)


        for i in range(self.n_stages):

            # exit_v = v2
            self.U_mean()
            self.Radius_mean()
            exit_v = self.v_abs_exit(self.U_meanline)
            axial_v = self.axial_outlet_velocity(exit_v)
            last_temp = self.stage_temperatures[-1]
            mach = self.Mach(last_temp, exit_v)
            next_temp = last_temp * self.T_ratio(mach)
            self.stage_temperatures.append(next_temp)
            last_pres = self.stage_pressures[-1]
            # stage_efficiency = self.stage_eff(exit_v)
            stage_efficiency = 0.985
            next_pres = last_pres * self.p_ratio(next_temp/last_temp, stage_efficiency)

            self.stage_pressures.append(next_pres)

            next_density = next_pres / (next_temp * self.R)
            self.stage_densities.append(next_density)
            blade_length = self.blade_length(next_density,axial_v)
            tip_radius = self.tip_radius(blade_length)
            self.tip_radii.append(tip_radius)

            tip_speed = self.Omega*tip_radius
            print('tip speed at stage', i, '=', '%.2f' % tip_speed)
            print('pressure ratio for stage', i, '=', '%.2f' % (next_pres/last_pres))
            if tip_speed > 450:
                print("the speed limit is reached in stage ", i)

        print('opr =', (self.stage_pressures[-1]/self.stage_pressures[0]))
    # --- QUESTIONS ---

    # - (can we assume) does the axial velocity component remain constant throughout the entire compressor?
    # - should we match compressor with turbine and nozzle?
    # - what is stage total to total efficiency?
    # - pressure ratio from temp ratio

viper_compressor = compressor()
viper_compressor.stage_loop()



