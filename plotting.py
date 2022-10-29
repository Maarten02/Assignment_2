import matplotlib.pyplot as plt
import numpy as np
import pyromat as pm


# Make meridional gas path plot - Maarten
def plot_mgas_path(stages, radii, x_locs):

    plt.figure()
    #plt.ylim([-1,1])
    plt.xlabel('Axial length [m]')
    plt.ylabel('Radius [m]')
    plt.grid()
    plt.title('Meridional Gas Path')

    fmts = ['r.-', 'b.-']

    for i in range(stages):

        xids = [[4 * i, 4 * i + 1], [4 * i + 2, 4 * i + 3]]
        yids = [[3 * i, 3 * i + 1], [3 * i + 1, 3 * i + 2]]

        for fmt, xid, yid in zip(fmts, xids, yids):

            # [[inner, outer], [inner, outer]]
            x_le = x_locs[xid[0]]
            x_te = x_locs[xid[1]]
            y_le = radii[yid[0]]
            y_te = radii[yid[1]]

            # plot rotor (red) or stator (blue)
            plt.vlines(x_le, y_le[1], y_le[0], colors=fmt[:-2])
            plt.vlines(x_te, y_te[1], y_te[0], colors=fmt[:-2])
            plt.plot([x_le, x_te],[y_le[0], y_te[0]], fmt)
            plt.plot([x_le, x_te],[y_le[1], y_te[1]], fmt)

    plt.savefig('figures/m_gas_path.pdf')


def plot_pressure_ratios(stage_pressure_ratios, n_stages_list):
    plt.figure(1)
    plt.plot(n_stages_list, stage_pressure_ratios, marker = '.')
    plt.title("Flow properties along Gas Path")
    plt.xlabel('Stages')
    plt.ylabel(r'Pressure Ratio, $\beta_{stage}$')
    plt.grid()
    plt.savefig('figures/pressure_ratio.pdf')

def plot_static_pres_over_total(inter_static_pres, first_total_pres, n_stages):
    inter_static_pres_array = np.array(inter_static_pres)
    static_total_ratio = inter_static_pres_array/first_total_pres
    plt.figure(2)
    x = list(range(n_stages * 3))
    x_ticks = np.zeros(len(x))
    m = 0
    for k in range(1,n_stages+1):
        x_ticks[m] = str(k) + '.1'
        x_ticks[m+1] = str(k) + '.2'
        x_ticks[m+2] = str(k) + '.3'
        m += 3
    plt.xticks(x, x_ticks)
    i = 0
    p = 0
    while i < n_stages*2:
        plt.plot(x[p:p+3], static_total_ratio[i:i+3], marker='.')
        i += 2
        p += 3
    plt.title(r"Flow properties along Gas Path, $\frac{P}{P_{t0}}$ ")
    plt.xlabel(r'Stages')
    plt.ylabel(r'$\frac{P}{P_{t0}}$')
    plt.grid()
    # plt.show()
    plt.savefig('figures/static_pres_over_total.pdf')

def plot_static_temp_over_total(inter_static_temp, first_total_temp, n_stages):
    inter_static_temp_array = np.array(inter_static_temp)
    static_total_ratio = inter_static_temp_array/first_total_temp
    plt.figure(3)
    x = list(range(n_stages*3))
    x_ticks = np.zeros(len(x))
    m = 0
    for k in range(1, n_stages + 1):
        x_ticks[m] = str(k) + '.1'
        x_ticks[m + 1] = str(k) + '.2'
        x_ticks[m + 2] = str(k) + '.3'
        m += 3
    plt.xticks(x, x_ticks)
    i = 0
    p = 0
    while i < n_stages*2:
        plt.plot(x[p:p+3], static_total_ratio[i:i+3], marker='.')
        i += 2
        p += 3
    plt.title(r"Flow properties along Gas Path, $\frac{T}{T_{t0}}$ ")
    plt.xlabel(r'Stages')
    plt.ylabel(r'$\frac{T}{T_{t0}}$')
    plt.grid()
    # plt.show()
    plt.savefig('figures/static_temp_over_total.pdf')

def plot_total_pres_over_total(inter_total_pres, first_total_pres, n_stages):
    inter_total_pres_array = np.array(inter_total_pres)
    total_total_ratio = inter_total_pres_array/first_total_pres
    plt.figure(4)
    x = list(range(n_stages*3))
    x_ticks = np.zeros(len(x))
    m = 0
    for k in range(1, n_stages + 1):
        x_ticks[m] = str(k) + '.1'
        x_ticks[m + 1] = str(k) + '.2'
        x_ticks[m + 2] = str(k) + '.3'
        m += 3
    plt.xticks(x, x_ticks)
    y_plot = np.zeros(3)
    i = 0
    p = 0
    while i < n_stages:
        y_plot[0] = total_total_ratio[i]
        y_plot[1] = total_total_ratio[i+1]
        y_plot[2] = total_total_ratio[i+1]
        plt.plot(x[p:p+3], y_plot, marker='.')
        i += 1
        p += 3
    plt.title(r"Flow properties along Gas Path, $\frac{P_t}{P_{t0}}$ ")
    plt.xlabel(r'Stages')
    plt.ylabel(r'$\frac{P_t}{P_{t0}}$')
    plt.grid()
    # plt.show()
    plt.savefig('figures/total_pres_over_total.pdf')

def plot_total_temp_over_total(inter_total_temp, first_total_temp, n_stages):
    inter_total_temp_array = np.array(inter_total_temp)
    total_total_ratio = inter_total_temp_array/first_total_temp
    plt.figure(5)
    x = list(range(n_stages*3))
    x_ticks = np.zeros(len(x))
    m = 0
    for k in range(1, n_stages + 1):
        x_ticks[m] = str(k) + '.1'
        x_ticks[m + 1] = str(k) + '.2'
        x_ticks[m + 2] = str(k) + '.3'
        m += 3
    plt.xticks(x, x_ticks)
    y_plot = np.zeros(3)
    i = 0
    p = 0
    while i < n_stages:
        y_plot[0] = total_total_ratio[i]
        y_plot[1] = total_total_ratio[i+1]
        y_plot[2] = total_total_ratio[i+1]
        plt.plot(x[p:p+3], y_plot, marker='.')
        i += 1
        p += 3
    plt.title(r"Flow properties along Gas Path, $\frac{T_t}{T_{t0}}$ ")
    plt.xlabel(r'Stages')
    plt.ylabel(r'$\frac{T_t}{T_{t0}}$')
    plt.grid()
    # plt.show()
    plt.savefig('figures/total_temp_over_total.pdf')



def plot_h_s_diagram(stage_temperatures_total, stage_pressures_total, stage_temperatures_static, stage_pressures_static, n_stages):

    plt.figure()
    air = pm.get('ig.air')

    # enthalpy_stage = air.h(T=stage_temperatures)
    # entropy_stage = air.s(T=stage_temperatures, p=stage_pressures)
    # plt.figure(6)
    # plt.plot(entropy_stage,enthalpy_stage, marker = '.')
    # for i in range(n_stages+1):
    #     T_plot = stage_temperatures[i]
    #     plt.plot(air.s(T=[T_plot-10, T_plot, T_plot+10], p=stage_pressures[i]), air.h(T =[T_plot-10, T_plot, T_plot+10]))

    fmts = ['r.--', 'g.--', 'b.--', 'c.--', 'm.--']
    for i, fmt in zip(range(n_stages), fmts):

        # plot lines between (h,s) evaluated at ((pt_i1, Tt_i1), (pt_i2, Tt_i2), (pt_i3, Tt_i3))
        # stage 0 --> Pt/Tt[0, 1, 1] || Ps [0, 1, 2]
        # stage 1 --> Pt/Tt[1, 2, 2] || Ps [2, 3, 4]
        # stage 2 -->                || Ps [4, 5, 6]
        t_ids = [i, i+1, i+1]
        s_ids = [2*i, 2*i+1, 2*i+2]

        enthalpies_total = []
        entropies_total = []
        enthalpies_static = []
        entropies_static = []
        for t_id, s_id in zip(t_ids, s_ids):
        # for t_id in t_ids:


            enthalpies_static.append(air.h(T=stage_temperatures_static[s_id]))
            entropies_static.append(air.s(T=stage_temperatures_static[s_id], p=stage_pressures_static[s_id]))
            enthalpies_total.append(air.h(T=stage_temperatures_total[t_id]))
            entropies_total.append(air.s(T=stage_temperatures_total[t_id], p=stage_pressures_total[t_id]))


        plt.plot(entropies_static, enthalpies_static, fmt)
        plt.plot(entropies_total, enthalpies_total, fmt[:-1])

        # plot lines between (h,s) evaluated at ((p_i1, T_i1),(p_i2, T_i2),(p_i3, T_i3))



        # enthalpy_stage = air.h(T=stage_temperatures)
        # entropy_stage = air.s(T=stage_temperatures, p=stage_pressures)

    plt.title("h-s Diagram")
    plt.xlabel(r'Entropy, $s$, [$J/kg \cdot K$]')
    plt.ylabel(r'Enthalpy, $h$, [$J/K$]')
    plt.grid()
    plt.savefig('figures/h_s_diagram.pdf')
    plt.show()
