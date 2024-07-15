from Solvers.OperatingPointSolver import engOperation
from Miscellaneous.AuxilliaryFunctions import componentPlot
from Components.DesignVariables import N_ref_c, N_ref_t

import warnings, numpy as np, matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

# OPERATING POINT ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Arguments: Flight Mach Number and Compressor Relative Corrected Spool Speed

plot = True

# Choose entries: ------------------------------------------------------------------------------------------------------------------------------------------

M0 = 0.8
N_c = 0.9965
nozzle = 'conv'

# Loop parameters: -----------------------------------------------------------------------------------------------------------------------------------------

num_iter0 = 15
relaxation_factor = 0.15

m_0, T2t_T0, p2t_p0, eta_d, m_2, T3t_T2t, p3t_p2t, eta_c, \
m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, \
choked, T9_T5t, p9_p5t, eta_n, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, load_param, fuel_param, N  = \
engOperation(M0, N_c, nozzle, num_iter0, relaxation_factor)

if np.isnan(m_0):
    M0 = np.NaN
    N_c = np.NaN

# Print the variables on screen: --------------------------------------------------------------------------------------------------------------------------- 

print(" ")
print("     OPERATING POINT REPORT    ")
print("-------------------------------")
print(" ")
print("Input: ")
print("-------------------------------")
print("M0             [-] = " + str(np.round(M0,5)))
print("N*/Nref* (C) [-] = " + str(np.round(N_c,5)))
print(" ")
print("Output: ")
print("-------------------------------")

print("m* (0)   [kg/s] = " + str(np.round(m_0,5)))
print("m* (2t)  [kg/s] = " + str(np.round(m_2,5)))
print("m* (3t)  [kg/s] = " + str(np.round(m_3,5)))
print("m* (4t)  [kg/s] = " + str(np.round(m_4,5)))
print("m* (41t) [kg/s] = " + str(np.round(m_41,5)))
print("m* (5t)  [kg/s] = " + str(np.round(m_5,5)))

print("------------------------------")

print("T2t/T0    [-] = " + str(np.round(T2t_T0,5)))
print("T3t/T2t  [-] = " + str(np.round(T3t_T2t,5)))
print("T4t/T3t   [-] = " + str(np.round(T4t_T3t,5)))
print("T41t/T4t  [-] = " + str(np.round(T41t_T4t,5)))
print("T5t/T41t [-] = " + str(np.round(T5t_T41t,5)))
print("T9/T5t    [-] = " + str(np.round(T9_T5t,5)))

print("-------------------------------")

print("p2t/p0    [-] = " + str(np.round(T2t_T0,5)))
print("p3t/p2t  [-] = " + str(np.round(p3t_p2t,5)))
print("p4t/p3t   [-] = " + str(np.round(p4t_p3t,5)))
print("p41t/p4t  [-] = " + str(np.round(p41t_p4t,5)))
print("p5t/p41t [-] = " + str(np.round(p5t_p41t,5)))
print("p9/p5t    [-] = " + str(np.round(p9_p5t,5)))

print("-------------------------------")

print("ηd   [-] = " + str(np.round(eta_d,5)))
print("ηc [-] = " + str(np.round(eta_c,5)))
print("ηt [-] = " + str(np.round(eta_t,5)))
print("ηn   [-] = " + str(np.round(eta_n,5)))

print("-------------------------------")

print("T9/T0 [-] = " + str(np.round(T9_T0,5)))
print("p9/p0 [-] = " + str(np.round(p9_p0,5)))
print("M9    [-] = " + str(np.round(M9,5)))
print("A9/A8 [-] = " + str(np.round(A9_A8,5)))
print("Choke conditions: " + str(choked))

print("-------------------------------")

print("E*    [kN]      = " + str(np.round(E/1000,5)))
print("Isp*  [kN/kg/s] = " + str(np.round(Isp,5)))
print("TSFC* [m/s]     = " + str(np.round(TSFC,5)))

print("-------------------------------")

print("p3t/p2t (OPR)    [-] = " + str(np.round(p3t_p2t,5)))
print("T4t/T0  (OTR)    [-] = " + str(np.round(T4t_T3t*T3t_T2t*T2t_T0,5)))
print("T4t/T3t - 1      [-] = " + str(np.round(load_param,5)))
print("ηcc·f·L/(Cpc·T0) [-] = " + str(np.round(fuel_param,5)))

print("-------------------------------")

print("N*/Nref* (T) [-] = " + str(np.round(N_t,5)))

print("-------------------------------")

print("N* (C) [rpm] = " + str(np.round(N_c*N_ref_c*np.sqrt(T2t_T0),5)))
print("N* (T) [rpm] = " + str(np.round(N_t*N_ref_t*np.sqrt(T41t_T4t*T4t_T3t*T3t_T2t*T2t_T0),5)))

print(" ")

if plot:

    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker

    fig, ax = plt.subplots(1,2,num=1, figsize=(16,8), edgecolor='k')
    plt1 = plt.subplot(1,2,1)

    componentPlot("C",True,'viridis',0.5)
    plt1.scatter(m_2,p3t_p2t,125,marker="x",color='r',linewidth=2)
    plt1.set_title(r"$\bf{COMPRESSOR}$", fontsize = 12)

    ybox1 = TextArea(r"$M_{0} = $",textprops=dict(color='k',size=12,ha='center',va='bottom'))
    ybox2 = TextArea(str(np.round(M0,2)),textprops=dict(size=10,ha='center',va='bottom'))
    ybox3 = TextArea(r"$, \ \frac{N}{\sqrt{T_{\rm 0}/T_{\rm ref}}} = $",textprops=dict(size=14,ha='center',va='center'))
    ybox4 = TextArea(str(np.round(N,1)) + " rpm",textprops=dict(size=10,ha='center',va='bottom'))
    ybox = HPacker(children=[ybox1, ybox2, ybox3, ybox4],align="center",pad=0,sep=5)
    anchored_ybox = AnchoredOffsetbox(loc=10, child=ybox, pad=0, frameon=False, bbox_to_anchor=(1.06,1.06),\
    bbox_transform = ax[0].transAxes, borderpad=0)
    ax[0].add_artist(anchored_ybox)

    plt.suptitle(r"$\rm{OPERATING \ POINT \ - \ COMPRESSOR}$", fontsize = 14, weight = "bold")
    
    plt2 = plt.subplot(1,2,2)

    componentPlot("T",True,'viridis',0.5)
    plt2.scatter(m_41,1/p5t_p41t,125,marker="x",color='r',linewidth=2)
    plt2.set_title(r"$\bf{TURBINE}$", fontsize = 12)

    ybox1 = TextArea(r"$M_{0} = $",textprops=dict(color='k',size=12,ha='center',va='bottom'))
    ybox2 = TextArea(str(np.round(M0,2)),textprops=dict(size=10,ha='center',va='bottom'))
    ybox3 = TextArea(r"$, \ \frac{N}{\sqrt{T_{\rm 0}/T_{\rm ref}}} = $",textprops=dict(size=14,ha='center',va='center'))
    ybox4 = TextArea(str(np.round(N,1)) + " rpm",textprops=dict(size=10,ha='center',va='bottom'))
    ybox = HPacker(children=[ybox1, ybox2, ybox3, ybox4],align="center",pad=0,sep=5)
    anchored_ybox = AnchoredOffsetbox(loc=10, child=ybox, pad=0, frameon=False, bbox_to_anchor=(1.06,1.06),\
    bbox_transform = ax[0].transAxes, borderpad=0)
    ax[0].add_artist(anchored_ybox)

    plt.suptitle(r"$\rm{OPERATING \ POINT \ -  \ TURBINE}$", fontsize = 14, weight = "bold")
    plt.show()

