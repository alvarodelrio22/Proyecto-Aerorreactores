from Solvers.FunctioningPointSolver import engCoupling
from Miscellaneous.AuxilliaryFunctions import componentPlot
from Components.DesignVariables import N_ref_HPT, N_ref_LPT

import warnings, numpy as np, matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

# FUNCTIONING POINT ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

plot = True

# Choose entries: ------------------------------------------------------------------------------------------------------------------------------------------

M0 = 0.8
N_LPC = 0.9945
nozzle = 'conv'

# Loop parameters: -----------------------------------------------------------------------------------------------------------------------------------------

num_iter0 = 15
relaxation_factor = 0.15

m_0, T2t_T0, p2t_p0, eta_d, m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, \
eta_LPT, N_LPT, m_5, choked, T9_T5t, p9_p5t, eta_n, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, load_param, fuel_param, N1, N2  = \
engCoupling(M0, N_LPC, nozzle, num_iter0, relaxation_factor)

# Print the variables on screen: --------------------------------------------------------------------------------------------------------------------------- 

print(" ")
print("     COUPLED POINT REPORT    ")
print("-------------------------------")

print("m* (0)   [kg/s] = " + str(np.round(m_0,5)))
print("m* (2t)  [kg/s] = " + str(np.round(m_2,5)))
print("m* (25t) [kg/s] = " + str(np.round(m_25,5)))
print("m* (3t)  [kg/s] = " + str(np.round(m_3,5)))
print("m* (4t)  [kg/s] = " + str(np.round(m_4,5)))
print("m* (41t) [kg/s] = " + str(np.round(m_41,5)))
print("m* (45t) [kg/s] = " + str(np.round(m_45,5)))
print("m* (5t)  [kg/s] = " + str(np.round(m_5,5)))

print("------------------------------")

print("T2t/T0    [-] = " + str(np.round(T2t_T0,5)))
print("T25t/T2t  [-] = " + str(np.round(T25t_T2t,5)))
print("T3t/T25t  [-] = " + str(np.round(T3t_T25t,5)))
print("T4t/T3t   [-] = " + str(np.round(T4t_T3t,5)))
print("T41t/T4t  [-] = " + str(np.round(T41t_T4t,5)))
print("T45t/T41t [-] = " + str(np.round(T45t_T41t,5)))
print("T5t/T45t  [-] = " + str(np.round(T5t_T45t,5)))
print("T9/T5t    [-] = " + str(np.round(T9_T5t,5)))

print("-------------------------------")

print("p2t/p0    [-] = " + str(np.round(T2t_T0,5)))
print("p25t/p2t  [-] = " + str(np.round(p25t_p2t,5)))
print("p3t/p25t  [-] = " + str(np.round(p3t_p25t,5)))
print("p4t/p3t   [-] = " + str(np.round(p4t_p3t,5)))
print("p41t/p4t  [-] = " + str(np.round(p41t_p4t,5)))
print("p45t/p41t [-] = " + str(np.round(p45t_p41t,5)))
print("p5t/p45t  [-] = " + str(np.round(p5t_p45t,5)))
print("p9/p5t    [-] = " + str(np.round(p9_p5t,5)))

print("-------------------------------")

print("ηd   [-] = " + str(np.round(eta_d,5)))
print("ηLPC [-] = " + str(np.round(eta_LPC,5)))
print("ηHPC [-] = " + str(np.round(eta_HPC,5)))
print("ηHPT [-] = " + str(np.round(eta_HPT,5)))
print("ηLPT [-] = " + str(np.round(eta_LPT,5)))
print("ηn   [-] = " + str(np.round(eta_n,5)))

print("-------------------------------")

print("T9/T0 [-] = " + str(np.round(T9_T0,5)))
print("p9/p0 [-] = " + str(np.round(p9_p0,5)))
print("M9    [-] = " + str(np.round(M9,5)))
print("A9/A8 [-] = " + str(np.round(A9_A8,5)))
print("Choke conditions: " + str(choked))

print("-------------------------------")

print("E*    [kN]      = " + str(np.round(E,5)))
print("Isp*  [kN/kg/s] = " + str(np.round(Isp,5)))
print("TSFC* [m/s]     = " + str(np.round(TSFC,5)))

print("-------------------------------")

print("p3t/p2t (OPR)    [-] = " + str(np.round(p3t_p25t*p25t_p2t,5)))
print("T4t/T0  (OTR)    [-] = " + str(np.round(T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0,5)))
print("T4t/T3t - 1      [-] = " + str(np.round(load_param,5)))
print("ηcc·f·L/(Cpc·T0) [-] = " + str(np.round(fuel_param,5)))

print("-------------------------------")

print("N*/Nref* (LPC) [-] = " + str(np.round(N_LPC,5)))
print("N*/Nref* (HPC) [-] = " + str(np.round(N_HPC,5)))
print("N*/Nref* (HPT) [-] = " + str(np.round(N_HPT,5)))
print("N*/Nref* (LPT) [-] = " + str(np.round(N_LPT,5)))

print("-------------------------------")

print("NLP* (LPC) [rpm] = " + str(np.round(N1,5)))
print("NHP* (HPC) [rpm] = " + str(np.round(N2,5)))
print("NHP* (HPT) [rpm] = " + str(np.round(N_HPT*N_ref_HPT*np.sqrt(T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0),5)))
print("NLP* (LPT) [rpm] = " + str(np.round(N_LPT*N_ref_LPT*np.sqrt(T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0),5)))

print(" ")

if plot:

    fig = plt.figure(num=1, figsize=(14,8), edgecolor='k')
    plt1 = plt.subplot(1,2,1)

    componentPlot("LPC",True,'viridis',0.5)
    plt1.scatter(m_2,p25t_p2t,100,marker="x",color='r',linewidth=2)
    plt1.set_title(r"$\bf{LOW \ PRESSURE \ COMPRESSOR}$", fontsize = 12)

    plt2 = plt.subplot(1,2,2)

    componentPlot("HPC",True,'viridis',0.5)
    plt2.scatter(m_25,p3t_p25t,100,marker="x",color='r',linewidth=2)
    plt2.set_title(r"$\bf{HIGH \ PRESSURE \ COMPRESSOR}$", fontsize = 12)

    plt.suptitle(r"$\rm{FUNCTIONING \ POINT \ - \ COMPRESSORS}$", fontsize = 14, weight = "bold")
    plt.show()

    fig = plt.figure(num=2, figsize=(14,8), edgecolor='k')
    plt3 = plt.subplot(1,2,1)

    componentPlot("HPT",True,'viridis',0.5)
    plt3.scatter(m_41,1/p45t_p41t,100,marker="x",color='r',linewidth=2)
    plt3.set_title(r"$\bf{HIGH \ PRESSURE \ TURBINE}$", fontsize = 12)

    plt4 = plt.subplot(1,2,2)

    componentPlot("LPT",True,'viridis',0.5)
    plt4.scatter(m_45,1/p5t_p45t,100,marker="x",color='r',linewidth=2)
    plt4.set_title(r"$\bf{LOW \ PRESSURE \ TURBINE}$", fontsize = 12)

    plt.suptitle(r"$\rm{FUNCTIONING \ POINT \ - \ TURBINES}$", fontsize = 14, weight = "bold")
    plt.show()

