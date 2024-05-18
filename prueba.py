from Solvers.InternalCouplingSolver import lpCoupling
from Components.DesignVariables import N_ref_LPC, N_ref_HPC, N_ref_HPT,N_ref_LPT, M0_design
import warnings
from AuxilliaryFunctions.MapPlotFunction import componentPlot
import numpy as np

warnings.filterwarnings("ignore")

plot = True

# BETA_LPC DESIGN:

beta_LPC = 0.7687747928270632
N_LPC = 0.9949540661360499

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

m_2, T25t_T2t, p25t_p2t, eta_LPC, N_LPC_, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, fuel_param  = lpCoupling(beta_LPC, N_LPC, 15, 0.15, False)

# Print the variables on screen:

print(" ")
print("COUPLED POINT REPORT")
print(" ")

print("m* (2t)  [kg/s] = " + str(np.round(m_2,5)))
print("m* (25t) [kg/s] = " + str(np.round(m_25,5)))
print("m* (3t)  [kg/s] = " + str(np.round(m_3,5)))
print("m* (4t)  [kg/s] = " + str(np.round(m_4,5)))
print("m* (41t) [kg/s] = " + str(np.round(m_41,5)))
print("m* (45t) [kg/s] = " + str(np.round(m_45,5)))
print("m* (5t)  [kg/s] = " + str(np.round(m_5,5)))

print("----------------------------")

print("T25t/T2t  [-] = " + str(np.round(T25t_T2t,5)))
print("T3t/T25t  [-] = " + str(np.round(T3t_T25t,5)))
print("T4t/T3t   [-] = " + str(np.round(T4t_T3t,5)))
print("T41t/T4t  [-] = " + str(np.round(T41t_T4t,5)))
print("T45t/T41t [-] = " + str(np.round(T45t_T41t,5)))
print("T5t/T45t  [-] = " + str(np.round(T5t_T45t,5)))

print("----------------------------")

print("p25t/p2t  [-] = " + str(np.round(p25t_p2t,5)))
print("p3t/p25t  [-] = " + str(np.round(p3t_p25t,5)))
print("p4t/p3t   [-] = " + str(np.round(p4t_p3t,5)))
print("p41t/p4t   [-] = " + str(np.round(p41t_p4t,5)))
print("p45t/p41t [-] = " + str(np.round(p45t_p41t,5)))
print("p5t/p45t  [-] = " + str(np.round(p5t_p45t,5)))

print("----------------------------")

print("N*/Nref* (LPC) [-] = " + str(np.round(N_LPC,5)))
print("N*/Nref* (HPC) [-] = " + str(np.round(N_HPC,5)))
print("N*/Nref* (HPT) [-] = " + str(np.round(N_HPT,5)))
print("N*/Nref* (LPT) [-] = " + str(np.round(N_LPT,5)))

print("----------------------------")

T2t_T0 = 1 + (1.4-1)/2*M0_design**2

# Assuming T0_Tref = 1:

print("N (LPC) [rpm] = " + str(np.round(N_LPC*N_ref_LPC*np.sqrt(T2t_T0),5)))
print("N (HPC) [rpm] = " + str(np.round(N_HPC*N_ref_HPC*np.sqrt(T25t_T2t*T2t_T0),5)))
print("N (HPT) [rpm] = " + str(np.round(N_HPT*N_ref_HPT*np.sqrt(T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0),5)))
print("N (LPT) [rpm] = " + str(np.round(N_LPT*N_ref_LPT*np.sqrt(T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0),5)))

print("----------------------------")

print("(ηcc·f·L)/(Cp·T3t) [-] = " + str(np.round(fuel_param,5)))

print(" ")

if plot:

    pltLPC = componentPlot("LPC", False)
    pltLPC.scatter(m_2,p25t_p2t,50,marker="o",color='r',linewidth=1.25)
    pltLPC.show()

    pltHPC = componentPlot("HPC", False)
    pltHPC.scatter(m_25,p3t_p25t,50,marker="o",color='r',linewidth=1.25)
    pltHPC.show()

    pltHPT = componentPlot("HPT", False)
    pltHPT.scatter(m_41,1/p45t_p41t,50,marker="o",color='r',linewidth=1.25)
    pltHPT.show()

    pltLPT = componentPlot("LPT", False)
    pltLPT.scatter(m_45,1/p5t_p45t,50,marker="o",color='r',linewidth=1.25)
    pltLPT.show()

