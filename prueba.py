from Solvers.HighPressureCouplingSolver import hpCoupling
from Solvers.LowPressureCouplingSolver import lpCoupling
from Components.ComponentMap import compressor, turbine
from Components.DesignVariables import gamma_c, gamma_e, Cp_c, Cp_e, N_ref_LPC, N_ref_LPT, b_25, eta_mLP, \
f_assumed, m_HPC_design, m_LPC_design, m_HPT_design, m_LPT_design, \
pi_LPC_design, pi_HPC_design, pi_HPT_design, pi_LPT_design 
import warnings
from matplotlib.pyplot import plot as plt
from Components.MapPlotFunction import componentPlot

warnings.filterwarnings("ignore")

# BETA_LPC DESIGN:

N_LPC_design = 0.9949538584005653
N_HPC_design = 0.9949538972546433

beta_LPC_design = compressor(m_LPC_design, N_LPC_design, "m", "N", "beta", "LPC")
beta_HPC_design = compressor(m_HPC_design, N_HPC_design, "m", "N", "beta", "HPC")

# m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, \
# m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, fuel_param = hpCoupling(beta_HPC_design, N_HPC_design, 100, 0.9, False)

# print(fuel_param)

m_2, T25t_T2t, p25t_p2t, eta_LPC, N_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
m_45, fuel_param, m_5, T5t_T45t, p5t_p45t, eta_LPT, N_LPT = lpCoupling(beta_LPC_design, N_LPC_design, 50, 0.65, False)

print(m_2, T25t_T2t, p25t_p2t, eta_LPC, N_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
m_45, fuel_param, m_5, T5t_T45t, p5t_p45t, eta_LPT, N_LPT)

pltLPC = componentPlot("LPC", False)
pltLPC.scatter(m_2,p25t_p2t,10,marker="*",linewidth=1.25)
pltLPC.show()

pltHPC = componentPlot("HPC", False)
pltHPC.scatter(m_25,p3t_p25t,10,marker="*",linewidth=1.25)
pltHPC.show()

pltHPT = componentPlot("HPT", False)
pltHPT.scatter(m_41,1/p45t_p41t,10,marker="*", linewidth=1.25)
pltHPT.show()

pltLPT = componentPlot("LPT", False)
pltLPT.scatter(m_45,1/p5t_p45t,10,marker="*",linewidth=1.25)
pltLPT.show()

