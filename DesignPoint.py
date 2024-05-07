from Components.ComponentMap import compressor, turbine
from scipy.optimize import newton
import numpy as np
import warnings

warnings.filterwarnings("ignore")

## Gas constants - Ambient conditions (Cold and Hot differentiation) ---------------------------------------------------------------------------------------

gamma_c = 1.4;    # [-]
gamma_e = 1.309;  # [-]
R = 287           # [J/kg/K]

Cp_c = 1004.5;    # [J/kg/K] (Room temperature)
Cp_e = 1210.0;    # [J/kg/K] (1500 K)

p_ref = 1.01325   # [bar]
T_ref = 288.15    # [K]

## Design parameters ---------------------------------------------------------------------------------------------------------------------------------------

m_0_design = 80  # [kg/s]
M0_design = 0.8  # [-]

pi_LPC_design = 2.74
pi_HPC_design = 3.28

N_ref_LPC = 972.0187*np.sqrt(T_ref)
N_ref_HPC = 778.5837*np.sqrt(T_ref)

b_25 = 0.025
b_3 = 0.05

eta_mLP = 0.98
eta_mHP = 0.99

fuel_param_design = 1.56

# Farfield conditions (0) -----------------------------------------------------------------------------------------------------------------------------------

m_0 = m_0_design
M0 = M0_design

p0 = p_ref
T0 = T_ref

# Diffuser outlet - LPC Inlet (2t) --------------------------------------------------------------------------------------------------------------------------

eta_d = 0.98

T2t_T0 = 1 + (gamma_c - 1)/2*M0_design**2
p2t_p0 = (1 + (gamma_c - 1)/2*eta_d*M0_design**2)**(gamma_c/(gamma_c - 1))
m_2 = m_0*np.sqrt(T2t_T0)/p2t_p0
m_LPC_design = m_2

# LPC Outlet - HPC Inlet (25t) ------------------------------------------------------------------------------------------------------------------------------

tolerance = 1e-8

def f(N):
    return abs(compressor(m_LPC_design, N, "m", "N", "pi", "LPC") - pi_LPC_design) - tolerance

N_LPC_design = newton(f,1)

eta_LPC = compressor(m_LPC_design, N_LPC_design, "m", "N", "eta", "LPC")
T25t_T2t = 1 + 1/eta_LPC*(pi_LPC_design**((gamma_c - 1)/gamma_c) - 1)
p25t_p2t = pi_LPC_design

m_25 = m_2*(1-b_25)*np.sqrt(T25t_T2t)/p25t_p2t
m_HPC_design = m_25

## HPC Outlet - Combustion Chamber Inlet (3t) ----------------------------------------------------------------------------------------------------------------

def f(N):
    return abs(compressor(m_HPC_design, N, "m", "N", "pi", "HPC") - pi_HPC_design)

N_HPC_design = newton(f,1)

eta_HPC = compressor(m_HPC_design, N_HPC_design, "m", "N", "eta","HPC")
T3t_T25t = 1 + 1/eta_HPC*(pi_HPC_design**((gamma_c - 1)/gamma_c) - 1)
p3t_p25t = pi_HPC_design

m_3 = m_25/(1-b_25)*(1-b_25-b_3)*np.sqrt(T3t_T25t)/p3t_p25t

## Combustion Chamber Outlet - NGV Bleed Injection (4t) -----------------------------------------------------------------------------------------------------

T4t_T3t = 1 + fuel_param_design
f_assumed = 0.025

PLF = 0.27 + 0.015*(T4t_T3t - 1)
p4t_p3t = 1 - PLF/2*(m_3*np.sqrt(R)*np.sqrt(T_ref)/(p_ref*1e5)/0.06)**2
    
m_4 =  m_3/(1-b_25-b_3)*(1+f_assumed-b_25-b_3)*np.sqrt(T4t_T3t)/p4t_p3t

## NGV Bleed Injection - HPT Inlet (41t) --------------------------------------------------------------------------------------------------------------------

p41t_p4t = 1
T41t_T4t = (m_4 + b_3*m_25/(1-b_25)*np.sqrt(T4t_T3t*T3t_T25t)*1/p3t_p25t*1/p4t_p3t*
Cp_c/Cp_e*(1/(T4t_T3t*T3t_T25t)))/(m_4 + b_3*m_25/(1-b_25)*np.sqrt(T4t_T3t*T3t_T25t)*1/p3t_p25t*1/p4t_p3t)
   
m_41 = m_4/(1+f_assumed-b_25-b_3)*(1+f_assumed-b_25)*np.sqrt(T41t_T4t)/p41t_p4t
m_HPT_design = m_41

# From power equilibrium and speed compatibility:

T41t_T25t = T41t_T4t*T4t_T3t*T3t_T25t
T45t_T41t = 1 - Cp_c*(1 - b_25)*(T3t_T25t - 1)/(eta_mHP*Cp_e*T41t_T25t*(1 + f_assumed - 0.5*b_3 - b_25))

N_HPT_design = 1
error = np.NaN

while np.isnan(error) or error >= 1e-8:

    eta_HPT = turbine(m_41, N_HPT_design, "m", "N", "eta", "HPT")
    pi_HPT_design = 1/(1 + 1/eta_HPT*(T45t_T41t - 1))**(gamma_e/(gamma_e - 1))

    def f(N):
        return abs(1+eta_HPT*(1/turbine(m_HPT_design, N, "m", "N", "pi", "HPT")**((gamma_e - 1)/gamma_e) - 1) - T45t_T41t)

    error = np.abs(newton(f,N_HPT_design) - N_HPT_design)
    
    N_HPT_design = newton(f,N_HPT_design)

N_ref_HPT = N_ref_HPC/np.sqrt(T41t_T25t)*N_HPC_design/N_HPT_design

p45t_p41t = 1/pi_HPT_design
m_45 = m_41*np.sqrt(T45t_T41t)/p45t_p41t
m_LPT_design = m_45

## HPT Outlet - LPT Inlet (45t) -----------------------------------------------------------------------------------------------------------------------------

# From power equilibrium and speed compatibility:

T45t_T2t = T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t
T5t_T45t = 1 - Cp_c*(T25t_T2t-1)/(eta_mLP*Cp_e*T45t_T2t*(1 + f_assumed - b_25))

N_LPT_design = 1
error = np.NaN

while np.isnan(error) or error >= 1e-8:

    eta_LPT = turbine(m_45, N_LPT_design, "m", "N", "eta", "LPT")
    pi_LPT_design = 1/(1 + 1/eta_LPT*(T5t_T45t - 1))**(gamma_e/(gamma_e - 1))

    def f(N):
        return abs(1+eta_LPT*(1/turbine(m_LPT_design, N, "m", "N", "pi", "LPT")**((gamma_e - 1)/gamma_e) - 1) - T5t_T45t)

    error = np.abs(newton(f,N_LPT_design) - N_LPT_design)
    
    N_LPT_design = newton(f,N_LPT_design)

N_ref_LPT = N_ref_LPC/np.sqrt(T45t_T2t)*N_LPC_design/N_LPT_design

p5t_p45t = 1/pi_LPT_design
m_5 = m_45*np.sqrt(T5t_T45t)/p5t_p45t

## LPT Outlet (5t) - Nozzle Outlet (9) ----------------------------------------------------------------------------------------------------------------------

# Designed to maximum expansion: choke conditions.

eta_n = 0.97

# As NPR > 5 a variable convergent - divergent segment is considered for the nozzle. A full expansion is considered for the design configuration.
# For NPR greater than 1.845, choke conditions happen for an engine embedded in ambient conditions. 
# No bigger NPR can be achieved without this phenomenon happening. It is assumed that the ambient conditions in the exit are p0, T0.

NPR_design = p5t_p45t*p45t_p41t*p41t_p4t*p4t_p3t*p3t_p25t*p25t_p2t*p2t_p0

def nozzle_eq(A8, M8):
    return A8*(p_ref*1e5)/np.sqrt(R*T_ref)*np.sqrt(gamma_e)*M8*(1+(1-1/eta_n)*(gamma_e-1)/2*M8**2)**(gamma_e/(gamma_e-1))* \
    (1+(gamma_e-1)/2*M8**2)**(-(gamma_e+1)/(2*(gamma_e-1)))

def f(A8):
    return np.abs(nozzle_eq(A8, 1) - m_5)

A8 = newton(f, 0.5)

p9_p0_design = 1
p9_p5t = p9_p0_design/NPR_design
M9 = np.sqrt(2/(gamma_e-1)*(1/(1-eta_n*(1-(p9_p5t)**((gamma_e-1)/gamma_e)))-1))

A9_A8 = (1+(1-1/eta_n)*(gamma_e-1)/2*M9**2)**(-gamma_e/(gamma_e-1))*(1/M9)*(2/(gamma_e+1)* \
(1+(gamma_e-1)/2*M9**2))**((gamma_e+1)/(2*(gamma_e-1)))

A9 = A9_A8*A8

# In case the design parameters want to be printed on screen:

print(" ")
print("DESIGN POINT REPORT")
print(" ")

print("m* (LPC) [kg/s] = " + str(m_LPC_design))
print("m* (HPC) [kg/s] = " + str(m_HPC_design))
print("m* (HPT) [kg/s] = " + str(m_HPT_design))
print("m* (LPT) [kg/s] = " + str(m_LPT_design))

print("------------------------------------------")

print("π (LPC) [-] = " + str(pi_LPC_design))
print("π (HPC) [-] = " + str(pi_HPC_design))
print("π (HPT) [-]] = " + str(pi_HPT_design))
print("π (LPT) [-] = " + str(pi_LPT_design))

print("------------------------------------------")

print("N*/Nref* (LPC) [-] = " + str(N_LPC_design))
print("N*/Nref* (HPC) [-] = " + str(N_HPC_design))
print("N*/Nref* (HPT) [-] = " + str(N_HPT_design))
print("N*/Nref* (LPT) [-] = " + str(N_LPT_design))

print("------------------------------------------")

print("ηcc·f·L/(Cp·T3t) [-] = " + str(fuel_param_design))

print("------------------------------------------")

print("A8 [m^2] = " + str(A8))
print("A9 [m^2] = " + str(A9))
print("NPR [-] = " + str(NPR_design))

print(" ")






