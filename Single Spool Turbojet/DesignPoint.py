from Components.ComponentMap import compressor, turbine
from scipy.optimize import newton
import numpy as np
import warnings

warnings.filterwarnings("ignore")

## DESIGN POINT DETERMINATION: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Iterating along component map definition

## Gas constants - Ambient conditions (Cold and Hot differentiation) ---------------------------------------------------------------------------------------

gamma_c = 1.4;    # [-]
gamma_e = 1.309;  # [-]
R = 287           # [J/kg/K]

Cp_c = 1004.5;    # [J/kg/K] (Reference Temperature)
Cp_e = 1210.0;    # [J/kg/K] (1500 K)

p_ref = 1.01325   # [bar]
T_ref = 288.15    # [K]

## Design parameters ----------------------------------------------------------------------------------------------------------------------------------------

m_0_design = 80            # [kg/s]
M0_design = 0.8            # [-]

pi_c_design = 10.52          # [-]
N_ref_c = 25000.0          # [rpm]

b_3 = 0.05                 # [-]
momentum_factor = 0.5      # [-]

eta_m = 0.99               # [-]
load_param_design = 1.56   # [-]

## Farfield conditions (0) ----------------------------------------------------------------------------------------------------------------------------------

m_0 = m_0_design
M0 = M0_design

p0 = p_ref
T0 = T_ref

## Diffuser Outlet - Compressor Inlet (2t) ------------------------------------------------------------------------------------------------------------------

eta_d = 0.98

T2t_T0 = 1 + (gamma_c - 1)/2*M0_design**2
p2t_p0 = (1 + (gamma_c - 1)/2*eta_d*M0_design**2)**(gamma_c/(gamma_c - 1))
m_2 = m_0*np.sqrt(T2t_T0)/p2t_p0

m_c_design = m_2

## Compressor Outlet - Combustion Chamber Inlet (3t) --------------------------------------------------------------------------------------------------------

tolerance = 1e-8

def f(N):
    return abs(compressor(m_c_design, N, "m", "N", "pi") - pi_c_design) - tolerance

N_c_design = newton(f,1)

eta_c = compressor(m_c_design, N_c_design, "m", "N", "eta")
T3t_T2t = 1 + 1/eta_c*(pi_c_design**((gamma_c - 1)/gamma_c) - 1)
p3t_p2t = pi_c_design

m_3 = m_2*(1-b_3)*np.sqrt(T3t_T2t)/p3t_p2t

## Combustion Chamber Outlet - NGV Bleed Injection (4t) -----------------------------------------------------------------------------------------------------

T4t_T3t = 1 + load_param_design
f_assumed = 0.025
A3 = 0.9  # [m^2]

PLF = 27 + 15*(T4t_T3t - 1)
p4t_p3t = 1 - PLF/2*(m_3*np.sqrt(R)*np.sqrt(T_ref)/(p_ref*1e5)/A3)**2

fuel_param_design = (T3t_T2t*T2t_T0)*(Cp_e/Cp_c*(1+f_assumed)*T4t_T3t - 1)
m_4 =  m_3*(1+f_assumed)*np.sqrt(T4t_T3t)/p4t_p3t

## NGV Bleed Injection - Turbine Inlet (41t) ----------------------------------------------------------------------------------------------------------------

p41t_p4t = 1
T41t_T4t = ((1+f_assumed)*(1-b_3) + b_3*(Cp_c/Cp_e)*(1/T4t_T3t))/((1+f_assumed)*(1-b_3) + b_3)
   
m_41 = m_4/((1+f_assumed)*(1-b_3))*((1+f_assumed)*(1-b_3)+momentum_factor*b_3)*np.sqrt(T41t_T4t)/p41t_p4t
m_t_design = m_41

# From power equilibrium and speed compatibility:

T5t_T41t = 1 - (T3t_T2t - 1)/(eta_m*(Cp_e/Cp_c)*(T41t_T4t*T4t_T3t*T3t_T2t)*((1+f_assumed)*(1-b_3)+momentum_factor*b_3))

N_t_design = 1
error = np.NaN

while np.isnan(error) or error >= 1e-8:

    eta_t = turbine(m_41, N_t_design, "m", "N", "eta")
    pi_t_design = 1/(1 + 1/eta_t*(T5t_T41t - 1))**(gamma_e/(gamma_e - 1))

    def f(N):
        return abs(1+eta_t*(1/turbine(m_t_design, N, "m", "N", "pi")**((gamma_e - 1)/gamma_e) - 1) - T5t_T41t)

    error = np.abs(newton(f,N_t_design) - N_t_design)
    
    N_t_design = newton(f,N_t_design)

N_ref_t = (N_ref_c*N_c_design)/np.sqrt(T41t_T4t*T4t_T3t*T3t_T2t)/N_t_design

## Turbine Outlet - Nozzle Inlet (5t) -----------------------------------------------------------------------------------------------------------------------

p5t_p41t = 1/pi_t_design
m_5 = m_41/((1+f_assumed)*(1-b_3)+momentum_factor*b_3)*((1+f_assumed)*(1-b_3)+b_3)*np.sqrt(T5t_T41t)/p5t_p41t
m_5_max = m_5

## Nozzle Outlet (9) ----------------------------------------------------------------------------------------------------------------------------------------

# Designed to maximum expansion: choke conditions.

eta_n = 0.97

# As NPR > 5 a variable convergent - divergent segment is considered for the nozzle. A full expansion is considered for the design configuration.
# For NPR greater than 1.845, choke conditions happen for an engine embedded in ambient conditions. 
# No bigger NPR can be achieved without this phenomenon happening. It is assumed that the ambient conditions in the exit are p0, T0.

NPR_design = p5t_p41t*p41t_p4t*p4t_p3t*p3t_p2t*p2t_p0

A8 = m_5/((p_ref*1e5)/np.sqrt(R*T_ref)*np.sqrt(gamma_e)*(1+(1-1/eta_n)*(gamma_e-1)/2)**(gamma_e/(gamma_e-1))* \
(2/(gamma_e+1))**((gamma_e+1)/(2*(gamma_e-1))))

p9_p0_design = 1
p9_p5t = p9_p0_design/NPR_design
M9 = np.sqrt(2/(gamma_e-1)*(1/(1-eta_n*(1-(p9_p5t)**((gamma_e-1)/gamma_e)))-1))

eta_n_exit = 1 - 0.01*(M9**2-1)
A9_A8 = (1+(1-1/eta_n_exit)*(gamma_e-1)/2*M9**2)**(-gamma_e/(gamma_e-1))*(1/M9)*(2/(gamma_e+1)* \
(1+(gamma_e-1)/2*M9**2))**((gamma_e+1)/(2*(gamma_e-1)))

A9 = A9_A8*A8

# Print the design parameters on screen:

print(" ")
print("DESIGN POINT REPORT")
print(" ")

print("m* (C) [kg/s] = " + str(np.round(m_c_design,5)))
print("m* (T) [kg/s] = " + str(np.round(m_t_design,5)))
print("m* (5t*) [kg/s] = " + str(np.round(m_5_max,5)))

print("------------------------------")

print("π (C) [-] = " + str(np.round(pi_c_design,5)))
print("π (T) [-] = " + str(np.round(pi_t_design,5)))

print("------------------------------")

print("N*/Nref* (C) [-] = " + str(np.round(N_c_design,5)))
print("N*/Nref* (T) [-] = " + str(np.round(N_t_design,5)))

print("------------------------------")

print("Nref* (C) [rpm] = " + str(np.round(N_ref_c,5)))
print("Nref* (T) [rpm] = " + str(np.round(N_ref_t,5)))

print("------------------------------")

print("T4t/T3t - 1 [-] = " + str(np.round(load_param_design,5)))
print("ηcc·f·L/(Cpc·T0) [-] = " + str(np.round(fuel_param_design,5)))

print("------------------------------")

print("A8  [m^2] = " + str(np.round(A8,5)))
print("A9  [m^2] = " + str(np.round(A9,5)))
print("NPR [-] = " + str(np.round(NPR_design,5)))

print(" ")





