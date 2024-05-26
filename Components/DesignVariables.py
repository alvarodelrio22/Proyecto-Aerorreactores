from scipy.optimize import newton
import numpy as np

# Obtained from the calculations in DesignPoint.py

## Gas constants - Ambient conditions (Cold and Hot differentiation) ---------------------------------------------------------------------------------------

gamma_c = 1.4;    # [-]
gamma_e = 1.309;  # [-]
R = 287           # [J/kg/K]

Cp_c = 1004.5;    # [J/kg/K] (Room temperature)
Cp_e = 1210.0;    # [J/kg/K] (1500 K)

p_ref = 1.01325   # [bar]
T_ref = 288.15    # [K]

## Design parameters ---------------------------------------------------------------------------------------------------------------------------------------

m_0_design = 80.0         # [kg/s]
M0_design = 0.8           # [-]

m_LPC_design = 56.18447   # [kg/s]
m_HPC_design = 23.51406   # [kg/s]
m_HPT_design = 13.85085   # [kg/s]
m_LPT_design = 22.56845   # [kg/s]

pi_LPC_design = 2.74000   # [-]
pi_HPC_design = 3.28000   # [-]
pi_HPT_design = 1.68169   # [-]
pi_LPT_design = 1.41094   # [-]

N_ref_LPC = 16500.00000   # [rpm]
N_ref_HPC = 13216.40000   # [rpm]
N_ref_HPT = 6863.91362    # [rpm]
N_ref_LPT = 7713.00937    # [rpm]

b_25 = 0.025              # [-]
b_3 = 0.05                # [-]
momentum_factor = 0.5     # [-]
f_assumed = 0.025         # [-]

eta_d = 0.98              # [-]
eta_n = 0.97              # [-]
eta_mLP = 0.98            # [-]
eta_mHP = 0.99            # [-]

fuel_param_design = 1.56  # [-]
NPR_design = 5.63071      # [-]
A8 = 0.13277              # [m^2]
A9 = 0.20344              # [m^2]
