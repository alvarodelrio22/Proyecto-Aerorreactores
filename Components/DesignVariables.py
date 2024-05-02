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

m_0_design = 80  # [kg/s]
M0_design = 0.8  # [-]

m_LPC_design = 56.18448
m_HPC_design = 23.51406
m_HPT_design = 14.27228
m_LPT_design = 22.71008

pi_LPC_design = 2.74000
pi_HPC_design = 3.28000
pi_HPT_design = 1.68481
pi_LPT_design = 1.41287

N_ref_LPC = 16499.99954
N_ref_HPC = 13216.44397
N_ref_HPT = 6882.020073
N_ref_LPT = 7735.030858

b_25 = 0.025
b_3 = 0.05

eta_mLP = 0.98
eta_mHP = 0.99

fuel_param_design = 1.56

NPR_design = 4.681519
A8 = 0.156473
A9 = 0.217675