from Components.ComponentMap import compressor, turbine
from Solvers.HighPressureCouplingSolver import hpCoupling
from Components.DesignVariables import gamma_c, gamma_e, Cp_c, Cp_e, N_ref_LPC, N_ref_LPT, b_25, eta_mLP as eta_m
import numpy as np

def lpCoupling(beta_LPC,N_LPC,num_iter,relaxation_factor,representing):

    max_iterations = 10*num_iter

## Diffuser Outlet - LPC Inlet (2t) -------------------------------------------------------------------------------------------------------------------------

    m_2 = compressor(beta_LPC,N_LPC,"beta","N","m","LPC")

## LPC Outlet - HPC Inlet (25t) -----------------------------------------------------------------------------------------------------------------------------
    
    p25t_p2t = compressor(beta_LPC,N_LPC,"beta","N","pi","LPC")
    eta_LPC = compressor(beta_LPC,N_LPC,"beta","N","eta","LPC")

    if np.isnan(m_2) or np.isnan(p25t_p2t) or np.isnan(eta_LPC):
        
        if representing == True:
            return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
        else:
            return np.NaN, np.NaN, np.NaN
            
    T25t_T2t = 1 + 1/eta_LPC*(p25t_p2t**((gamma_c-1)/gamma_c) - 1)
    m_25 = m_2*(1-b_25)*np.sqrt(T25t_T2t)/p25t_p2t
    f_assumed = 0.025

## Solution to the high pressure coupling problem - LPT Inlet (45t) -----------------------------------------------------------------------------------------
# An iterative process is requiered (using beta achieves better ortogonality than N):

    beta_HPC = 0      # Initial guess starts at zero
    error = np.NaN
    iterations = 0

    while np.isnan(error) or error >= 1e-3:

        iterations = iterations + 1

        if iterations > max_iterations:

            print("Limit time reached, no convergence: Modify the relaxation factor")
    
            if representing == True:
                return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
            else:
                return np.NaN, np.NaN, np.NaN

        N_HPC = compressor(m_25, beta_HPC, "m", "beta", "N", "HPC")

#        if np.isnan(N_HPC):

#            beta_HPC = iterations/max_iterations
#            continue

#        else:

        m_25, p3t_p25t, eta_HPC, m_41, p45t_p41t, eta_HPT, m_45, T45t_T25t, p45t_p25t, fuel_param = hpCoupling(beta_HPC, N_HPC, 100, 0.9, False)

#        if np.isnan(m_45) or np.isnan(T45t_T25t):

#            beta_HPC = iterations/max_iterations
#            continue

#        else:

        N_LPT = (N_LPC*N_ref_LPC)/np.sqrt(T45t_T25t)/np.sqrt(T25t_T2t)/N_ref_LPT

        p5t_p45t = 1/turbine(m_45,N_LPT,"m","N","pi","LPT")
        eta_LPT = turbine(m_45,N_LPT,"m","N","eta","LPT")

        T5t_T45t = 1 + eta_LPT*(p5t_p45t**((gamma_e-1)/gamma_e) - 1)
        m_5 = m_45*np.sqrt(T5t_T45t)/p5t_p45t

            ## Power balance (HPC-HPT) and solution
            # Determination of a new value for beta_HPC

        if iterations <= num_iter:

            if np.isnan(T5t_T45t) or np.isnan(T45t_T25t):

                beta_HPC = iterations/num_iter

            else:

                # current_iterations = iterations - 1

                error = abs(T25t_T2t - (1/(1 - eta_m*Cp_e/Cp_c*T45t_T25t*(1 + f_assumed - b_25)*(1 - T5t_T45t))))

                if error > 1e-1:
                    progression_factor = 5e-2
                elif error > 1e-2 and error <= 1e-1:
                    progression_factor = 2.5e-2
                elif error <= 1e-2 :
                    progression_factor = 5e-3

                if T25t_T2t > (1/(1 - eta_m*Cp_e/Cp_c*T45t_T25t*(1 + f_assumed - b_25)*(1 - T5t_T45t))):

                    beta_HPC = beta_HPC + progression_factor*(1-relaxation_factor)

                else: 

                    beta_HPC = beta_HPC - progression_factor*(1-relaxation_factor)             

                print(error)

        elif iterations > num_iter:
            
            if np.isnan(T5t_T45t) or np.isnan(T45t_T25t):

                if representing == True:
            
                    return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
                
                else:

                    return np.NaN, np.NaN, np.NaN, np.NaN
            
            else:

                error = abs(T25t_T2t - (1/(1 - eta_m*Cp_e/Cp_c*T45t_T25t*(1 + f_assumed - b_25)*(1 - T5t_T45t))))

                if error > 1e-1:
                    progression_factor = 5e-2
                elif error > 1e-2 and error <= 1e-1:
                    progression_factor = 2.5e-2
                elif error <= 1e-2:
                    progression_factor = 5e-3

                if T25t_T2t > (1/(1 - eta_m*Cp_e/Cp_c*T45t_T25t*(1 + f_assumed - b_25)*(1 - T5t_T45t))):

                    beta_HPC = beta_HPC + progression_factor*(1-relaxation_factor)

                else: 

                    beta_HPC = beta_HPC - progression_factor*(1-relaxation_factor) 

                print(error)
          
    if representing == True:
        
        return m_2, m_25, p25t_p2t, p3t_p25t, eta_LPC, eta_HPC, m_41, m_45, p45t_p41t, p5t_p45t, eta_LPT, eta_HPT
        
    else:

        return fuel_param, m_5, T5t_T45t, p5t_p45t   

