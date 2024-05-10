from Components.ComponentMap import compressor, turbine
from Solvers.HighPressureCouplingSolver import hpCoupling
from Components.DesignVariables import gamma_c, gamma_e, Cp_c, Cp_e, N_ref_LPC, N_ref_LPT, b_25, eta_mLP, f_assumed
import numpy as np

# Using a Newton - Rhapson method with Succesive Over - Relaxation Method

def lpCoupling(beta_LPC,N_LPC,num_iter0,relaxation_factor,representing):

    max_iter = 10*num_iter0

## Diffuser Outlet - LPC Inlet (2t) -------------------------------------------------------------------------------------------------------------------------

    m_2 = compressor(beta_LPC,N_LPC,"beta","N","m","LPC")

## LPC Outlet - HPC Inlet (25t) -----------------------------------------------------------------------------------------------------------------------------
    
    p25t_p2t = compressor(beta_LPC,N_LPC,"beta","N","pi","LPC")
    eta_LPC = compressor(beta_LPC,N_LPC,"beta","N","eta","LPC")

    if np.isnan(m_2) or np.isnan(p25t_p2t) or np.isnan(eta_LPC):
        
        if representing:
        
            return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
            np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
                
        else:

            return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
            np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,        \
            np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,        \
            np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN  
            
    T25t_T2t = 1 + 1/eta_LPC*(p25t_p2t**((gamma_c-1)/gamma_c) - 1)
    m_25 = m_2*(1-b_25)*np.sqrt(T25t_T2t)/p25t_p2t

## Solution to the high pressure coupling problem - LPT Inlet (45t) -----------------------------------------------------------------------------------------
# An iterative process is requiered (using beta achieves better ortogonality than N):

    def hpIterations(beta_HPC, iterate):

        N_HPC = compressor(m_25, beta_HPC, "m", "beta", "N", "HPC")

        T3t_T25t, p3t_p25t, eta_HPC, m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, \
        m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, fuel_param = hpCoupling(beta_HPC, N_HPC, 100, 0.9, False)
        
        N_LPT = (N_LPC*N_ref_LPC)/np.sqrt(T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t)/N_ref_LPT

        p5t_p45t = 1/turbine(m_45,N_LPT,"m","N","pi","LPT")
        eta_LPT = turbine(m_45,N_LPT,"m","N","eta","LPT")
        T5t_T45t = 1 + eta_LPT*(p5t_p45t**((gamma_e-1)/gamma_e) - 1)

        m_5 = m_45*np.sqrt(T5t_T45t)/p5t_p45t

        if iterate:

            # Power balance (HPC-HPT) and solution
            # Determination of a new value for beta_HPC

            T25t_T2t_est = 1/(1 - eta_mLP*Cp_e/Cp_c*T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*(1 + f_assumed - b_25)*(1 - T5t_T45t))
            signed_error = T25t_T2t_est - T25t_T2t

            return signed_error
        
        else:

            return T3t_T25t, p3t_p25t, eta_HPC, N_HPC, m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, \
            m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, fuel_param  
    
    # Number of initial guesses is determined by the number of iterations:
    
    beta_HPC = 0         # Start from an initial value 0
    dBeta = 1e-5         # Increment for numerical derivative

    error = np.NaN
    iterations = 0

    while np.isnan(error) or error >= 1e-7:

        iterations = iterations + 1

        if iterations > max_iter:

            print("Limit time reached, no convergence: Modify the relaxation factor")

            if representing:
        
                return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
        
            else:

                return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,        \
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,        \
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN

        signed_error = hpIterations(beta_HPC,True)

        if np.isnan(signed_error) and iterations <= num_iter0:

            beta_HPC = iterations/num_iter0

        elif np.isnan(signed_error) and iterations > num_iter0:

            if representing:
            
                return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
        
            else:

                return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,        \
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,        \
                np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN 
        
        else:

            # Numerical derivative:

            dtau_dbeta =  (hpIterations(beta_HPC + dBeta,True) - signed_error)/dBeta

            # Newton - Rhapson method (second order) implementation:

            beta_HPC_star = beta_HPC - 1/dtau_dbeta*signed_error
            beta_HPC = beta_HPC_star*(1-relaxation_factor) + beta_HPC*relaxation_factor
            error = abs(signed_error)
        
    T3t_T25t, p3t_p25t, eta_HPC, N_HPC, m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, \
    m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, fuel_param = hpIterations(beta_HPC, False)

    if representing:

        beta_HPT = turbine(m_41,N_HPT,"m","N","beta","HPT")
        beta_LPT = turbine(m_45,N_LPT,"m","N","beta","LPT")

        return m_2, p25t_p2t, eta_LPC, m_25, p3t_p25t, eta_HPC, m_41, \
        p45t_p41t, eta_HPT, m_45, p5t_p45t, eta_LPT
        
    else:

        return m_2, T25t_T2t, p25t_p2t, eta_LPC, N_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
        m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, fuel_param  




