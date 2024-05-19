
from Solvers.InternalCouplingSolver import lpCoupling
from Components.DesignVariables import p_ref, T_ref, gamma_c, gamma_e, Cp_c, Cp_e, R, A8, eta_n, eta_d, m_5_max
from scipy.optimize import newton

import time, numpy as np
            
# Solver for non-choke conditions given Flight Mach Number -------------------------------------------------------------------------------------------------

def nonChokedExternalCoupling(M0, N_LPC, num_iter0, relaxation_factor, representing):

    # The intake is calculated according to the Mach Number:

    T2t_T0 = 1 + (gamma_c-1)/2*M0**2
    p2t_p0 = (1 + (gamma_c-1)/2*eta_d*M0**2)**(gamma_c/(gamma_c-1))

    ## Flight Mach loop ----------------------------------------------------------------------------------------------------------------------------------------

    def M0Loop(beta_LPC, iterate):

        m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
        m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, fuel_param  = lpCoupling(beta_LPC, N_LPC, 25, 0.45, False)

        if np.isnan(m_2) or np.isnan(p25t_p2t) or np.isnan(T25t_T2t):

            if iterate:

                return np.NaN
            
            else:

                return np.NaN*np.empty(30)

        if m_5 <= m_5_max:

            # Boundary condition (subsonic exit, Pascal's law):

            p8_p0 = 1

            ## Nozzle equation ----------------------------------------------------------------------------------------------------------------------------------------

            def nozzle_eq(M8):

                return np.abs(A8*(p_ref*1e5)/np.sqrt(R*T_ref)*np.sqrt(gamma_e)*M8*(1+(1-1/eta_n)*(gamma_e-1)/2*M8**2)**(gamma_e/(gamma_e-1))* \
                (1+(gamma_e-1)/2*M8**2)**(-(gamma_e+1)/(2*(gamma_e-1))) - m_5)

            # The nozzle is not critical:
            
            M8 = newton(nozzle_eq, 0.5)
            p8_p5t = (1-1/eta_n*(1-1/(1+(gamma_e-1)/2*M8**2)))**(gamma_e/(gamma_e-1))
            T8_T5t = 1/(1+(gamma_e-1)/2*M8**2)

            # The intake is calculated according to the boundary condition:
            
            p2t_p0_est = p8_p0/(p8_p5t*p5t_p45t*p45t_p41t*p41t_p4t*p4t_p3t*p3t_p25t*p25t_p2t)

            if p2t_p0_est < 1:

                if iterate:

                    return np.NaN
                
                else:

                    return np.NaN*np.empty(30)
            
            else:
            
                signed_error = p2t_p0_est - p2t_p0

                if iterate:
                
                    return signed_error
            
                else:

                    return m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
                    m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
                    m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, T8_T5t, p8_p5t, M8, fuel_param              

        else:

            if iterate:

                return np.NaN
            
            else:

                return np.NaN*np.empty(30)

    beta_LPC = 0.75   # First initial condition similar to beta in design point
    dBeta = 1e-7      # Increment for numerical derivative

    error = np.NaN
    tolerance = 1e-5
    iterations = 0

    start = time.time()
    timeOut = 60

    while np.isnan(error) or error > tolerance:

        if time.time() - start > timeOut:

            print("(EXT) ---> Limit time reached, no convergence: Modify the relaxation factor")
            print(" ")

            if representing:

                return np.NaN*np.empty(12) 
            
            else:

                return np.NaN*np.empty(37)
            
        signed_error = M0Loop(beta_LPC, True)
        error = abs(signed_error)

        if np.isnan(signed_error) and iterations < num_iter0:

            beta_LPC = 1 - iterations/(num_iter0-1)             # Initialitation is critical due to the
            iterations = iterations + 1                         # potential existence of multiple solutions
            continue

        elif np.isnan(signed_error) and iterations == num_iter0:

            if representing:

                return np.NaN*np.empty(12) 
        
            else:

                return np.NaN*np.empty(37) 
            
        elif ~np.isnan(signed_error) and error > tolerance:

            # Numerical derivative:

            dpi_dbeta =  (M0Loop(beta_LPC + dBeta, True) - signed_error)/dBeta

            # Newton - Rhapson method (second order) implementation:

            beta_LPC_star = beta_LPC - 1/dpi_dbeta*signed_error
            beta_LPC = beta_LPC_star*(1-relaxation_factor) + beta_LPC*relaxation_factor
    
    m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
    m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
    m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, T8_T5t, p8_p5t, M8, fuel_param = M0Loop(beta_LPC, False)

    m_0 = m_2/np.sqrt(T2t_T0)*p2t_p0

    if representing:

        return m_2, p25t_p2t, eta_LPC, m_25, p3t_p25t, eta_HPC, m_41, \
        p45t_p41t, eta_HPT, m_45, p5t_p45t, eta_LPT
        
    else:

        return  m_0, T2t_T0, p2t_p0, eta_d, M0, m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
        m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, T8_T5t, p8_p5t, eta_n, M8, fuel_param