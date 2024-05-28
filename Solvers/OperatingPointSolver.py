
from Solvers.CouplingSolver import lpCoupling
from Miscellaneous.AuxilliaryFunctions import relaxationFactor
from Components.DesignVariables import p_ref, T_ref, gamma_c, gamma_e, Cp_c, Cp_e, f_assumed, b_25, b_3, R, A8, eta_n, eta_d, N_ref_LPC, N_ref_HPC
from scipy.optimize import newton

import time, numpy as np

## FUNCTIONING POINT SOLVER: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
# Iterative method (Newton - Rhapson with Succesive Over - Relaxation)
# Arguments: Flight Mach Number and LPC Relative Corrected Spool Speed

# Relaxation factor map for the LP coupling solver: --------------------------------------------------------------------------------------------------------

num_iter0 = 15     

relaxation_B_lim = [1, 0]
relaxation_N_lim = [0.45,0.7,1.08]
relaxation_matrix = [[0.15,0.00]]

## Function definition: ------------------------------------------------------------------------------------------------------------------------------------

def engOperation(M0, N_LPC, nozzle, num_iter0, relaxation_factor):

    # The intake is calculated according to the Mach Number:

    T2t_T0 = 1 + (gamma_c-1)/2*M0**2
    p2t_p0 = (1 + (gamma_c-1)/2*eta_d*M0**2)**(gamma_c/(gamma_c-1))

    ## Flight Mach loop ------------------------------------------------------------------------------------------------------------------------------------

    def engLoop(beta_LPC, iterate):

        m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
        m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, load_param, fuel_param_uncorrected = \
        lpCoupling(beta_LPC, N_LPC, num_iter0, relaxationFactor(beta_LPC,N_LPC,relaxation_B_lim,relaxation_N_lim,relaxation_matrix), False)

        m_0 = m_2/np.sqrt(T2t_T0)*p2t_p0
        fuel_param = T3t_T25t*T25t_T2t*T2t_T0*fuel_param_uncorrected        

        # The fuel parameter and shaft speeds are corrected with exterior conditions:

        N1 = N_LPC*N_ref_LPC*np.sqrt(T2t_T0)                                
        N2 = N_HPC*N_ref_HPC*np.sqrt(T25t_T2t*T2t_T0)

        # The maximum corrected mass flow according to A8 is calculated:

        m_5_max = (A8*p_ref*1e5)/np.sqrt(R*T_ref)*np.sqrt(gamma_e)*(1+(1-1/eta_n)*(gamma_e-1)/2)**(gamma_e/(gamma_e-1))* \
        (2/(gamma_e+1))**((gamma_e+1)/(2*(gamma_e-1)))                      

        if m_5 <= m_5_max:

            choked = False

            ## Nozzle equation:

            def nozzle_eq(M9):

                return A8*(p_ref*1e5)/np.sqrt(R*T_ref)*np.sqrt(gamma_e)*M9*(1+(1-1/eta_n)*(gamma_e-1)/2*M9**2)**(gamma_e/(gamma_e-1))* \
                (1+(gamma_e-1)/2*M9**2)**(-(gamma_e+1)/(2*(gamma_e-1))) - m_5

            # The nozzle is either sub-critical or critical:
            
            M9 = newton(nozzle_eq, 0.5)
            p9_p5t = (1-(1/eta_n)*(1-1/(1+(gamma_e-1)/2*M9**2)))**(gamma_e/(gamma_e-1))

            # Boundary condition (subsonic exit, Pascal's law):

            p9_p0 = 1

            # The intake is calculated according to the boundary condition:
            
            p2t_p0_est = p9_p0/(p9_p5t*p5t_p45t*p45t_p41t*p41t_p4t*p4t_p3t*p3t_p25t*p25t_p2t)
            signed_error = p2t_p0_est - p2t_p0

            if iterate:
                
                return signed_error
            
            else:

                # Exterior static temperature relation:

                T9_T5t = 1/(1+(gamma_e-1)/2*M9**2)
                T9_T0 = T9_T5t*T5t_T45t*T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0

                if nozzle == "conv":

                    A9_A8 = np.NaN

                elif nozzle == "conv-div":

                    A9_A8 = 1
                
                # Characteristic parameters:

                E = m_0*(((1+f_assumed)*(1-b_25-b_3)+b_3)*np.sqrt(gamma_e*R*T_ref)*M9*np.sqrt(T9_T0) - np.sqrt(gamma_c*R*T_ref)*M0) + \
                (p9_p0 - 1)*(p_ref*1e5*A8)
                Isp = E/m_0
                TSFC = ((1+f_assumed)*(Cp_e/Cp_c)*T4t_T3t -1)*(Cp_c*T_ref)*(T3t_T25t*T25t_T2t*T2t_T0)*(1-b_25-b_3)/Isp

                return m_0, T2t_T0, p2t_p0, eta_d, m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
                m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, \
                eta_LPT, N_LPT, m_5, choked, T9_T5t, p9_p5t, eta_n, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, load_param, fuel_param, N1, N2             

        elif np.isnan(m_5):

            choked = True

            if iterate:

                return np.NaN
            
            else:

                return np.NaN*np.empty(45)
            
        else:

            choked = True

            # The corrected mass flow in (5t) can not exceed the maximum value imposed by the throat:

            signed_error = m_5 - m_5_max

            if nozzle == "conv":

                T9_T5t = 2/(gamma_e+1)
                p9_p5t = (1-(1/eta_n)*(gamma_e-1)/(gamma_e+1))**(gamma_e/(gamma_e-1))
                M9 = 1

                # Exterior static conditions relation:

                T9_T0 = T9_T5t*T5t_T45t*T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0
                p9_p0 = p9_p5t*p5t_p45t*p45t_p41t*p41t_p4t*p4t_p3t*p3t_p25t*p25t_p2t*p2t_p0

                if p9_p0 < 1:

                    if iterate:

                        return np.NaN
                    
                    else:

                        return np.NaN*np.empty(45)
                    
                A9_A8 = np.NaN

                # Characteristic parameters:

                E = m_0*(((1+f_assumed)*(1-b_25-b_3)+b_3)*np.sqrt(gamma_e*R*T_ref)*M9*np.sqrt(T9_T0) - np.sqrt(gamma_c*R*T_ref)*M0) + \
                (p9_p0 - 1)*(p_ref*1e5*A8)
                Isp = E/m_0
                TSFC = ((1+f_assumed)*(Cp_e/Cp_c)*T4t_T3t -1)*(Cp_c*T_ref)*(T3t_T25t*T25t_T2t*T2t_T0)*(1-b_25-b_3)/Isp

            elif nozzle == "conv-div":

                # The nozzle area law is calculated to ensure a full expansion for all Mach numbers:

                p9_p0 = 1

                p9_p5t = p9_p0/(p5t_p45t*p45t_p41t*p41t_p4t*p4t_p3t*p3t_p25t*p25t_p2t*p2t_p0)
                M9 = np.sqrt(2/(gamma_e-1)*(1/(1-eta_n*(1-(p9_p5t)**((gamma_e-1)/gamma_e)))-1))
                
                # Exterior static temperature relation:
                
                if M9 < 1:

                    if iterate:

                        return np.NaN
                    
                    else:

                        return np.NaN*np.empty(45)
                    
                T9_T5t = 1/(1+(gamma_e-1)/2*M9**2)
                T9_T0 = T9_T5t*T5t_T45t*T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0

                # An extra nozzle efficiency between 8 and 9 is considered to be quadratic with M9 
                # due to a gradual apperture of the nozzle petals. Favorable pressure gradient.
                # No internal shockwaves due to underexpansion.

                eta_n_exit = 1 - 0.01*(M9**2-1)
                A9_A8 = (1+(1-1/eta_n_exit)*(gamma_e-1)/2*M9**2)**(-gamma_e/(gamma_e-1))*(1/M9)*(2/(gamma_e+1)* \
                (1+(gamma_e-1)/2*M9**2))**((gamma_e+1)/(2*(gamma_e-1)))

                A9 = A8*A9_A8

                # Characteristic parameters:

                E = m_0*(((1+f_assumed)*(1-b_25-b_3)+b_3)*np.sqrt(gamma_e*R*T_ref)*M9*np.sqrt(T9_T0) - np.sqrt(gamma_c*R*T_ref)*M0) + \
                (p9_p0 - 1)*(p_ref*1e5*A9)
                Isp = E/m_0
                TSFC = ((1+f_assumed)*(Cp_e/Cp_c)*T4t_T3t -1)*(Cp_c*T_ref)*(T3t_T25t*T25t_T2t*T2t_T0)*(1-b_25-b_3)/Isp

            else:

                print("Unknown type of nozzle. Options: 'conv' / 'conv-div'")
                print(" ")
                exit()

            if iterate:
                
                return signed_error
            
            else:

                return m_0, T2t_T0, p2t_p0, eta_d, m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
                m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, \
                eta_LPT, N_LPT, m_5, choked, T9_T5t, p9_p5t, eta_n, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, load_param, fuel_param, N1, N2            


    beta_LPC = 0        # Start from an initial value
    dBeta = 1e-7        # Increment for numerical derivative

    error = np.NaN
    tolerance = 1e-5
    iterations = 0

    start = time.time()
    timeOut = 25        # [s]

    while np.isnan(error) or error > tolerance:

        if time.time() - start > timeOut:

            print("(EN) -> Limit time reached, no convergence. Modify function parameters")
            print(" ")

            return np.NaN*np.empty(45)
            
        signed_error = engLoop(beta_LPC, True)
        error = abs(signed_error)

        if np.isnan(signed_error) and iterations < num_iter0:

            beta_LPC = iterations/(num_iter0-1)                 # Initialitation is critical due to the
            iterations = iterations + 1                         # potential existence of multiple solutions

        elif np.isnan(signed_error) and iterations == num_iter0:

            return np.NaN*np.empty(45) 
            
        elif ~np.isnan(signed_error) and error > tolerance:

            # Numerical derivative:

            dpi_dbeta =  (engLoop(beta_LPC + dBeta, True) - signed_error)/dBeta

            # Newton - Rhapson method (second order) implementation:

            beta_LPC_star = beta_LPC - 1/dpi_dbeta*signed_error
            beta_LPC = beta_LPC_star*(1-relaxation_factor) + beta_LPC*relaxation_factor
    
    return  engLoop(beta_LPC, False)