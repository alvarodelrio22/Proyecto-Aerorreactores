from Components.ComponentMap import compressor, turbine
from Components.DesignVariables import p_ref, T_ref, gamma_c, gamma_e, Cp_c, Cp_e, eta_mHP, eta_mLP, \
f_assumed, b_25, b_3, momentum_factor, R, A8, N_ref_LPC, N_ref_HPC, N_ref_HPT, N_ref_LPT, eta_d, eta_n
from scipy.optimize import newton

import time, numpy as np

## TRANSIENT OPERATION SOLVER :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Iterative method (Newton - Rhapson with Succesive Over - Relaxation) + Initial Value Problem Numerical Solver (with choice of propagator)
# Arguments: Initial State Conditions and Corrected Fuel Parameter Time Evolution

## Function definition: -------------------------------------------------------------------------------------------------------------------------------------

def transientOperation(w_0, time_fuel_param, propagator, delta_t, relaxation_factor, num_iter0, tolerance, nozzle):

    # Constants and Initial Conditions: 

    M0, I1, I2, N1_0, N2_0 = w_0

    # Operation Function Definition: 

    def transientOperationPoint(beta0, N1, N2, fuel_param, num_iter0, relaxation_factor, tolerance):

## Farfield Conditions (0) ----------------------------------------------------------------------------------------------------------------------------------

        T2t_T0 = 1 + (gamma_c-1)/2*M0**2
        p2t_p0 = (1 + (gamma_c-1)/2*eta_d*M0**2)**(gamma_c/(gamma_c-1))

        ## Flight Mach loop -------------------------------------------------------------------------------------------------------------------------------------

        def transientLoop(beta_LPC, iterate):

            # Shaft speeds are corrected with inlet conditions:

            N_LPC = N1/N_ref_LPC/np.sqrt(T2t_T0)                          

## Diffuser Outlet - LPC Inlet (2t) -------------------------------------------------------------------------------------------------------------------------

            m_2 = compressor(beta_LPC,N_LPC,"beta","N","m","LPC")

            # Farfield corrected mass flow is calculated:

            m_0 = m_2/np.sqrt(T2t_T0)*p2t_p0

## LPC Outlet - HPC Inlet (25t) -----------------------------------------------------------------------------------------------------------------------------

            p25t_p2t = compressor(beta_LPC,N_LPC,"beta","N","pi","LPC")
            eta_LPC = compressor(beta_LPC,N_LPC,"beta","N","eta","LPC") 
                
            T25t_T2t = 1 + 1/eta_LPC*(p25t_p2t**((gamma_c-1)/gamma_c) - 1)
            m_25 = m_2*(1-b_25)*np.sqrt(T25t_T2t)/p25t_p2t

            N_HPC = N2/N_ref_HPC/np.sqrt(T25t_T2t*T2t_T0)

## HPC Outlet - Combustion Chamber Inlet (3t) ---------------------------------------------------------------------------------------------------------------

            p3t_p25t = compressor(m_25,N_HPC,"m","N","pi","HPC")
            eta_HPC = compressor(m_25,N_HPC,"m","N","eta","HPC")

            T3t_T25t = 1 + 1/eta_HPC*(p3t_p25t**((gamma_c-1)/gamma_c)-1)
            m_3 = m_25/(1-b_25)*(1-b_25-b_3)*np.sqrt(T3t_T25t)/p3t_p25t

## Combustion Chamber Outlet - NGV Bleed Injection (4t) -----------------------------------------------------------------------------------------------------

            T4t_T3t = (Cp_c/Cp_e)/(1+f_assumed)*(fuel_param/(T3t_T25t*T25t_T2t*T2t_T0) + 1)

            A3 = 0.9  # [m^2]
            PLF = 27 + 15*(T4t_T3t - 1)
            p4t_p3t = 1 - PLF/2*(m_3*np.sqrt(R)*np.sqrt(T_ref)/(p_ref*1e5)/A3)**2
    
            m_4 =  m_3*(1+f_assumed)*np.sqrt(T4t_T3t)/p4t_p3t

## Combustion Chamber Outlet - NGV Bleed Injection (4t) -----------------------------------------------------------------------------------------------------

            p41t_p4t = 1
            T41t_T4t = ((1+f_assumed)*(1-b_25-b_3) + b_3*(Cp_c/Cp_e)*(1/T4t_T3t))/((1+f_assumed)*(1-b_25-b_3) + b_3)

            m_41 = m_4/((1+f_assumed)*(1-b_25-b_3))*((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)*np.sqrt(T41t_T4t)/p41t_p4t

## HPT Outlet - LPT Inlet (45t) -----------------------------------------------------------------------------------------------------------------------------

            N_HPT = (N_HPC*N_ref_HPC)/np.sqrt(T41t_T4t*T4t_T3t*T3t_T25t)/N_ref_HPT
            
            p45t_p41t = 1/turbine(m_41,N_HPT,"m","N","pi","HPT")
            eta_HPT = turbine(m_41,N_HPT,"m","N","eta","HPT")
    
            T45t_T41t = 1 + eta_HPT*((p45t_p41t)**((gamma_e-1)/gamma_e) - 1)
            m_45 = m_41/((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)*((1+f_assumed)*(1-b_25-b_3)+b_3)*np.sqrt(T45t_T41t)/p45t_p41t

## LPT Outlet - Nozzle Inlet (5t) ---------------------------------------------------------------------------------------------------------------------------
   
            N_LPT = (N_LPC*N_ref_LPC)/np.sqrt(T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t)/N_ref_LPT

            p5t_p45t = 1/turbine(m_45,N_LPT,"m","N","pi","LPT")
            eta_LPT = turbine(m_45,N_LPT,"m","N","eta","LPT")
            T5t_T45t = 1 + eta_LPT*(p5t_p45t**((gamma_e-1)/gamma_e) - 1)

            m_5 = m_45*np.sqrt(T5t_T45t)/p5t_p45t

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
                signed_error = (p2t_p0_est - p2t_p0)/p2t_p0

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

                    return m_0, T2t_T0, p2t_p0, m_2, T25t_T2t, p25t_p2t, eta_LPC, N_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
                    m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, \
                    eta_LPT, N_LPT, m_5, choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC             

            elif np.isnan(m_5):

                choked = True

                if iterate:
                    
                    return np.NaN
                
                else:

                    return np.NaN*np.empty(40)
                
            else:

                choked = True

                # The corrected mass flow in (5t) can not exceed the maximum value imposed by the throat:

                signed_error = (m_5 - m_5_max)/m_5_max

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

                            return np.NaN*np.empty(40)
                        
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

                            return np.NaN*np.empty(40)
                        
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

                    return m_0, T2t_T0, p2t_p0, m_2, T25t_T2t, p25t_p2t, eta_LPC, N_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
                    m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, \
                    eta_LPT, N_LPT, m_5, choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC


        beta_LPC = beta0    # Start from an initial value
        dBeta = 1e-8        # Increment for numerical derivatives

        error = np.NaN
        iterations = 0

        start = time.time()
        timeOut = 15       # [s]

        while np.isnan(error) or error > tolerance:

            if time.time() - start > timeOut:

                print("(TransEN) -> Limit time reached, no convergence. Modify function parameters")

                return np.NaN*np.empty(40)
            
            elif time.time() - start > timeOut/2 and error < 5e-2:

                print("(TransEN) -> Warning. Low rate of convergence. Assumed low error restarting point.")

                break
        
            signed_error = transientLoop(beta_LPC, True)
            error = abs(signed_error)

            if np.isnan(signed_error) and iterations < num_iter0:
                
                beta_LPC = iterations/(num_iter0-1)        
                iterations = iterations + 1                      

            elif np.isnan(signed_error) and iterations == num_iter0:

                return np.NaN*np.empty(40) 
                
            elif ~np.isnan(signed_error) and error > tolerance:
                
                # Functions of numerical derivatives:
                
                F = dBeta*(signed_error - transientLoop(beta_LPC - dBeta, True))/(signed_error - 2*transientLoop(beta_LPC - dBeta, True) + \
                transientLoop(beta_LPC - 2*dBeta, True))
                G = (signed_error*dBeta**2)/(signed_error - 2*transientLoop(beta_LPC - dBeta, True) + \
                transientLoop(beta_LPC - 2*dBeta, True))
            
                # Halley's method (cubic order of convergence):

                beta_LPC_star = beta_LPC + np.min([-F-np.sqrt(F**2-2*G),-F+np.sqrt(F**2-2*G)])
                beta_LPC = beta_LPC_star*(1-relaxation_factor) + beta_LPC*relaxation_factor

        return  transientLoop(beta_LPC, False)

    # Initial Value Problem Solver:
    
    N1 = N1_0
    N2 = N2_0

    n = np.size(time_fuel_param)
    iteration, beta0 = 0, 0.8
    w = np.empty([n,40])
    time_N1, time_N2, current_time = [], [], []

    while iteration in range(n):

        # Update fuel parameter and time:

        fuel_param = time_fuel_param[iteration]
        current_time.append(iteration*delta_t)

        m_0, T2t_T0, p2t_p0, m_2, T25t_T2t, p25t_p2t, eta_LPC, N_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, \
        eta_LPT, N_LPT, m_5, choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC = \
        transientOperationPoint(beta0, N1, N2, fuel_param, num_iter0, relaxation_factor[iteration], tolerance[iteration])

        w[iteration,:] = [m_0, T2t_T0, p2t_p0, m_2, T25t_T2t, p25t_p2t, eta_LPC, N_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, \
        eta_LPT, N_LPT, m_5, choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC]

        # Coupling Excess Power:

        if propagator == 'Euler':

            N1 = N1 + delta_t*(m_0*Cp_c*T2t_T0*T_ref*((Cp_e/Cp_c)*eta_mLP*((1+f_assumed)*(1-b_25-b_3)+b_3)*(1-T5t_T45t)*\
            T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t - (T25t_T2t-1)))/((2*np.pi/60)**2*I1*N1)

            time_N1.append(N1)

            N2 = N2 + delta_t*(m_0*Cp_c*T25t_T2t*T2t_T0*T_ref*((Cp_e/Cp_c)*eta_mHP*((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)*(1-T45t_T41t)*\
            T41t_T4t*T4t_T3t*T3t_T25t - (1-b_25)*(T3t_T25t-1)))/((2*np.pi/60)**2*I2*N2)

            time_N2.append(N2)

        # Update the iteration counter:

        iteration = iteration + 1
        beta0 = compressor(p25t_p2t,N_LPC,"pi","N","beta","LPC")    # Initialize according to last point

        print("Iteration: "+str(iteration)+"/"+str(n)+": "+str(np.round(iteration/n*100,2))+"%")

    return w, time_N1, time_N2, current_time
