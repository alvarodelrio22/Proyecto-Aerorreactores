from Components.ComponentMap import compressor, turbine
from Components.DesignVariables import p_ref, T_ref, gamma_c, gamma_e, Cp_c, Cp_e, eta_m, \
f_assumed, b_3, momentum_factor, R, A8, N_ref_c, N_ref_t, eta_d
from scipy.optimize import newton

import time, numpy as np, csv

## TRANSIENT OPERATION SOLVER :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Iterative method (Newton - Rhapson with Succesive Over - Relaxation) + Initial Value Problem Numerical Solver (with choice of propagator)
# Arguments: Initial State Conditions and Corrected Fuel Parameter Time Evolution

## Function definition: -------------------------------------------------------------------------------------------------------------------------------------

def transientOperation(w_0, time_fuel_param, propagator, delta_t, relaxation_factor, num_iter0, tolerance, nozzle):

    # Constants and Initial Conditions: 

    M0, I, N_0 = w_0

    # Operation Function Definition: 

    def transientOperationPoint(beta0, N, fuel_param, num_iter0, relaxation_factor, tolerance):

## Farfield Conditions (0) ----------------------------------------------------------------------------------------------------------------------------------

        T2t_T0 = 1 + (gamma_c-1)/2*M0**2
        p2t_p0 = (1 + (gamma_c-1)/2*eta_d*M0**2)**(gamma_c/(gamma_c-1))

        ## Flight Mach loop -------------------------------------------------------------------------------------------------------------------------------------

        def transientLoop(beta_c, iterate):

            # Shaft speeds are corrected with inlet conditions:

            N_c = N/N_ref_c/np.sqrt(T2t_T0)                          

## Diffuser Outlet - Compressor Inlet (2t) ------------------------------------------------------------------------------------------------------------------

            m_2 = compressor(beta_c,N_c,"beta","N","m")

            # Farfield corrected mass flow is calculated:

            m_0 = m_2/np.sqrt(T2t_T0)*p2t_p0

## Compressor Outlet -  Combustion Chamber Inlet (3t)  ------------------------------------------------------------------------------------------------------

            p3t_p2t = compressor(beta_c,N_c,"beta","N","pi")
            eta_c = compressor(beta_c,N_c,"beta","N","eta") 
                
            T3t_T2t = 1 + 1/eta_c*(p3t_p2t**((gamma_c-1)/gamma_c) - 1)
            m_3 = m_2*(1-b_3)*np.sqrt(T3t_T2t)/p3t_p2t

## Combustion Chamber Outlet - NGV Bleed Injection (4t) -----------------------------------------------------------------------------------------------------

            T4t_T3t = (Cp_c/Cp_e)/(1+f_assumed)*(fuel_param/(T3t_T2t*T2t_T0*m_0) + 1)

            A3 = 0.0025  # [m^2]
            PLF = 1/3 + 5/27*(T4t_T3t - 1)
            p4t_p3t = 1 - PLF/2*(m_2*(1-b_3)*np.sqrt(R)*np.sqrt(T_ref)/(p_ref*1e5)/A3)**2
    
            m_4 =  m_3*(1+f_assumed)*np.sqrt(T4t_T3t)/p4t_p3t

## NGV Bleed Injection - Turbine Inlet (41t) ----------------------------------------------------------------------------------------------------------------

            p41t_p4t = 1
            T41t_T4t = ((1+f_assumed)*(1-b_3) + b_3*(Cp_c/Cp_e)*(1/T4t_T3t))/((1+f_assumed)*(1-b_3) + b_3)
    
            m_41 = m_4/((1+f_assumed)*(1-b_3))*((1+f_assumed)*(1-b_3)+momentum_factor*b_3)*np.sqrt(T41t_T4t)/p41t_p4t

## Turbine Outlet - Nozzle Inlet (5t) -----------------------------------------------------------------------------------------------------------------------
   
            N_t = N_c*N_ref_c/np.sqrt(T41t_T4t*T4t_T3t*T3t_T2t)/N_ref_t
            
            p5t_p41t = 1/turbine(m_41,N_t,"m","N","pi")
            eta_t = turbine(m_41,N_t,"m","N","eta")*0.9
    
            T5t_T41t = 1 + eta_t*((p5t_p41t)**((gamma_e-1)/gamma_e) - 1)
            m_5 = m_41/((1+f_assumed)*(1-b_3)+momentum_factor*b_3)*((1+f_assumed)*(1-b_3)+b_3)*np.sqrt(T5t_T41t)/p5t_p41t

            k1, k2 = 4.75, 1.475

            def cd(p9t_p9):

                CD_min, CD_max = 0.68, 0.974
                return (CD_max-CD_min)*(1-k2**(-k1*(p9t_p9-1))) + CD_min

            def m(p9t_p9):

                M9 = np.sqrt(2/(gamma_e-1)*(p9t_p9**((gamma_e-1)/gamma_e)-1))
                return cd(p9t_p9)*A8*(p_ref*1e5)/np.sqrt(R*T_ref)*np.sqrt(gamma_e)*M9* \
                (1+(gamma_e-1)/2*M9**2)**(-(gamma_e+1)/(2*(gamma_e-1)))

            # The maximum corrected mass flow according to A8 is calculated:

            m_5_max = (A8*p_ref*1e5)/np.sqrt(R*T_ref)*np.sqrt(gamma_e)*(2/(gamma_e+1))**((gamma_e+1)/(2*(gamma_e-1)))           

            if m_5 <= m_5_max:

                choked = False

                ## Nozzle equation:

                def nozzle_eq(p9t_p9):

                    return m(p9t_p9) - m_5

                # The nozzle is either sub-critical or critical:
                
                p9_p5t = 1/newton(nozzle_eq, np.linspace(1,1.4,100))
                p9_p5t = p9_p5t[0]
                M9 = np.sqrt(2/(gamma_e-1)*(((1/p9_p5t)**((gamma_e-1)/gamma_e)-1)))

                # Boundary condition (subsonic exit, Pascal's law):

                p9_p0 = 1

                # The intake is calculated according to the boundary condition:
                
                p2t_p0_est = p9_p0/(p9_p5t*p5t_p41t*p41t_p4t*p4t_p3t*p3t_p2t)
                signed_error = (p2t_p0_est - p2t_p0)/p2t_p0

                if iterate:
                    
                    return signed_error
                
                else:

                    # Exterior static temperature relation:

                    T9_T5t = 1/(1+(gamma_e-1)/2*M9**2)
                    T9_T0 = T9_T5t*T5t_T41t*T41t_T4t*T4t_T3t*T3t_T2t*T2t_T0

                    if nozzle == "conv":

                        A9_A8 = np.NaN

                    elif nozzle == "conv-div":

                        A9_A8 = 1
                    
                    # Characteristic parameters:

                    E = m_0*(((1+f_assumed)*(1-b_3)+b_3)*np.sqrt(gamma_e*R*T_ref)*M9*np.sqrt(T9_T0) - np.sqrt(gamma_c*R*T_ref)*M0) + \
                    (p9_p0 - 1)*(p_ref*1e5*A8)
                    Isp = E/m_0
                    TSFC = ((1+f_assumed)*(Cp_e/Cp_c)*T4t_T3t -1)*(Cp_c*T_ref)*(T3t_T2t*T2t_T0)*(1-b_3)/Isp

                    return m_0, T2t_T0, p2t_p0, m_2, T3t_T2t, p3t_p2t, eta_c, N_c, \
                    m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, \
                    choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, fuel_param, N              

            elif np.isnan(m_5):

                choked = True

                if iterate:
                    
                    return np.NaN
                
                else:

                    return np.NaN*np.empty(32)
                
            else:

                choked = True

                # The corrected mass flow in (5t) can not exceed the maximum value imposed by the throat:

                signed_error = (m_5 - m_5_max)/m_5_max

                if nozzle == "conv":

                    T9_T5t = 2/(gamma_e+1)
                    p9_p5t = (2/(gamma_e+1))**(gamma_e/(gamma_e-1))
                    M9 = 1

                    # Exterior static conditions relation:

                    T9_T0 = T9_T5t*T5t_T41t*T41t_T4t*T4t_T3t*T3t_T2t*T2t_T0
                    p9_p0 = p9_p5t*p5t_p41t*p41t_p4t*p4t_p3t*p3t_p2t*p2t_p0

                    if p9_p0 < 1:

                        if iterate:

                            return np.NaN
                        
                        else:

                            return np.NaN*np.empty(32)
                        
                    A9_A8 = np.NaN

                    # Characteristic parameters:

                    E = m_0*(((1+f_assumed)*(1-b_3)+b_3)*np.sqrt(gamma_e*R*T_ref)*M9*np.sqrt(T9_T0) - np.sqrt(gamma_c*R*T_ref)*M0) + \
                    (p9_p0 - 1)*(p_ref*1e5*A8*cd(1/p9_p5t))
                    Isp = E/m_0
                    TSFC = ((1+f_assumed)*(Cp_e/Cp_c)*T4t_T3t -1)*(Cp_c*T_ref)*(T3t_T2t*T2t_T0)*(1-b_3)/Isp

                elif nozzle == "conv-div":

                    # The nozzle area law is calculated to ensure a full expansion for all Mach numbers:

                    p9_p0 = 1

                    p9_p5t = p9_p0/(p5t_p41t*p41t_p4t*p4t_p3t*p3t_p2t*p2t_p0)
                    M9 = np.sqrt(2/(gamma_e-1)*((1/p9_p5t)**((gamma_e-1)/gamma_e))-1)
                    
                    # Exterior static temperature relation:
                    
                    if M9 < 1:

                        if iterate:

                            return np.NaN
                        
                        else:

                            return np.NaN*np.empty(32)
                        
                    T9_T5t = 1/(1+(gamma_e-1)/2*M9**2)
                    T9_T0 = T9_T5t*T5t_T41t*T41t_T4t*T4t_T3t*T3t_T2t*T2t_T0

                    # An extra nozzle efficiency between 8 and 9 is considered to be quadratic with M9 
                    # due to a gradual apperture of the nozzle petals. Favorable pressure gradient.
                    # No internal shockwaves due to underexpansion.

                    eta_n_exit = 1 - 0.01*(M9**2-1)
                    A9_A8 = (1+(1-1/eta_n_exit)*(gamma_e-1)/2*M9**2)**(-gamma_e/(gamma_e-1))*(1/M9)*(2/(gamma_e+1)* \
                    (1+(gamma_e-1)/2*M9**2))**((gamma_e+1)/(2*(gamma_e-1)))

                    A9 = A8*A9_A8

                    # Characteristic parameters:

                    E = m_0*(((1+f_assumed)*(1-b_3)+b_3)*np.sqrt(gamma_e*R*T_ref)*M9*np.sqrt(T9_T0) - np.sqrt(gamma_c*R*T_ref)*M0) + \
                    (p9_p0 - 1)*(p_ref*1e5*A9)
                    Isp = E/m_0
                    TSFC = ((1+f_assumed)*(Cp_e/Cp_c)*T4t_T3t -1)*(Cp_c*T_ref)*(T3t_T2t*T2t_T0)*(1-b_3)/Isp

                else:

                    print("Unknown type of nozzle. Options: 'conv' / 'conv-div'")
                    print(" ")
                    exit()

                if iterate:
                    
                    return signed_error
                
                else:

                    return m_0, T2t_T0, p2t_p0, m_2, T3t_T2t, p3t_p2t, eta_c, N_c, \
                    m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, \
                    choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, fuel_param, N

        beta_c = beta0    # Start from an initial value
        dBeta = 1e-30     # Increment for numerical derivatives

        error = np.NaN
        iterations = 0

        start = time.time()
        timeOut = 35      # [s]
        i = 0
        
        while np.isnan(error) or error > tolerance:

            if time.time() - start > timeOut:

                print("(TransEN) -> Limit time reached, no convergence. Modify function parameters")

                return np.NaN*np.empty(32)
            
            elif time.time() - start > timeOut/2 and error < 3e-2:

                print("(TransEN) -> Warning. Low rate of convergence. Assumed low error restarting point.")

                break
        
            signed_error = transientLoop(beta_c, True)
            error = abs(signed_error)

            if np.isnan(signed_error) and iterations < num_iter0:

                if iterations <= np.floor(num_iter0/2):
                    beta_c = beta0 + iterations/np.floor(num_iter0/2)*(1-beta0)
                elif iterations > np.floor(num_iter0/2):
                    beta_c = beta0 - (iterations-np.floor(num_iter0/2))/np.floor(num_iter0/2)*beta0

                iterations = iterations + 1

            elif np.isnan(signed_error) and iterations == num_iter0:

                return np.NaN*np.empty(32) 
                
            elif ~np.isnan(signed_error) and error > tolerance:

                # Functions of numerical derivatives:
                
                F = dBeta*(signed_error - transientLoop(beta_c - dBeta, True))/(signed_error - 2*transientLoop(beta_c - dBeta, True) + \
                transientLoop(beta_c - 2*dBeta, True))
                G = (signed_error*dBeta**2)/(signed_error - 2*transientLoop(beta_c - dBeta, True) + \
                transientLoop(beta_c - 2*dBeta, True))
            
                # Halley's method (cubic order of convergence):

                beta_c_star = beta_c + np.min([-F-np.sqrt(F**2-2*G),-F+np.sqrt(F**2-2*G)])
                beta_c = beta_c_star*(1-relaxation_factor) + beta_c*relaxation_factor

        return  transientLoop(beta_c, False)

    # Initial Value Problem Solver:
    
    N = N_0

    n = np.size(time_fuel_param)
    nm = 2                          # Number of measured variables

    iteration, beta0 = 0, 0.6
    time_N, current_time = [], []

    w = np.empty([n,32])
    dN = 1e-2

    P = 0
    Phi_s = 1000
    Rk = np.zeros([nm,nm])

    Rk[0,0] = 2e-4
    Rk[1,1] = 1e-4

    file = open("Custom/CSV Files/Validation/Measurements.csv", 'r')
    data = np.double(np.array(list(csv.reader(file))))

    while iteration in range(n):

        # Update fuel parameter and time:

        fuel_param = time_fuel_param[iteration]
        current_time.append(iteration*delta_t)

        # Differentially advanced solution:

        m_0, T2t_T0, p2t_p0, m_2, T3t_T2t, p3t_p2t, eta_c, N_c, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, \
        choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, fuel_param, N  = \
        transientOperationPoint(beta0, N + dN, fuel_param, num_iter0, relaxation_factor, tolerance)

        f_df = (m_0*Cp_c*T2t_T0*T_ref*((Cp_e/Cp_c)*eta_m*((1+f_assumed)*(1-b_3)+momentum_factor*b_3)*(1-T5t_T41t)*\
        T41t_T4t*T4t_T3t*T3t_T2t - (T3t_T2t-1)))/((2*np.pi/60)**2*I*N)

        h_dh = np.array([m_0, p3t_p2t*p2t_p0])

        # Change back to the desired speed:

        N = N - dN

        # Solution:

        m_0, T2t_T0, p2t_p0, m_2, T3t_T2t, p3t_p2t, eta_c, N_c, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, \
        choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, fuel_param, N  = \
        transientOperationPoint(beta0, N, fuel_param, num_iter0, relaxation_factor, tolerance)

        f = (m_0*Cp_c*T2t_T0*T_ref*((Cp_e/Cp_c)*eta_m*((1+f_assumed)*(1-b_3)+momentum_factor*b_3)*(1-T5t_T41t)*\
        T41t_T4t*T4t_T3t*T3t_T2t - (T3t_T2t-1)))/((2*np.pi/60)**2*I*N)

        h = np.array([m_0, p3t_p2t*p2t_p0])

        # Save the solution:

        w[iteration,:] = np.array([m_0, T2t_T0, p2t_p0, m_2, T3t_T2t, p3t_p2t, eta_c, N_c, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, \
        choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, fuel_param, N])

        # Kalman Filter Parameters:
        
        H = (h_dh - h)/dN     # Measurement Jacobian
        F = (f_df - f)/dN     # State Jacobian 
        Phi = 1 + F*delta_t   # State Propagation Matrix

        Q = Phi_s*((F**2)*delta_t**3/3  + F*delta_t**2 + delta_t)   # Process Noise Covariance Matrix
        M = Phi*P*Phi + Q                                           # Non-updated Filter Covariance Matrix

        K = M*np.dot(H,np.linalg.inv(M*np.outer(H,H) + Rk))         # Kalman Gain
        P = (1 - np.dot(K,H))*M                                     # Updated covariance matrix

        # Coupling Excess Power:

        if propagator == 'Euler':

            N_hat = N + delta_t*f
            time_N.append(N)

        m_0, T2t_T0, p2t_p0, m_2, T3t_T2t, p3t_p2t, eta_c, N_c, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, \
        choked, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, fuel_param, N  = \
        transientOperationPoint(beta0, N_hat, fuel_param, num_iter0, relaxation_factor, tolerance)

        # Kalman Filter implementation:

        h = np.array([m_0, p3t_p2t*p2t_p0])

        z = [data[iteration,0],data[iteration,2]]
        N = N_hat + np.dot(K,z - h)

        # Update the iteration counter:

        iteration = iteration + 1
        beta0 = compressor(p3t_p2t,N_c,"pi","N","beta")    # Initialize according to last point

        print("Iteration: "+str(iteration)+"/"+str(n)+": "+str(np.round(iteration/n*100,2))+"%")

    return w, time_N, current_time
