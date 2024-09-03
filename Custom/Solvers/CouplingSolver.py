from Components.ComponentMap import compressor, turbine
from Components.DesignVariables import p_ref, T_ref, R, gamma_c, gamma_e, \
Cp_c, Cp_e, N_ref_c, N_ref_t, b_3, eta_m, f_assumed, momentum_factor

import time, numpy as np

## COUPLING SOLVER :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Iterative method (Succesive Over - Relaxation)
# Arguments: Compressor beta line and Compressor Relative Corrected Spool Speed

## Function definition: ------------------------------------------------------------------------------------------------------------------------------------

def coupling(beta_c,N_c,num_iter0,relaxation_factor,representing):

    # Bleed values are all referred to the intake air.

    if np.isnan(beta_c) or np.isnan(N_c):

        if representing:

            return np.NaN*np.empty(6) 
        
        else:

            return np.NaN*np.empty(18) 
        
## Diffuser Outlet - Compressor Inlet (2t) ------------------------------------------------------------------------------------------------------------------
    
    m_2 = compressor(beta_c,N_c,"beta","N","m")

## Compressor Outlet - Combustion Chamber Inlet (3t) --------------------------------------------------------------------------------------------------------
  
    p3t_p2t = compressor(beta_c,N_c,"beta","N","pi")
    eta_c = compressor(beta_c,N_c,"beta","N","eta")

    if np.isnan(m_2) or np.isnan(p3t_p2t) or np.isnan(eta_c):
        
        if representing:

            return np.NaN*np.empty(6) 
        
        else:

            return np.NaN*np.empty(18) 
            
    T3t_T2t = 1 + 1/eta_c*(p3t_p2t**((gamma_c-1)/gamma_c)-1)
    m_3 = m_2*(1-b_3)*np.sqrt(T3t_T2t)/p3t_p2t
  
## Combustion Chamber Outlet - NGV Bleed Injection (4t) -----------------------------------------------------------------------------------------------------

    T4t_T2t = 3.5       # Start from an initial value
    T4t_T2t_max = 5     # Maximum initial guess

    error = np.NaN
    tolerance = 1e-6
    iterations = 0

    start = time.time()
    timeOut = 1.0        # [s]

    while np.isnan(error) or error >= tolerance:

        if time.time() - start > timeOut:
        
            if representing:

                print("Limit time reached, no convergence. Modify function parameters")

                return np.NaN*np.empty(6) 
            
            else:

                return np.NaN*np.empty(18) 
      
        T4t_T3t = T4t_T2t/T3t_T2t

        A3 = 0.0025  # [m^2]
        PLF = 1/3 + 5/27*(T4t_T3t - 1)
        p4t_p3t = 1 - PLF/2*(m_2*(1-b_3)*np.sqrt(R)*np.sqrt(T_ref)/(p_ref*1e5)/A3)**2

    # A FAR value is assumed, with respect to the combustion chamber intake air. 
    # Gross temperature magnitude is requiered to obtain its actual form. 
    # The problem is almost independent of the value of FAR and bleeds, 
    # making a reduction in DOFs possible.
  
        m_4 =  m_3*(1+f_assumed)*np.sqrt(T4t_T3t)/p4t_p3t

## NGV Bleed Injection - Turbine Inlet (41t) ----------------------------------------------------------------------------------------------------------------

        p41t_p4t = 1
        T41t_T4t = ((1+f_assumed)*(1-b_3) + b_3*(Cp_c/Cp_e)*(1/T4t_T3t))/((1+f_assumed)*(1-b_3) + b_3)
   
        m_41 = m_4/((1+f_assumed)*(1-b_3))*((1+f_assumed)*(1-b_3)+momentum_factor*b_3)*np.sqrt(T41t_T4t)/p41t_p4t

## Turbine Outlet - Nozzle Inlet (5t) -----------------------------------------------------------------------------------------------------------------------

        N_t = N_c*N_ref_c/np.sqrt(T41t_T4t*T4t_T2t)/N_ref_t
          
        p5t_p41t = 1/turbine(m_41,N_t,"m","N","pi")
        eta_t = turbine(m_41,N_t,"m","N","eta")*0.9
   
        T5t_T41t = 1 + eta_t*((p5t_p41t)**((gamma_e-1)/gamma_e) - 1)
        m_5 = m_41/((1+f_assumed)*(1-b_3)+momentum_factor*b_3)*((1+f_assumed)*(1-b_3)+b_3)*np.sqrt(T5t_T41t)/p5t_p41t
   
## Power balance (C-T) and solution
# Determination of a new value for T4t_T2t. The bleed is considered to be injected here.
# Momentum
      
        T4t_T2t_last = T4t_T2t

        T4t_T2t = 1/T41t_T4t*(T3t_T2t - 1)/(eta_m*(Cp_e/Cp_c)*(1 - T5t_T41t)*
        ((1+f_assumed)*(1-b_3)+momentum_factor*b_3))*(1 - relaxation_factor) + T4t_T2t*relaxation_factor

        error = abs(T4t_T2t - T4t_T2t_last)
   
## Try new initial values of  T4t_T2t that might converge initially
   
        if np.isnan(error) and iterations < num_iter0:

            T4t_T2t = T3t_T2t*(num_iter0 - iterations)/num_iter0 + T4t_T2t_max*iterations/num_iter0
            iterations = iterations + 1

        elif np.isnan(error) and iterations == num_iter0:

            if representing:

                return np.NaN*np.empty(6) 
            
            else:

                return np.NaN*np.empty(18) 

    load_param = T4t_T3t - 1
    fuel_param_uncorrected = (1+f_assumed)*Cp_e/Cp_c*(load_param + 1) - 1

    if representing:
        
        return m_2, p3t_p2t, eta_c, m_41, p5t_p41t, eta_t
    
    else:

        return m_2, T3t_T2t, p3t_p2t, eta_c, m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, \
        m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, load_param, fuel_param_uncorrected
        