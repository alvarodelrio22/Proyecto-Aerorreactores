from Components.ComponentMap import compressor, turbine
from Miscellaneous.AuxilliaryFunctions import relaxationFactor
from Components.DesignVariables import p_ref, T_ref, R, gamma_c, gamma_e, \
Cp_c, Cp_e, N_ref_LPC, N_ref_HPC, N_ref_HPT, N_ref_LPT, b_25, b_3, eta_mHP, eta_mLP, f_assumed, momentum_factor

import time, numpy as np

## ENGINE INTERNAL COUPLING ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## HIGH PRESSURE COUPLING SOLVER ---------------------------------------------------------------------------------------------------------------------------
# Iterative method (Succesive Over - Relaxation)
# Arguments: LPC beta line and LPC Relative Corrected Spool Speed

## Function definition: ------------------------------------------------------------------------------------------------------------------------------------

def hpCoupling(beta_HPC,N_HPC,num_iter0,relaxation_factor,representing):

    # Bleed values are all referred to the intake air.

    if np.isnan(beta_HPC) or np.isnan(N_HPC):

        if representing:

            return np.NaN*np.empty(6) 
        
        else:

            return np.NaN*np.empty(18) 
        
## LPC Outlet - HPC Inlet (25t) --------------------------------------------------------------------------------------------------------------------------
    
    m_25 = compressor(beta_HPC,N_HPC,"beta","N","m","HPC")

## HPC Outlet - Combustion Chamber Inlet (3t) ------------------------------------------------------------------------------------------------------------
  
    p3t_p25t = compressor(beta_HPC,N_HPC,"beta","N","pi","HPC")
    eta_HPC = compressor(beta_HPC,N_HPC,"beta","N","eta","HPC")

    if np.isnan(m_25) or np.isnan(p3t_p25t) or np.isnan(eta_HPC):
        
        if representing:

            return np.NaN*np.empty(6) 
        
        else:

            return np.NaN*np.empty(18) 
            
    T3t_T25t = 1 + 1/eta_HPC*(p3t_p25t**((gamma_c-1)/gamma_c)-1)
    m_3 = m_25/(1-b_25)*(1-b_25-b_3)*np.sqrt(T3t_T25t)/p3t_p25t
  
## Combustion Chamber Outlet - NGV Bleed Injection (4t) -----------------------------------------------------------------------------------------------------

    T4t_T25t = 3.5       # Start from an initial value
    T4t_T25t_max = 5     # Maximum initial guess

    error = np.NaN
    tolerance = 1e-6
    iterations = 0

    start = time.time()
    timeOut = 0.5        # [s]

    while np.isnan(error) or error >= tolerance:

        if time.time() - start > timeOut:
        
            if representing:

                print("(HP) -> Limit time reached, no convergence. Modify function parameters")

                return np.NaN*np.empty(6) 
            
            else:

                return np.NaN*np.empty(18) 
      
        T4t_T3t = T4t_T25t/T3t_T25t

        A3 = 0.9  # [m^2]
        PLF = 27 + 15*(T4t_T3t - 1)
        p4t_p3t = 1 - PLF/2*(m_3*np.sqrt(R)*np.sqrt(T_ref)/(p_ref*1e5)/A3)**2
  
    # A FAR value is assumed, with respect to the combustion chamber intake air. 
    # Gross temperature magnitude is requiered to obtain its actual form. 
    # The problem is almost independent of the value of FAR and bleeds, 
    # making a reduction in DOFs possible.
  
        m_4 =  m_3*(1+f_assumed)*np.sqrt(T4t_T3t)/p4t_p3t

## Combustion Chamber Outlet - NGV Bleed Injection (4t) -----------------------------------------------------------------------------------------------------

        p41t_p4t = 1
        T41t_T4t = ((1+f_assumed)*(1-b_25-b_3) + b_3*(Cp_c/Cp_e)*(1/T4t_T3t))/((1+f_assumed)*(1-b_25-b_3) + b_3)
   
        m_41 = m_4/((1+f_assumed)*(1-b_25-b_3))*((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)*np.sqrt(T41t_T4t)/p41t_p4t

## HPT Outlet - LPT Inlet (45t) -----------------------------------------------------------------------------------------------------------------------------

        N_HPT = (N_HPC*N_ref_HPC)/np.sqrt(T41t_T4t*T4t_T25t)/N_ref_HPT
          
        p45t_p41t = 1/turbine(m_41,N_HPT,"m","N","pi","HPT")
        eta_HPT = turbine(m_41,N_HPT,"m","N","eta","HPT")
   
        T45t_T41t = 1 + eta_HPT*((p45t_p41t)**((gamma_e-1)/gamma_e) - 1)
        m_45 = m_41/((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)*((1+f_assumed)*(1-b_25-b_3)+b_3)*np.sqrt(T45t_T41t)/p45t_p41t
   
## Power balance (HPC-HPT) and solution
# Determination of a new value for T4t_T25t. The bleed is considered to be injected here.
# Momentum
      
        T4t_T25t_last = T4t_T25t

        T4t_T25t = 1/T41t_T4t*(T3t_T25t - 1)*(1 - b_25)/(eta_mHP*(Cp_e/Cp_c)*(1 - T45t_T41t)*
        ((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3))*(1 - relaxation_factor) + T4t_T25t*relaxation_factor

        error = abs(T4t_T25t - T4t_T25t_last)
   
## Try new initial values of  T4t_T25t that might converge initially
   
        if np.isnan(error) and iterations < num_iter0:

            T4t_T25t = T3t_T25t*(num_iter0 - iterations)/num_iter0 + T4t_T25t_max*iterations/num_iter0
            iterations = iterations + 1

        elif np.isnan(error) and iterations == num_iter0:

            if representing:

                return np.NaN*np.empty(6) 
            
            else:

                return np.NaN*np.empty(18) 

    load_param = T4t_T3t - 1
    fuel_param_uncorrected = (1+f_assumed)*Cp_e/Cp_c*(load_param + 1) - 1
    
    if representing:
        
        return m_25, p3t_p25t, eta_HPC, m_41, p45t_p41t, eta_HPT
    
    else:

        return m_25, T3t_T25t, p3t_p25t, eta_HPC, m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, \
        m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, load_param, fuel_param_uncorrected

    
## LOW PRESSURE COUPLING SOLVER ----------------------------------------------------------------------------------------------------------------------------
# Iterative method (Newton - Rhapson with Succesive Over - Relaxation)
# Arguments: LPC beta line and LPC Relative Corrected Spool Speed

# Relaxation factor map for the HP coupling solver: --------------------------------------------------------------------------------------------------------

num_iter0_HPC = 100

relaxation_B_lim = [1, 0.90, 0.75, 0]
relaxation_N_lim = [0.45, 0.75, 0.95, 1.08]
relaxation_matrix = \
[[0.65,0.75,0.85],  \
[0.65,0.75,0.90],   \
[0.70,0.85,0.90]]

## Function definition: ------------------------------------------------------------------------------------------------------------------------------------

def lpCoupling(beta_LPC,N_LPC,num_iter0,relaxation_factor,representing):

## Diffuser Outlet - LPC Inlet (2t) ------------------------------------------------------------------------------------------------------------------------

    m_2 = compressor(beta_LPC,N_LPC,"beta","N","m","LPC")

## LPC Outlet - HPC Inlet (25t) ----------------------------------------------------------------------------------------------------------------------------

    p25t_p2t = compressor(beta_LPC,N_LPC,"beta","N","pi","LPC")
    eta_LPC = compressor(beta_LPC,N_LPC,"beta","N","eta","LPC") 

    if np.isnan(m_2) or np.isnan(p25t_p2t) or np.isnan(eta_LPC):
        
        if representing:
        
            return np.NaN*np.empty(12) 
                
        else:

            return np.NaN*np.empty(28) 
           
    T25t_T2t = 1 + 1/eta_LPC*(p25t_p2t**((gamma_c-1)/gamma_c) - 1)
    m_HPC = m_2*(1-b_25)*np.sqrt(T25t_T2t)/p25t_p2t

## Solution to the high pressure coupling problem - LPT Inlet (45t) ----------------------------------------------------------------------------------------
# An iterative process is requiered (using beta achieves better ortogonality than N):

    def hpLoop(beta_HPC, iterate):

        N_HPC = compressor(m_HPC, beta_HPC, "m", "beta", "N", "HPC")

        if np.isnan(N_HPC):

            if iterate:

                return np.NaN
        
            else:

                 return np.NaN*np.empty(23)
        
        else:

            m_25, T3t_T25t, p3t_p25t, eta_HPC, m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, \
            m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, load_param, fuel_param_uncorrected = hpCoupling(beta_HPC,N_HPC, \
            num_iter0_HPC,relaxationFactor(beta_HPC,N_HPC,relaxation_B_lim,relaxation_N_lim,relaxation_matrix), False)

## LPT Outlet - Nozzle Inlet (5t) ------------------------------------------------------------------------------------------------------------------------

            N_LPT = (N_LPC*N_ref_LPC)/np.sqrt(T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t)/N_ref_LPT

            p5t_p45t = 1/turbine(m_45,N_LPT,"m","N","pi","LPT")
            eta_LPT = turbine(m_45,N_LPT,"m","N","eta","LPT")
            T5t_T45t = 1 + eta_LPT*(p5t_p45t**((gamma_e-1)/gamma_e) - 1)

            m_5 = m_45*np.sqrt(T5t_T45t)/p5t_p45t

            if iterate:

                # Power balance (HPC-HPT) and solution
                # Determination of a new value for beta_HPC

                T25t_T2t_est = 1/(1 - eta_mLP*(Cp_e/Cp_c)*T45t_T41t*T41t_T4t*T4t_T3t*T3t_T25t*((1+f_assumed)*(1-b_25-b_3)+b_3)*(1 - T5t_T45t))
                signed_error = T25t_T2t_est - T25t_T2t

                return signed_error
            
            else:

                return m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, \
                m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, load_param, fuel_param_uncorrected  
    
    # Number of initial guesses is determined by the number of iterations:
    
    beta_HPC = 0         # Start from an initial value (critical, for external coupling solution recommended: 0)
    dBeta = 1e-7         # Increment for numerical derivative

    error = np.NaN
    tolerance = 1e-5
    iterations = 0

    start = time.time()
    timeOut = 1.5        # [s]

    while np.isnan(error) or error >= tolerance:

        if time.time() - start > timeOut:
            
            if representing:

                print("(LP) -> Limit time reached, no convergence. Modify function parameters")

                return np.NaN*np.empty(12) 
            
            else:

                return np.NaN*np.empty(28) 
            
        signed_error = hpLoop(beta_HPC, True)
        error = abs(signed_error)

        if np.isnan(signed_error) and iterations < num_iter0:

            beta_HPC = iterations/(num_iter0-1)                 # Initialization is critical due to the
            iterations = iterations + 1                         # potential existence of multiple solutions

        elif np.isnan(signed_error) and iterations == num_iter0:

            if representing:

                return np.NaN*np.empty(12) 
        
            else:

                return np.NaN*np.empty(28) 
            
        elif ~np.isnan(signed_error) and error > tolerance:

            # Numerical derivative:

            dtau_dbeta =  (hpLoop(beta_HPC + dBeta, True) - signed_error)/dBeta

            # Newton - Rhapson method (second order) implementation:

            beta_HPC_star = beta_HPC - 1/dtau_dbeta*signed_error
            beta_HPC = beta_HPC_star*(1-relaxation_factor) + beta_HPC*relaxation_factor

    m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
    m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, load_param, fuel_param_uncorrected = hpLoop(beta_HPC, False)
    
    if representing:

        return m_2, p25t_p2t, eta_LPC, m_25, p3t_p25t, eta_HPC, m_41, \
        p45t_p41t, eta_HPT, m_45, p5t_p45t, eta_LPT
        
    else:

        return m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, \
        m_45, T5t_T45t, p5t_p45t, eta_LPT, N_LPT, m_5, load_param, fuel_param_uncorrected 
        