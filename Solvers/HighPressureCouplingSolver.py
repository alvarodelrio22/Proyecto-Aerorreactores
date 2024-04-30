from Components.ComponentMap import compressor, turbine
from Components.DesignVariables import p_ref, T_ref, gamma_c, gamma_e, Cp_c, Cp_e, N_ref_HPC, N_ref_HPT, b_25, b_3, eta_mHP as eta_m
import numpy as np

def hpCoupling(beta_HPC,N_HPC,max_iterations,relaxation_factor,representing):

## LPC Outlet - HPC Inlet (25t) -----------------------------------------------------------------------------------------------------------------------------
    
    m_25 = compressor(beta_HPC,N_HPC,"beta","N","m","HPC")

## HPC Outlet - Combustion Chamber Inlet (3t) ---------------------------------------------------------------------------------------------------------------
  
    p3t_p25t = compressor(beta_HPC,N_HPC,"beta","N","pi","HPC")
    eta_HPC = compressor(beta_HPC,N_HPC,"beta","N","eta","HPC")

    if np.isnan(m_25) or np.isnan(p3t_p25t) or np.isnan(eta_HPC):
        
        if representing == True:
            return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
        else:
            return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
            
    T3t_T25t = 1 + 1/eta_HPC*(p3t_p25t**((gamma_c-1)/gamma_c)-1)
    m_3 = m_25/(1-b_25)*(1-b_25-b_3)*np.sqrt(T3t_T25t)/p3t_p25t
  
## Combustion Chamber Outlet - HPT Inlet (4t) ---------------------------------------------------------------------------------------------------------------

    T4t_T25t = 3.5
    num_iter0 = 50

    error = []
    iterations = 0
    f_assumed = 0.025

    while np.size(error) == 0 or error >= 1e-4:

        iterations = iterations + 1
      
        T4t_T3t = T4t_T25t/T3t_T25t

        PLF = 0.27 + 0.015*(T4t_T25t/T3t_T25t - 1)
        p4t_p3t = 1 - PLF/2*(m_3*np.sqrt(287)*np.sqrt(T_ref)/(p_ref*1e5)/0.06)**2
  
    # A FAR value is assumed. Gross temperature magnitude is requiered to
    # obtain its actual form. The problem is almost independent of the
    # value of FAR and bleeds, making a reduction in DOF possible.
  
        m_4 =  m_3/(1-b_25-b_3)*(1+f_assumed-b_25-b_3)*np.sqrt(T4t_T3t)/p4t_p3t

## NGV Bleed Injection (41t) -------------------------------------------------------------------------------------------------------------------------------

        p41t_p4t = 1
        T41t_T4t = (m_4 + b_3*m_25/(1-b_25)*np.sqrt(T4t_T25t)*1/p3t_p25t*1/p4t_p3t*
              Cp_c/Cp_e*(1/T4t_T25t))/(m_4 + b_3*m_25*np.sqrt(T4t_T25t)*1/p3t_p25t*1/p4t_p3t)
   
        m_41 = m_4/(1+f_assumed-b_25-b_3)*(1+f_assumed-b_25)*np.sqrt(T41t_T4t)/p41t_p4t

## HPT Outlet - LPT Inlet (45t) ----------------------------------------------------------------------------------------------------------------------------

        N_HPT = N_HPC*N_ref_HPC/np.sqrt(T4t_T25t)/N_ref_HPT
          
        p45t_p41t = 1/turbine(m_41,N_HPT,"m","N","pi","HPT")
        eta_HPT = turbine(m_41,N_HPT,"m","N","eta","HPT")
   
        T45t_T41t = 1 + eta_HPT*((p45t_p41t)**((gamma_e-1)/gamma_e) - 1)
        m_45 = m_41*np.sqrt(T45t_T41t)/p45t_p41t
   
## Power balance (HPC-HPT) and solution
# Determination of a new value for T4t_T25t. The bleed is considered to be injected here.
# Momentum
      
        T4t_T25t_last = T4t_T25t

        T4t_T25t = Cp_c*(T3t_T25t - 1)*(1 - b_25)/(eta_m*Cp_e*(1 - T45t_T41t)*
        T41t_T4t*(1 + f_assumed - 0.5*b_3 - b_25))*(1 - relaxation_factor) + T4t_T25t*relaxation_factor
   
## Try new initial values of  T4t_T25t that might converge initially
   
        if np.isnan(T4t_T25t) and iterations <= num_iter0:

            T4t_T25t = T3t_T25t*(num_iter0 - iterations)/num_iter0 + 7*iterations/num_iter0
     
        elif iterations > max_iterations:
        
            if representing == True:
                return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
            else:
                return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
     
        else:
    
            error = abs(T4t_T25t - T4t_T25t_last)

    T45t_T25t = T45t_T41t*T41t_T4t*T4t_T25t
    p45t_p25t = p45t_p41t*p41t_p4t*p4t_p3t*p3t_p25t
    fuel_param = (1+f_assumed-b_25-b_3)*Cp_e/Cp_c*T4t_T3t - 1
    
    if representing == True:
        
        if np.isnan(T4t_T25t):
            return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN
        else:
            return m_25, p3t_p25t, eta_HPC, m_41, p45t_p41t, eta_HPT
    
    else:

        return m_25, p3t_p25t, eta_HPC, m_41, p45t_p41t, eta_HPT, m_45, T45t_T25t, p45t_p25t, fuel_param
        