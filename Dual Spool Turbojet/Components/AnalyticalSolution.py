from Components.DesignVariables import p_ref, T_ref, gamma_c, gamma_e, Cp_c, Cp_e, b_25, \
f_assumed, b_3, momentum_factor, R, A8, eta_d, eta_n, eta_LPC, eta_HPC, eta_HPT, eta_LPT, k_LPT, k_HPT, k_n, pi_CC, pi_NGV, \
pi_HPT, pi_LPT, tau_NGV, tau_HPT, tau_LPT, eta_mHP, eta_mLP
import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt

def analyticalSolution(M0, pi_LPC, nozzle):

    T2t_T0 = 1 + (gamma_c - 1)/2*M0**2
    p2t_p0 = (1 + (gamma_c - 1)/2*eta_d*M0**2)**(gamma_c/(gamma_c - 1))

    T25t_T2t = 1 + 1/eta_LPC*(pi_LPC**((gamma_c-1)/gamma_c)-1)
    p25t_p2t = pi_LPC

    p4t_p3t = pi_CC

    T41t_T4t = tau_NGV
    p41t_p4t = pi_NGV

    T45t_T41t = tau_HPT
    p45t_p41t = 1/pi_HPT

    def p3tp25t(p45t_p5t):
    
        return (1+(eta_mHP*eta_HPC)/(eta_mLP*eta_LPC)* \
        ((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)/((1-b_25)*((1+f_assumed)*(1-b_25-b_3)+b_3)) * \
        (1-tau_HPT)/(tau_HPT*eta_LPT*(1-(1/p45t_p5t)**((gamma_e-1)/gamma_e)))* \
        (pi_LPC**((gamma_c-1)/gamma_c)-1)/(1+(1/eta_LPC)*(pi_LPC**((gamma_c-1)/gamma_c)-1)))**(gamma_c/(gamma_c-1))

    # Calculation of unlocking conditions:

    pi_LPC_unlock = 1/(p2t_p0*p3tp25t(pi_LPT)*pi_CC*pi_NGV*(1/pi_HPT)*(1/pi_LPT)*(1-(1/eta_n)*(gamma_e-1)/(gamma_e+1))**(gamma_e/(gamma_e-1)))

    # Only when the pressure equalizes:

    def p9p5t(p45t_p5t):

        return 1/(p2t_p0*pi_LPC*p3tp25t(p45t_p5t)*pi_CC*pi_NGV*(1/pi_HPT)*(1/p45t_p5t))

    if pi_LPC <= pi_LPC_unlock:
        
        def f(p45t_p5t):

            return k_LPT*p45t_p5t*np.sqrt(1-eta_LPT*(1-(1/p45t_p5t)**((gamma_e-1)/gamma_e))) - \
            A8*(p_ref*10**5)/np.sqrt(R*T_ref)*p9p5t(p45t_p5t)*(1-eta_n*(1-p9p5t(p45t_p5t)**((gamma_e-1)/gamma_e)))**(-1) * \
            np.sqrt((2*gamma_e*eta_n)/(gamma_e-1)*(1-p9p5t(p45t_p5t)**((gamma_e-1)/gamma_e)))
        
        p5t_p45t = 1/newton(f,np.linspace(1,3,10000))
        p5t_p45t = p5t_p45t[~np.isnan(p5t_p45t)][-1]
        T5t_T45t = 1-eta_LPT*(1-(p5t_p45t)**((gamma_e-1)/gamma_e))

        p9_p5t = p9p5t(1/p5t_p45t)
        M9 = np.sqrt(2/(gamma_e-1)*(1/(1-eta_n*(1-(p9_p5t)**((gamma_e-1)/gamma_e)))-1))
        T9_T5t = 1/(1+(gamma_e-1)/2*M9**2)     

        choke = False

        A9 = A8
        A9_A8 = 1

    elif pi_LPC > pi_LPC_unlock:

        p5t_p45t = 1/pi_LPT
        T5t_T45t = tau_LPT
        choke = True

        if nozzle == 'conv':
        
            T9_T5t = 2/(gamma_e+1)
            p9_p5t = (1-(1/eta_n)*(gamma_e-1)/(gamma_e+1))**(gamma_e/(gamma_e-1))
            M9 = 1

            A9 = A8
            A9_A8 = 1
        
        elif nozzle == 'conv-div':

            p9_p5t = p9p5t(1/p5t_p45t)
            M9 = np.sqrt(2/(gamma_e-1)*(1/(1-eta_n*(1-(p9_p5t)**((gamma_e-1)/gamma_e)))-1))
            T9_T5t = 1/(1+(gamma_e-1)/2*M9**2)

            eta_n_exit = 1 - 0.01*(M9**2-1)
            A9_A8 = (1+(1-1/eta_n_exit)*(gamma_e-1)/2*M9**2)**(-gamma_e/(gamma_e-1))*(1/M9)*(2/(gamma_e+1)* \
            (1+(gamma_e-1)/2*M9**2))**((gamma_e+1)/(2*(gamma_e-1)))

            A9 = A8*A9_A8

    # Now corrected mass flows are determined:

    m_2 = k_HPT*pi_NGV*pi_CC*np.sqrt((Cp_e/Cp_c)*eta_LPC*eta_mLP*tau_HPT*(1-T5t_T45t)*((1+f_assumed)*(1-b_25-b_3)+b_3))*\
    1/((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)*pi_LPC/np.sqrt(pi_LPC**((gamma_c-1)/gamma_c)-1)*\
    (1+(eta_mHP*eta_HPC)/(eta_mLP*eta_LPC)* \
    ((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)/((1-b_25)*((1+f_assumed)*(1-b_25-b_3)+b_3)) * \
    (1-tau_HPT)/(tau_HPT*eta_LPT*(1-p5t_p45t**((gamma_e-1)/gamma_e)))* \
    (pi_LPC**((gamma_c-1)/gamma_c)-1)/(1+(1/eta_LPC)*(pi_LPC**((gamma_c-1)/gamma_c)-1)))**(gamma_c/(gamma_c-1))

    T41t_T25t = (T25t_T2t - 1)/(eta_mLP*(Cp_e/Cp_c)*((1+f_assumed)*(1-b_25-b_3)+b_3)*tau_HPT*(1-T5t_T45t))/(T25t_T2t)
    
    # The inlet corrected mass flow is determined:

    m_0 = m_2/np.sqrt(T2t_T0)*p2t_p0

    # The HPC:

    m_25 = m_2*(1-b_25)*np.sqrt(T25t_T2t)/p25t_p2t
    p3t_p25t = p3tp25t(1/p5t_p45t)
    T3t_T25t = 1 + (1/eta_HPC)*(p3t_p25t**((gamma_c-1)/gamma_c)-1)
    
    # The combustion chamber:

    m_3 = m_25/(1-b_25)*(1-b_25-b_3)*np.sqrt(T3t_T25t)/p3t_p25t
    T4t_T3t = 1/T41t_T4t*T41t_T25t/T3t_T25t
    fuel_parameter = (T3t_T25t*T25t_T2t*T2t_T0)*(Cp_e/Cp_c*(1+f_assumed)*T4t_T3t - 1)

    # The NGV Cooling bleed:

    m_4 =  m_3*(1+f_assumed)*np.sqrt(T4t_T3t)/p4t_p3t
    m_41 = m_4/((1+f_assumed)*(1-b_25-b_3))*((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)*np.sqrt(T41t_T4t)/p41t_p4t

    # The HPT:

    m_45 = m_41/((1+f_assumed)*(1-b_25-b_3)+momentum_factor*b_3)*((1+f_assumed)*(1-b_25-b_3)+b_3)*np.sqrt(T45t_T41t)/p45t_p41t

    # The last thermodynamic variables:

    m_5 = m_45*np.sqrt(T5t_T45t)/p5t_p45t

    p9_p0 = p9_p5t*p5t_p45t*p45t_p41t*p41t_p4t*p41t_p4t*p4t_p3t*p3t_p25t*p25t_p2t*p2t_p0
    T9_T0 = T9_T5t*T5t_T45t*T45t_T41t*T41t_T4t*T41t_T4t*T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0

    if nozzle == 'conv':

        E = m_0*(((1+f_assumed)*(1-b_25-b_3)+b_3)*np.sqrt(gamma_e*R*T_ref)*M9*np.sqrt(T9_T0) - np.sqrt(gamma_c*R*T_ref)*M0) + \
        (p9_p0 - 1)*(p_ref*1e5*A8)

    elif nozzle == 'conv-div':
        
        E = m_0*(((1+f_assumed)*(1-b_25-b_3)+b_3)*np.sqrt(gamma_e*R*T_ref)*M9*np.sqrt(T9_T0) - np.sqrt(gamma_c*R*T_ref)*M0) + \
        (p9_p0 - 1)*(p_ref*1e5*A9)

    Isp = E/m_0
    TSFC = ((1+f_assumed)*(Cp_e/Cp_c)*T4t_T3t -1)*(Cp_c*T_ref)*(T3t_T25t*T25t_T2t*T2t_T0)*(1-b_25-b_3)/Isp

    Isp = E/m_0
    TSFC = ((1+f_assumed)*(Cp_e/Cp_c)*T4t_T3t -1)*(Cp_c*T_ref)*(T3t_T25t*T25t_T2t*T2t_T0)*(1-b_25-b_3)/Isp

    return T2t_T0, p2t_p0, m_0, T25t_T2t, p25t_p2t, m_2, T3t_T25t, p3t_p25t, m_25, T4t_T3t, p4t_p3t, m_3, T41t_T4t, p41t_p4t, m_4, \
    T45t_T41t, p45t_p41t, m_41, T5t_T45t, p5t_p45t, m_45, m_5, choke, T9_T5t, p9_p5t, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, fuel_parameter


