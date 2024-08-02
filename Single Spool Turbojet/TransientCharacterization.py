from Solvers.OperatingPointSolver import engOperation
from Solvers.TransientOperationSolver import transientOperation
from Miscellaneous.AuxilliaryFunctions import componentPlot
import numpy as np, matplotlib.pyplot as plt
from Components.DesignVariables import fuel_param_design

## TRANSIENT CHARACTERIZATION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Iterative method (Halley Method) + Initial Value Problem Numerical Solver (with choice of propagator)
# Arguments: Initial State Conditions and Corrected Fuel Parameter Time Evolution

## Simulation parameters ------------------------------------------------------------------------------------------------------------------------------------

M0 = 0
N_c = 0.65

simulation_time = 30  # [s]
num_iter0 = 500
n = 200

I = 23.43500/25  # [kg·m^2]

nozzle = 'conv'

# Initial state definition: ---------------------------------------------------------------------------------------------------------------------------------

m_0, T2t_T0, p2t_p0, eta_d, m_2, T3t_T2t, p3t_p2t, eta_c, \
m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, \
choked, T9_T5t, p9_p5t, eta_n, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, load_param, fuel_param_0, N_0   = \
engOperation(M0,N_c,nozzle,15,0.15)

if np.isnan(N_0):

    print("Error: Non-existing starting point.")
    exit()

w_0 = M0, I, N_0

# Fuel parameter function definition: -----------------------------------------------------------------------------------------------------------------------

time_fuel_param = np.empty(n)
tolerance = np.empty(n)
relaxation_factor = np.empty(n)

for i in range(n):

    if i <= np.floor(n/20):

        time_fuel_param[i] = (fuel_param_design*i + fuel_param_0*(np.floor(n/20)-i))/np.floor(n/20)

        tolerance[i] = 5e-4
        relaxation_factor[i] = 0.95

    elif np.floor(n/20) < i < np.floor(n/2) :

        time_fuel_param[i] = fuel_param_design

        tolerance[i] = 5e-4
        relaxation_factor[i] = 0

    elif np.floor(n/2) <= i <= np.floor(11*n/20):

        time_fuel_param[i] = (fuel_param_0*(i-np.floor(n/2)) + \
        fuel_param_design*((np.floor(11*n/20)-np.floor(n/2))-(i-np.floor(n/2))))/(np.floor(11*n/20)-np.floor(n/2))

        tolerance[i] = 5e-4
        relaxation_factor[i] = 0.5

    else:

        time_fuel_param[i] = fuel_param_0

        tolerance[i] = 5e-4

        if np.floor(14*n/20) < i <= np.floor(17*n/20):

            relaxation_factor[i] = 0.95

        else:

            relaxation_factor[i] = 0.95

## Run Simulation ------------------------------------------------------------------------------------------------------------------------------------------

w, time_N, current_time = transientOperation(w_0, time_fuel_param, 'Euler', simulation_time/n, relaxation_factor, num_iter0, tolerance, 'conv')

## Shaft Speed Time Response -------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 1, figsize = (14,8), edgecolor = 'k')

plt1 = plt.subplot(1,2,1) 

plt.plot(current_time,time_N,'b')

plt1.set_xlabel(r"$ t \ [s]$",loc = 'right',fontsize = 18)
plt1.set_ylabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='center',fontsize=20)
plt1.set_xlim([0,simulation_time])
plt1.set_ylim([12000,30000])

plt1.set_title(r"$\bf{CORRECTED \ SHAFT \ SPEED}$", fontsize = 12)

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

plt2 = plt.subplot(1,2,2) 

plt.plot(current_time,time_fuel_param,'k')

plt2.set_xlabel(r"$ t \ [s]$",loc='right',fontsize=18)
plt2.set_ylabel(r"$\frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='center',fontsize=20)
plt2.set_xlim([0,simulation_time])
plt2.set_ylim([2,6])

plt2.set_title(r"$\bf{FUEL \ PARAMETER}$", fontsize = 12)

plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

plt.suptitle(r"$\rm{TRANSIENT \ SIMULATION - SHAFT \ SPEED}$", fontsize = 14)
plt.show()

## Operating Lines ------------------------------------------------------------------------------------------------------------------------------------------
            
# Compressor -------------------------------------------------------------------------------------------------------------------------------------------------------

fig = plt.figure(num=1, figsize=(14,8), edgecolor='k')

choke = w[:,20]
choke = choke == 1

plt1 = plt.subplot(1,2,1)
componentPlot("C",True,'viridis',0.5)

plt1.plot(w[:,3],w[:,5],color="r",linestyle="--")
plt1.plot(w[:,3][choke],w[:,5][choke],color="r")

for i in range(n):

    if i%np.floor(n/20) == 0:
        plt1.text(w[i,3],w[i,5]-0.025,str(np.round(current_time[i],2)),size = 8,color = 'r')

plt1.set_title(r"$\bf{COMPRESSOR}$", fontsize = 12)

# Turbine ---------------------------------------------------------------------------------------------------------------------------------------------------

choke = w[:,20]
choke = choke == 1

plt2 = plt.subplot(1,2,2)
componentPlot("T",True,'viridis',0.5)

plt2.plot(w[:,14],1/w[:,16],color="r",linestyle="--")
plt2.plot(w[:,14][choke],1/w[:,16][choke],color="r")
plt2.set_xlim([12,12.6])
plt2.set_ylim([2.30,2.50])

for i in range(n):

    if i%np.floor(n/20) == 0:
        plt2.text(w[i,14],1/w[i,16],str(np.round(current_time[i],2)),size = 8,color = 'r')

plt2.set_title(r"$\bf{TURBINE}$", fontsize = 12)

plt.suptitle(r"$\rm{OPERATING \ LINES}$", fontsize = 14)
plt.show()
