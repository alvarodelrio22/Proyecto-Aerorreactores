from Solvers.TransientOperationSolver import transientOperation
from Miscellaneous.AuxilliaryFunctions import componentPlot
import numpy as np, matplotlib.pyplot as plt, csv

## TRANSIENT CHARACTERIZATION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Iterative method (Halley Method) + Initial Value Problem Numerical Solver (with choice of propagator)
# Arguments: Initial State Conditions and Corrected Fuel Parameter Time Evolution

## Simulation parameters ------------------------------------------------------------------------------------------------------------------------------------

M0 = 0
I = 2e-5    #[kg·m^2]
num_iter0 = 1500

nozzle = 'conv'

# Initial state definition: ---------------------------------------------------------------------------------------------------------------------------------

N_0 = 60000
w_0 = M0, I, N_0

# Fuel parameter function definition: -----------------------------------------------------------------------------------------------------------------------

file = open("Custom/CSV Files/Validation/TransVal.csv", 'r')
data = np.double(np.array(list(csv.reader(file))))

simulation_time = data[-1,0] - data[0,0]
n = np.size(data[:,0])
time_fuel_param = data[:,1]

tolerance = 1e-2
relaxation_factor = 0.95

## Run Simulation ------------------------------------------------------------------------------------------------------------------------------------------

w, time_N, current_time = transientOperation(w_0, time_fuel_param, 'Euler', simulation_time/n, relaxation_factor, num_iter0, tolerance, 'conv')

## Shaft Speed Time Response -------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 1, figsize = (14,8), edgecolor = 'k')

plt1 = plt.subplot(1,2,1) 

plt.plot(current_time,time_N,'b')

plt1.set_xlabel(r"$ t \ [s]$",loc = 'right',fontsize = 18)
plt1.set_ylabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='center',fontsize=20)
plt1.set_xlim([0,simulation_time])
#plt1.set_ylim([12000,30000])

plt1.set_title(r"$\bf{CORRECTED \ SHAFT \ SPEED}$", fontsize = 12)

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

plt2 = plt.subplot(1,2,2) 

plt.plot(current_time,time_fuel_param,'k')

plt2.set_xlabel(r"$ t \ [s]$",loc='right',fontsize=18)
plt2.set_ylabel(r"$\frac{\eta_{\rm cc}\dot{m_{\rm f}}L}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='center',fontsize=20)
plt2.set_xlim([0,simulation_time])
#plt2.set_ylim([2,6])

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

plt1.plot(w[:,3],w[:,5],color="r")
plt1.scatter(data[:,3],data[:,2],20,c='k',marker='^')

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

for i in range(n):

    if i%np.floor(n/20) == 0:
        plt2.text(w[i,14],1/w[i,16],str(np.round(current_time[i],2)),size = 8,color = 'r')

plt2.set_title(r"$\bf{TURBINE}$", fontsize = 12)

plt.suptitle(r"$\rm{OPERATING \ LINES}$", fontsize = 14)
plt.show()
