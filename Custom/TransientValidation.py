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

file = open("Custom/CSV Files/Validation/Measurements.csv", 'r')
data2 = np.double(np.array(list(csv.reader(file))))

simulation_time = data[-1,0] - data[0,0]
n = np.size(data[:,0])
time_fuel_param = data[:,1]

tolerance = 1e-2
relaxation_factor = 0.95

## Run Simulation ------------------------------------------------------------------------------------------------------------------------------------------

w, time_N, current_time = transientOperation(w_0, time_fuel_param, 'Euler', simulation_time/n, relaxation_factor, num_iter0, tolerance, 'conv')

## Fuel Input ----------------------------------------------------------------------------------------------------------------------------------------------

plt.figure(figsize = (14,8), edgecolor = 'k')
plt.plot(current_time,time_fuel_param,'b')

plt.xlabel(r"$ t \ [s]$",loc='right',fontsize=18)
plt.ylabel(r"$\frac{\eta_{\rm cc}\dot{m_{\rm f}}L}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='center',fontsize=20)
plt.xlim([0,simulation_time])
plt.ylim([0.5,1.3])

plt.title(r"$\bf{FUEL \ PARAMETER}$", fontsize = 12)

plt.grid(True,linewidth=0.25,color='k',which="major")
plt.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt.minorticks_on()

plt.show()

## Shaft Speed Time Response -------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 1, figsize = (14,8), edgecolor = 'k')

plt1 = plt.subplot(1,2,1) 

plt.plot(current_time,time_N,'m')

plt1.set_xlabel(r"$ t \ [s]$",loc = 'right',fontsize = 18)
plt1.set_ylabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}} [\rm \frac{rev}{min}]$",loc='center',fontsize=20)
plt1.set_xlim([0,simulation_time])
plt1.set_ylim([60000,110000])

plt1.set_title(r"$\bf{CORRECTED \ SHAFT \ SPEED}$", fontsize = 12)

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

## Thrust Response ---------------------------------------------------------------------------------------------------------------------------------------

plt2 = plt.subplot(1,2,2) 

plt2.plot(current_time,w[:,27],'m')
plt2.scatter(current_time,data2[:,6],20,c='k',marker='^')

plt2.set_xlabel(r"$ t \ [s]$",loc = 'right',fontsize = 18)
plt2.set_ylabel(r"$\frac{E_{\rm}}{p_{\rm 0}/p_{\rm ref}}$",loc='center',fontsize=20)
plt2.set_xlim([0,simulation_time])
plt2.set_ylim([100,220])

plt2.set_title(r"$\bf{CORRECTED \ THRUST}$", fontsize = 12)

#plt2.grid(True,linewidth=0.25,color='k',which="major")
#plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
#plt2.minorticks_on()

## Covariance Matrix -------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 2, figsize = (14,8), edgecolor = 'k')

plt1 = plt.subplot(1,2,1) 

plt1.plot(current_time,w[:,32],'m')

plt1.set_xlabel(r"$ t \ [s]$",loc='right',fontsize=18)
plt1.set_ylabel(r"$\left(\mathbf{P}\right)_{11} \ [\rm \left(\frac{rev}{min})\right)^2]$",loc='center',fontsize=16)
plt1.set_xlim([0,simulation_time])

plt1.set_title(r"$\bf{COVARIANCE \ MATRIX}$", fontsize = 12)

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

## Kalman Gain Evolution ----------------------------------------------------------------------------------------------------------------------------------

plt2 = plt.subplot(1,2,2) 

plt2.plot(current_time,w[:,33]/np.max(w[:,33]),'m')

plt2.set_xlabel(r"$ t \ [s]$",loc='right',fontsize=18)
plt2.set_ylabel(r"$\frac{|\mathbf{K}|}{|\mathbf{K}|_{\rm max}} \ [-]$",loc='center',fontsize=20)
plt2.set_xlim([0,simulation_time])

plt2.set_title(r"$\bf{RELATIVE \ KALMAN \ GAIN \ NORM}$", fontsize = 12)

plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

plt.show()

## Operating Lines ------------------------------------------------------------------------------------------------------------------------------------------
            
# Compressor -------------------------------------------------------------------------------------------------------------------------------------------------------

fig = plt.figure(num=1, figsize=(14,8), edgecolor='k')

choke = w[:,20]
choke = choke == 1

plt1 = plt.subplot(1,2,1)
#componentPlot("C",True,'viridis',0.5)

plt1.plot(w[:,3],w[:,5],color="m")
plt1.scatter(data[:,3],data[:,2],20,c='k',marker='^')

plt1.set_xlim(0,0.6)
plt1.set_ylim(1,5)

plt1.set_title(r"$\bf{COMPRESSOR}$", fontsize = 12)

# Turbine ---------------------------------------------------------------------------------------------------------------------------------------------------

choke = w[:,20]
choke = choke == 1

plt2 = plt.subplot(1,2,2)
componentPlot("T",True,'viridis',0.5)

plt2.plot(w[:,14],1/w[:,16],color="m",linestyle="--")
plt2.plot(w[:,14][choke],1/w[:,16][choke],color="m")

plt2.set_title(r"$\bf{TURBINE}$", fontsize = 12)

plt.suptitle(r"$\rm{OPERATING \ LINES}$", fontsize = 14)
plt.show()
