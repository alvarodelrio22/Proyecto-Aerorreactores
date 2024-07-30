from Solvers.OperatingPointSolver import engOperation
from Solvers.TransientOperationSolver import transientOperation
from Miscellaneous.AuxilliaryFunctions import componentPlot
import numpy as np, matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from Components.DesignVariables import fuel_param_design

## TRANSIENT CHARACTERIZATION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Iterative method (Halley Method) + Initial Value Problem Numerical Solver (with choice of propagator)
# Arguments: Initial State Conditions and Corrected Fuel Parameter Time Evolution

## Simulation parameters ------------------------------------------------------------------------------------------------------------------------------------

M0 = 0
N_LPC = 0.65

simulation_time = 20  # [s]
n = 300
num_iter0 = 150

I1 = 23.43500/15  # [kg·m^2]
I2 = 12.17470/15  # [kg·m^2]

nozzle = 'conv'

# Initial state definition: ---------------------------------------------------------------------------------------------------------------------------------

m_0, T2t_T0, p2t_p0, eta_d, m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, \
eta_LPT, N_LPT, m_5, choked, T9_T5t, p9_p5t, eta_n, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, load_param, fuel_param_0, N1_0, N2_0  = \
engOperation(M0,N_LPC,nozzle,15,0.15)

w_0 = M0, I1, I2, N1_0, N2_0

# Fuel parameter function definition: -----------------------------------------------------------------------------------------------------------------------

time_fuel_param = np.empty(n)
tolerance = np.empty(n)
relaxation_factor = np.empty(n)

for i in range(n):

    if i <= np.floor(n/20):

        time_fuel_param[i] = (fuel_param_design*i + fuel_param_0*(np.floor(n/20)-i))/np.floor(n/20)

        tolerance[i] = 5e-4
        relaxation_factor[i] = 0.85

    elif np.floor(n/20) < i < np.floor(n/2) :

        time_fuel_param[i] = fuel_param_design

        tolerance[i] = 5e-4
        relaxation_factor[i] = 0.85


    elif np.floor(n/2) <= i <= np.floor(11*n/20):

        time_fuel_param[i] = (fuel_param_0*(i-np.floor(n/2)) + \
        fuel_param_design*((np.floor(11*n/20)-np.floor(n/2))-(i-np.floor(n/2))))/(np.floor(11*n/20)-np.floor(n/2))

        tolerance[i] = 5e-4
        relaxation_factor[i] = 0.90

    else:

        time_fuel_param[i] = fuel_param_0

        tolerance[i] = 5e-4

        if np.floor(11*n/20) < i <= np.floor(11.25*n/20):

            relaxation_factor[i] = 0.90

        else:

            relaxation_factor[i] = 0.80

## Run Simulation ------------------------------------------------------------------------------------------------------------------------------------------

w, time_N1, time_N2, current_time = transientOperation(w_0, time_fuel_param, 'Euler', simulation_time/n, relaxation_factor, num_iter0, tolerance, 'conv')

## Shaft Speed Time Response -------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 1, figsize = (14,8), edgecolor = 'k')

plt1 = plt.subplot(1,2,1) 

plt.plot(current_time,time_N1,'b')
plt.plot(current_time,time_N2,'r')

plt1.set_xlabel(r"$ t \ [s]$",loc = 'right',fontsize = 18)
plt1.set_xlim([0,simulation_time])
plt1.set_ylim([9000,18000])

plt1.set_title(r"$\bf{CORRECTED \ SHAFT \ SPEEDS}$", fontsize = 12)

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

ybox1 = TextArea(r"$\ \ \ \left[\rm \, \frac{rev}{min} \right]$",textprops=dict(color='k',size=18,rotation=90,ha='left',va='bottom'))
ybox2 = TextArea(r"$\frac{N_{\rm HP}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",textprops=dict(color='r',size=18,rotation=90,ha='right',va='bottom'))
ybox3 = TextArea(r"$ \ \ \ \ \ \ , $",textprops=dict(color="k",size=14,rotation=90,ha='left',va='bottom'))
ybox4 = TextArea(r"$\frac{N_{\rm LP}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",textprops=dict(color='b',size=18,rotation=90,ha='left',va='bottom'))
ybox = VPacker(children=[ybox1, ybox2, ybox3, ybox4],align="bottom",pad=0,sep=5)
anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0, frameon=False, bbox_to_anchor=(-0.15, 0.25),\
bbox_transform = ax[0].transAxes, borderpad=0)
ax[0].add_artist(anchored_ybox)

plt2 = plt.subplot(1,2,2) 

plt.plot(current_time,time_fuel_param,'k')

plt2.set_xlabel(r"$ t \ [s]$",loc='right',fontsize=18)
plt2.set_ylabel(r"$\frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='center',fontsize=20)
plt2.set_xlim([0,simulation_time])
plt2.set_ylim([1,6])

plt2.set_title(r"$\bf{FUEL \ PARAMETER}$", fontsize = 12)

plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

plt.suptitle(r"$\rm{TRANSIENT \ SIMULATION - SHAFT \ SPEED}$", fontsize = 14)
plt.show()

## Operating Lines ------------------------------------------------------------------------------------------------------------------------------------------
            
# LPC -------------------------------------------------------------------------------------------------------------------------------------------------------

fig = plt.figure(num=1, figsize=(14,8), edgecolor='k')

choke = w[:,30]
choke = choke == 1

plt1 = plt.subplot(1,2,1)
componentPlot("LPC",True,'viridis',0.5)

plt1.plot(w[:,3],w[:,5],color="r",linestyle="--")
plt1.plot(w[:,3][choke],w[:,5][choke],color="r")

for i in range(n):

    if i%np.floor(n/20) == 0:
        plt1.text(w[i,3],w[i,5]-0.025,str(np.round(current_time[i],2)),size = 8,color = 'r')

plt1.set_title(r"$\bf{LOW \ PRESSURE \ COMPRESSOR}$", fontsize = 12)

# HPC -------------------------------------------------------------------------------------------------------------------------------------------------------

plt2 = plt.subplot(1,2,2)
componentPlot("HPC",True,'viridis',0.5)

plt2.plot(w[:,8],w[:,10],color="r",linestyle="--")
plt2.plot(w[:,8][choke],w[:,10][choke],color="r")

for i in range(n):

    if i%np.floor(n/20) == 0:
        plt2.text(w[i,8],w[i,10]-0.025,str(np.round(current_time[i],2)),size = 8,color = 'r')

plt2.set_title(r"$\bf{HIGH \ PRESSURE \ COMPRESSOR}$", fontsize = 12)

plt.suptitle(r"$\rm{OPERATING \ LINES \ - \ COMPRESSORS}$", fontsize = 14)
plt.show()

# HPT -------------------------------------------------------------------------------------------------------------------------------------------------------

fig = plt.figure(num=2, figsize=(14,8), edgecolor='k')

choke = w[:,30]
choke = choke == 1

plt1 = plt.subplot(1,2,1)
componentPlot("HPT",True,'viridis',0.5)

plt1.plot(w[:,19],1/w[:,21],color="r",linestyle="--")
plt1.plot(w[:,19][choke],1/w[:,21][choke],color="r")
plt1.set_xlim([13.5,14.1])
plt1.set_ylim([1.62,1.72])

for i in range(n):

    if i%np.floor(n/20) == 0:
        plt1.text(w[i,19],1/w[i,21],str(np.round(current_time[i],2)),size = 8,color = 'r')

plt1.set_title(r"$\bf{HIGH \ PRESSURE \ TURBINE}$", fontsize = 12)

# LPT -------------------------------------------------------------------------------------------------------------------------------------------------------

plt2 = plt.subplot(1,2,2)
componentPlot("LPT",True,'viridis',0.5)

plt2.plot(w[:,24],1/w[:,26],color="r",linestyle="--")
plt2.plot(w[:,24][choke],1/w[:,26][choke],color="r")
plt2.set_xlim([21.4,23.0])
plt2.set_ylim([1.25,1.45])

for i in range(n):

    if i%np.floor(n/20) == 0:
        plt2.text(w[i,24],1/w[i,26],str(np.round(current_time[i],2)),size = 8,color = 'r')

plt2.set_title(r"$\bf{LOW \ PRESSURE \ TURBINE}$", fontsize = 12)

plt.suptitle(r"$\rm{OPERATING \ LINES \ - \ TURBINES}$", fontsize = 14)
plt.show()
