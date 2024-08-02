from Solvers.OperatingPointSolver import engOperation
from Miscellaneous.AuxilliaryFunctions import componentPlot
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from alive_progress import alive_bar

from matplotlib import pyplot as plt
import numpy as np

## ENGINE CHARACTERIZATION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Operating Lines and Characteristic Curves

# Choose nozzle type ("conv" / "conv-div") -----------------------------------------------------------------------------------------------------------------

nozzle_type = "conv"

# Choose accuracy and Mach numbers to plot -----------------------------------------------------------------------------------------------------------------

N_min, N_max =  0.45, 1.08
Mach_min, Mach_max = 0, 1

Num_points = 25
Num_Mach = 5

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

N_LPC = np.linspace(N_min,N_max,Num_points)
M0 = np.linspace(Mach_min,Mach_max,Num_Mach)

map = {}
var = ["m0", "M0", "mLPC", "piLPC", "N1", "mLPT", "piLPT", \
"mHPC", "piHPC", "N2", "T4tT0", "OPR", "mHPT", "piHPT", "M9", \
"A9A8", "p9p0", "T9T0", "E*", "Isp*", "TSFC*", "fuelParameter*", "choke"]

for x in var:
    map[x] = np.empty([Num_points,Num_Mach])

with alive_bar(Num_points*Num_Mach) as bar:

    for i in range(Num_points):

        for j in range(Num_Mach):

            m_0, T2t_T0, p2t_p0, eta_d, m_2, T25t_T2t, p25t_p2t, eta_LPC, m_25, T3t_T25t, p3t_p25t, eta_HPC, N_HPC, \
            m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T45t_T41t, p45t_p41t, eta_HPT, N_HPT, m_45, T5t_T45t, p5t_p45t, \
            eta_LPT, N_LPT, m_5, choked, T9_T5t, p9_p5t, eta_n, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, load_param, fuel_param, N1, N2   = \
            engOperation(M0[j], N_LPC[i], nozzle_type, 15, 0.15)

            # Distinguish between points in choke or non-choke conditions:

            map["m0"][i,j] = m_0
            map["M0"][i,j] = M0[j]

            map["mLPC"][i,j] = m_2
            map["piLPC"][i,j] = p25t_p2t
            map["N1"][i,j] = N1

            map["mHPC"][i,j] = m_25
            map["piHPC"][i,j] = p3t_p25t
            map["OPR"][i,j] = p3t_p25t*p25t_p2t
            map["N2"][i,j] = N2

            map["T4tT0"][i,j] = T4t_T3t*T3t_T25t*T25t_T2t*T2t_T0
            map["mHPT"][i,j] = m_41
            map["piHPT"][i,j] = 1/p45t_p41t

            map["mLPT"][i,j] = m_45
            map["piLPT"][i,j] = 1/p5t_p45t

            map["p9p0"][i,j] = p9_p0
            map["T9T0"][i,j] = T9_T0
            map["M9"][i,j] = M9
            map["A9A8"][i,j] = A9_A8

            # Characteristic curves:

            map["E*"][i,j] = E
            map["Isp*"][i,j] = Isp
            map["TSFC*"][i,j] = TSFC

            map["fuelParameter*"][i,j] = fuel_param
            map["choke"][i,j] = choked
            bar()

## Operating Lines ------------------------------------------------------------------------------------------------------------------------------------------
            
# LPC -------------------------------------------------------------------------------------------------------------------------------------------------------

fig = plt.figure(num=1, figsize=(14,8), edgecolor='k')

plt1 = plt.subplot(1,2,1)
componentPlot("LPC",True,'viridis',0.5)

for j in range(Num_Mach):

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1 

    plt1.plot(map["mLPC"][:,j][~np.isnan(map["mLPC"][:,j])],map["piLPC"][:,j][~np.isnan(map["piLPC"][:,j])],color="r",linestyle="--")
    plt1.plot(map["mLPC"][:,j][~np.isnan(map["mLPC"][:,j])][choke],map["piLPC"][:,j][~np.isnan(map["piLPC"][:,j])][choke],color="r")

plt1.set_title(r"$\bf{LOW \ PRESSURE \ COMPRESSOR}$", fontsize = 12)

# HPC -------------------------------------------------------------------------------------------------------------------------------------------------------

plt2 = plt.subplot(1,2,2)
componentPlot("HPC",True,'viridis',0.5)

for j in range(Num_Mach):

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1  

    plt2.plot(map["mHPC"][:,j][~np.isnan(map["mHPC"][:,j])],map["piHPC"][:,j][~np.isnan(map["piHPC"][:,j])],color="r",linestyle="--")
    plt2.plot(map["mHPC"][:,j][~np.isnan(map["mHPC"][:,j])][choke],map["piHPC"][:,j][~np.isnan(map["piHPC"][:,j])][choke],color="r")

plt2.set_title(r"$\bf{HIGH \ PRESSURE \ COMPRESSOR}$", fontsize = 12)

plt.suptitle(r"$\rm{OPERATING \ LINES \ - \ COMPRESSORS}$", fontsize = 14)
plt.show()

# HPT -------------------------------------------------------------------------------------------------------------------------------------------------------

fig = plt.figure(num=2, figsize=(14,8), edgecolor='k')

plt1 = plt.subplot(1,2,1)
componentPlot("HPT",True,'viridis',0.5)

for j in range(Num_Mach):

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1  

    plt1.plot(map["mHPT"][:,j][~np.isnan(map["mHPT"][:,j])],map["piHPT"][:,j][~np.isnan(map["piHPT"][:,j])],\
    color="r",linestyle="--",linewidth=1.5)
    plt1.plot(map["mHPT"][:,j][~np.isnan(map["mHPT"][:,j])][choke],map["piHPT"][:,j][~np.isnan(map["piHPT"][:,j])][choke],\
    color="r",linewidth=1.5)

plt1.set_title(r"$\bf{HIGH \ PRESSURE \ TURBINE}$", fontsize = 12)

# LPT -------------------------------------------------------------------------------------------------------------------------------------------------------

plt2 = plt.subplot(1,2,2)
componentPlot("LPT",True,'viridis',0.5)

for j in range(Num_Mach):

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1  

    plt2.plot(map["mLPT"][:,j][~np.isnan(map["mLPT"][:,j])],map["piLPT"][:,j][~np.isnan(map["piLPT"][:,j])],\
    color="r",linestyle="--",linewidth=1.5)
    plt2.plot(map["mLPT"][:,j][~np.isnan(map["mLPT"][:,j])][choke],map["piLPT"][:,j][~np.isnan(map["piLPT"][:,j])][choke],\
    color="r",linewidth=1.5)

plt2.set_title(r"$\bf{LOW \ PRESSURE \ TURBINE}$", fontsize = 12)

plt.suptitle(r"$\rm{OPERATING \ LINES \ - \ TURBINES}$", fontsize = 14, weight = "bold")
plt.show()

## Characteristic Curves ------------------------------------------------------------------------------------------------------------------------------------

# (1) Corrected Thrust -  Corrected Inlet Mass Flow and Fuel Parameter --------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 3, figsize = (14,8), edgecolor = 'k')
plt1 = plt.subplot(1,2,1) 
plt2 = plt.subplot(1,2,2) 

pos = -1
halign = 'left'
valign = 'bottom'
shift11 = -1
shift12 = 0.5
shift21 = 0.05
shift22 = -0.45

for j in range(Num_Mach):

    if nozzle_type == 'conv-div' and j==4:
        shift11 = -5.5
        shift12 = -1
        shift21 = 0.05
        shift22 = -1

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1

    plt1 = plt.subplot(1,2,1)  

    plt1.plot(map["m0"][:,j][~np.isnan(map["m0"][:,j])],map["E*"][:,j][~np.isnan(map["E*"][:,j])]/1000,\
    color="b",linestyle="--")
    plt1.plot(map["m0"][:,j][~np.isnan(map["m0"][:,j])][choke],map["E*"][:,j][~np.isnan(map["E*"][:,j])][choke]/1000,\
    color="b")

    plt2 = plt.subplot(1,2,2)

    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["E*"][:,j][~np.isnan(map["E*"][:,j])]/1000,\
    color="b",linestyle="--")
    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["E*"][:,j][~np.isnan(map["E*"][:,j])][choke]/1000,\
    color="b")

    plt1.text(map["m0"][:,j][~np.isnan(map["m0"][:,j])][pos]+shift11,map["E*"][:,j][~np.isnan(map["E*"][:,j])][pos]/1000+shift12,str(M0[j]),\
    color='#FF007F',size=10,horizontalalignment=halign,verticalalignment=valign)
    plt2.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos]+shift21,map["E*"][:,j][~np.isnan(map["E*"][:,j])][pos]/1000+shift22,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign,verticalalignment=valign)

plt1.set_title(r"$\bf{CORRECTED \ THRUST}$", fontsize = 12)
plt2.set_title(r"$\bf{CORRECTED \ THRUST}$", fontsize = 12)

plt1.set_xlim(0,120)
plt2.set_xlim(1,7)

plt1.set_ylim(0,90)
plt2.set_ylim(0,90)

plt1.set_xlabel(r"$\frac{\it \dot m_{\rm 0} \sqrt{\!T_{\rm 0}/T_{\rm ref}}}{p_{\rm 0}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$",\
loc='right',fontsize=20)
plt1.set_ylabel(r"$\frac{E}{p_{\rm 0}/p_{\rm ref}} \ \left[\,\rm _{kN}\,\right]$",loc='center',fontsize=20)

plt2.set_xlabel(r"$\frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)
plt2.set_ylabel(r"$\frac{E}{p_{\rm 0}/p_{\rm ref}} \ \left[\,\rm _{kN}\,\right]$",loc='center',fontsize=20)

ax[0].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[0].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')
ax[1].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[1].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()
plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

fig.suptitle(r"$\rm{CHARACTERISTIC \ CURVES \ - (1)}$", fontsize = 14)
plt.show()

# (2) Corrected Specific Impulse -  Corrected Inlet Mass Flow and Fuel Parameter ----------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 4, figsize = (14,8), edgecolor = 'k')
plt1 = plt.subplot(1,2,1) 
plt2 = plt.subplot(1,2,2)   


if nozzle_type == 'conv-div':

    pos1 = 1
    pos2 = -1
    halign1 = 'center'
    valign1 = 'top'
    halign2 = 'left'
    valign2 = 'bottom'
    shift11 = 0
    shift12 = -4
    shift21 = 0
    shift22 = 1

else:

    pos1 = -1
    pos2 = -1
    halign1 = 'left'
    valign1 = 'bottom'
    halign2 = 'left'
    valign2 = 'bottom'
    shift11 = -2
    shift12 = 4
    shift21 = 0.05
    shift22 = 5

for j in range(Num_Mach):

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1

    plt1.plot(map["m0"][:,j][~np.isnan(map["m0"][:,j])],map["Isp*"][:,j][~np.isnan(map["Isp*"][:,j])],\
    color="b",linestyle="--")
    plt1.plot(map["m0"][:,j][~np.isnan(map["m0"][:,j])][choke],map["Isp*"][:,j][~np.isnan(map["Isp*"][:,j])][choke],\
    color="b")

    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["Isp*"][:,j][~np.isnan(map["Isp*"][:,j])],\
    color="b",linestyle="--")
    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["Isp*"][:,j][~np.isnan(map["Isp*"][:,j])][choke],\
    color="b")

    plt1.text(map["m0"][:,j][~np.isnan(map["m0"][:,j])][pos1]+shift11,map["Isp*"][:,j][~np.isnan(map["Isp*"][:,j])][pos1]+shift12,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign1,verticalalignment=valign1)
    plt2.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos2]+shift21,map["Isp*"][:,j][~np.isnan(map["Isp*"][:,j])][pos2]+shift22,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign2,verticalalignment=valign2)
    
plt1.set_title(r"$\bf{SPECIFIC \ IMPULSE}$", fontsize = 12)
plt2.set_title(r"$\bf{SPECIFIC \ IMPULSE}$", fontsize = 12)

plt1.set_xlim(0,120)
plt2.set_xlim(1,7)

plt1.set_ylim(0,1000)
plt2.set_ylim(0,1000)

plt1.set_xlabel(r"$\frac{\it \dot m_{\rm 0} \sqrt{\!T_{\rm 0}/T_{\rm ref}}}{p_{\rm 0}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$",\
loc='right',fontsize=20)
plt1.set_ylabel(r"$\frac{I_{\rm sp}}{\sqrt{T_{\rm 0}/T_{\rm ref}}} \ \left[\,\frac{\rm N}{\rm kg \cdot s} \,\right]$",\
loc='center',fontsize=18)

plt2.set_xlabel(r"$\frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)
plt2.set_ylabel(r"$\frac{I_{\rm sp}}{\sqrt{T_{\rm 0}/T_{\rm ref}}} \ \left[\,\frac{\rm N}{\rm kg \cdot s} \,\right]$",\
loc='center',fontsize=18)

ax[0].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[0].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')
ax[1].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[1].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()
plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

fig.suptitle(r"$\rm{CHARACTERISTIC \ CURVES \ - (2)}$", fontsize = 14)
plt.show()

# (3) Corrected Thrust Specific Fuel Consumption -  Corrected Inlet Mass Flow and Fuel Parameter ------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 5, figsize = (14,8), edgecolor = 'k')
plt1 = plt.subplot(1,2,1) 
plt2 = plt.subplot(1,2,2)

pos = -1
halign = 'left'
valign = 'bottom'
shift11 = 0
shift12 = -200
shift21 = 0.05
shift22 = -50

if nozzle_type == 'conv-div':
    shift22 = -100

for j in range(Num_Mach):

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1

    plt1 = plt.subplot(1,2,1)  

    plt1.plot(map["m0"][:,j][~np.isnan(map["m0"][:,j])],map["TSFC*"][:,j][~np.isnan(map["TSFC*"][:,j])],\
    color="b",linestyle="--")
    plt1.plot(map["m0"][:,j][~np.isnan(map["m0"][:,j])][choke],map["TSFC*"][:,j][~np.isnan(map["TSFC*"][:,j])][choke],\
    color="b")

    plt2 = plt.subplot(1,2,2)  

    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["TSFC*"][:,j][~np.isnan(map["TSFC*"][:,j])],\
    color="b",linestyle="--")
    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["TSFC*"][:,j][~np.isnan(map["TSFC*"][:,j])][choke],\
    color="b")

    plt1.text(map["m0"][:,j][~np.isnan(map["m0"][:,j])][pos]+shift11,map["TSFC*"][:,j][~np.isnan(map["TSFC*"][:,j])][pos]+shift12,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign,verticalalignment=valign)
    plt2.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos]+shift21,map["TSFC*"][:,j][~np.isnan(map["TSFC*"][:,j])][pos]+shift22,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign,verticalalignment=valign)

plt1.set_title(r"$\bf{CORRECTED \ THRUST \ SPECIFIC \ FUEL \ CONSUMPTION}$", fontsize = 12)
plt2.set_title(r"$\bf{CORRECTED \ THRUST \ SPECIFIC \ FUEL \ CONSUMPTION}$", fontsize = 12)

plt1.set_xlim(0,120)
plt2.set_xlim(1,7)

plt1.set_ylim(1000,8000)
plt2.set_ylim(1000,8000)

plt1.set_xlabel(r"$\frac{\it \dot m_{\rm 0} \sqrt{\!T_{\rm 0}/T_{\rm ref}}}{p_{\rm 0}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$",\
loc='right',fontsize=20)
plt1.set_ylabel(r"$\frac{\eta_{\rm cc}c_{\rm E}L}{\sqrt{T_{\rm 0}/T_{\rm ref}}} \ \left[\,\frac{\rm m}{\rm s} \,\right]$",\
loc='center',fontsize=18)

plt2.set_xlabel(r"$\frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)
plt2.set_ylabel(r"$\frac{\eta_{\rm cc}c_{\rm E}L}{\sqrt{T_{\rm 0}/T_{\rm ref}}} \ \left[\,\frac{\rm m}{\rm s} \,\right]$",\
loc='center',fontsize=18)

ax[0].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[0].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')
ax[1].text(0.925,0.96,r"$M_{0} \ [-]$",transform=ax[1].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()
plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

fig.suptitle(r"$\rm{CHARACTERISTIC \ CURVES \ - (3)}$", fontsize = 14)
plt.show()

# (4) Corrected Shaft Speeds  -  Corrected Inlet Mass Flow and Fuel Parameter -------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 6, figsize = (14,8), edgecolor = 'k')
plt1 = plt.subplot(1,2,1) 
plt2 = plt.subplot(1,2,2)

pos11 = 1
pos12 = -1
halign11 = 'right'
valign11 = 'bottom'
halign12 = 'left'
valign12 = 'top'
shift111 = 0.125
shift112 = -400
shift121 = 0.05
shift122 = 125

pos2 = -1
halign2 = 'left'
valign2 = 'bottom'
shift21 = 0.05
shift22 = -1.5

for j in range(Num_Mach):

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1

    plt1 = plt.subplot(1,2,1)  

    plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["N1"][:,j][~np.isnan(map["N1"][:,j])],\
    color="b",linestyle="--")
    plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["N1"][:,j][~np.isnan(map["N1"][:,j])][choke],\
    color="b")

    plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["N2"][:,j][~np.isnan(map["N2"][:,j])],\
    color="r",linestyle="--")
    plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["N2"][:,j][~np.isnan(map["N2"][:,j])][choke],
    color="r")

    plt2 = plt.subplot(1,2,2)

    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["m0"][:,j][~np.isnan(map["m0"][:,j])],\
    color="b",linestyle="--",)
    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["m0"][:,j][~np.isnan(map["m0"][:,j])][choke],
    color="b")

    plt1.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos11]+shift111,map["N1"][:,j][~np.isnan(map["N1"][:,j])][pos11]+shift112,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign11,verticalalignment=valign11)
    plt1.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos12]+shift121,map["N2"][:,j][~np.isnan(map["N2"][:,j])][pos12]+shift122,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign12,verticalalignment=valign12)
    plt2.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos2]+shift21,map["m0"][:,j][~np.isnan(map["m0"][:,j])][pos2]+shift22,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign2,verticalalignment=valign2)

plt1.set_title(r"$\bf{CORRECTED \ SHAFT \ SPEEDS}$", fontsize = 12)
plt2.set_title(r"$\bf{CORRECTED \ INLET \ MASS \ FLOW}$", fontsize = 12)

plt1.set_xlim(1,7)
plt2.set_xlim(1,7)

plt1.set_ylim(6000,20000)
plt2.set_ylim(0,120)

plt1.set_xlabel(r"$\frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)
ybox1 = TextArea(r"$\ \ \ \left[\rm \, \frac{rev}{min} \right]$",textprops=dict(color='k',size=18,rotation=90,ha='left',va='bottom'))
ybox2 = TextArea(r"$\frac{N_{\rm HP}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",textprops=dict(color='r',size=18,rotation=90,ha='right',va='bottom'))
ybox3 = TextArea(r"$ \ \ \ \ \ \ , $",textprops=dict(color="k",size=14,rotation=90,ha='left',va='bottom'))
ybox4 = TextArea(r"$\frac{N_{\rm LP}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",textprops=dict(color='b',size=18,rotation=90,ha='left',va='bottom'))
ybox = VPacker(children=[ybox1, ybox2, ybox3, ybox4],align="bottom",pad=0,sep=5)
anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0, frameon=False, bbox_to_anchor=(-0.15, 0.25),\
bbox_transform = ax[0].transAxes, borderpad=0)
ax[0].add_artist(anchored_ybox)

plt2.set_xlabel(r"$\frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)
plt2.set_ylabel(r"$\frac{\it \dot m_{\rm 0} \sqrt{\!T_{\rm 0}/T_{\rm ref}}}{p_{\rm 0}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$",\
loc='center',fontsize=18)

ax[0].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[0].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')
ax[1].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[1].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()
plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

fig.suptitle(r"$\rm{CHARACTERISTIC \ CURVES \ - (4)}$", fontsize = 14)
plt.show()

# (5) Inlet/Outlet Pressure Temperature Relation and Cycle OPR, Temperature Difference - Fuel Parameter -----------------------------------------------------

fig, ax = plt.subplots(1,2, num = 7, figsize = (14,8), edgecolor = 'k')
plt1 = plt.subplot(1,2,1) 
plt2 = plt.subplot(1,2,2)

pos11 = -1
pos12 = -1
pos22 = -1
halign11 = 'left'
valign11 = 'center'
halign12 = 'left'
valign12 = 'center'
halign2 = 'center'
valign2 = 'bottom'
shift111 = 0.05
shift112 = 0
shift121 = 0.05
shift122 = 0
shift211 = 0
shift221 = 0.22
shift222 = -0.1

if nozzle_type == "conv-div":
    shift121 = 0
    shift122 = 0.075

for j in range(Num_Mach):

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1

    if j%2 == 0:
        pos21 = -1
        shift212 = 0.05
    else:
        pos21 = 1
        shift212 = -0.25

    plt1 = plt.subplot(1,2,1)  

    plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["p9p0"][:,j][~np.isnan(map["p9p0"][:,j])],\
    color="b",linestyle="--")
    plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["p9p0"][:,j][~np.isnan(map["p9p0"][:,j])][choke],\
    color="b")

    plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["T9T0"][:,j][~np.isnan(map["T9T0"][:,j])],\
    color="r",linestyle="--")
    plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["T9T0"][:,j][~np.isnan(map["T9T0"][:,j])][choke],\
    color="r")

    plt2 = plt.subplot(1,2,2)

    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["OPR"][:,j][~np.isnan(map["OPR"][:,j])],\
    color="b",linestyle="--")
    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["OPR"][:,j][~np.isnan(map["OPR"][:,j])][choke],\
    color="b")

    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["T4tT0"][:,j][~np.isnan(map["T4tT0"][:,j])],\
    color="r",linestyle="--")
    plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["T4tT0"][:,j][~np.isnan(map["T4tT0"][:,j])][choke],
    color="r")

    if nozzle_type != "conv-div":
        plt1.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos11]+shift111,\
        map["p9p0"][:,j][~np.isnan(map["p9p0"][:,j])][pos11]+shift112,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign11,verticalalignment=valign11)
    plt1.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos12]+shift121,\
    map["T9T0"][:,j][~np.isnan(map["T9T0"][:,j])][pos12]+shift122,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign12,verticalalignment=valign12)
    plt2.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos21]+shift211,\
    map["OPR"][:,j][~np.isnan(map["OPR"][:,j])][pos21]+shift212,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign2,verticalalignment=valign2)
    plt2.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos22]+shift221,\
    map["T4tT0"][:,j][~np.isnan(map["T4tT0"][:,j])][pos22]+shift222,\
    str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign2,verticalalignment=valign2)

plt1.set_title(r"$\bf{INLET/OUTLET \ PRESSURE \ AND \ TEMPERATURE \ RATIO}$", fontsize = 12)
plt2.set_title(r"$\bf{CYCLE \ OPR \ AND \ MAXIMUM \ TEMPERATURE \ RATIO}$", fontsize = 12)

plt1.set_xlim(1,7)
plt2.set_xlim(1,7)

plt1.set_ylim(0,6)
plt2.set_ylim(0,12)

plt1.set_xlabel(r"$\frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)
plt2.set_xlabel(r"$\frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)

ybox1 = TextArea(r"$\ [-]$",textprops=dict(color='k',size=18,rotation=90,ha='left',va='bottom'))
ybox2 = TextArea(r"$\frac{T_{\rm 9}}{T_{\rm 0}}$",textprops=dict(color='r',size=18,rotation=90,ha='left',va='bottom'))
ybox3 = TextArea(r"$ \ , $",textprops=dict(color="k",size=14,rotation=90,ha='left',va='bottom'))
ybox4 = TextArea(r"$\frac{p_{\rm 9}}{p_{\rm 0}} \ $",textprops=dict(color='b',size=20,rotation=90,ha='left',va='bottom'))
ybox = VPacker(children=[ybox1, ybox2, ybox3, ybox4],align="bottom",pad=0,sep=5)
anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0, frameon=False, bbox_to_anchor=(-0.08, 0.4),\
bbox_transform=ax[0].transAxes, borderpad=0)
ax[0].add_artist(anchored_ybox)

ybox1 = TextArea(r"$\ [-]$",textprops=dict(color='k',size=18,rotation=90,ha='left',va='bottom'))
ybox2 = TextArea(r"$\frac{T_{\rm 4t}}{T_{\rm 0}}$",textprops=dict(color='r',size=18,rotation=90,ha='center',va='bottom'))
ybox3 = TextArea(r"$ \ \ , $",textprops=dict(color="k",size=14,rotation=90,ha='left',va='bottom'))
ybox4 = TextArea(r"$\frac{p_{\rm 3t}}{p_{\rm 2t}}$",textprops=dict(color='b',size=20,rotation=90,ha='left',va='bottom'))
ybox = VPacker(children=[ybox1, ybox2, ybox3, ybox4],align="bottom",pad=0,sep=5)
anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0, frameon=False, bbox_to_anchor=(-0.08, 0.4),\
bbox_transform=ax[1].transAxes, borderpad=0)
ax[1].add_artist(anchored_ybox)

ax[0].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[0].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')
ax[1].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[1].transAxes,color='#FF007F',size=14,\
bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()
plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

fig.suptitle(r"$\rm{CYCLE \ PROPERTIES \ - (1)}$", fontsize = 14)
plt.show()

# (6) Outlet Mach Number - Nozzle Area Relation -------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 8, figsize = (14,8), edgecolor = 'k')
plt1 = plt.subplot(1,2,1) 
plt2 = plt.subplot(1,2,2)

if nozzle_type == 'conv-div':

    pos1 = -1
    pos2 = -1
    halign1 = 'left'
    valign1 = 'bottom'
    halign2 = 'left'
    valign2 = 'center'
    shift11 = 0.05
    shift12 = -0.025
    shift21 = 0.05
    shift22 = 0

else:

    pos1 = 1
    pos2 = -1
    halign1 = 'center'
    valign1 = 'top'
    halign2 = 'left'
    valign2 = 'bottom'
    shift11 = 0
    shift12 = -0.005
    shift21 = 0
    shift22 = 0

for j in range(Num_Mach):

    choke = map["choke"][:,j][~np.isnan(map["choke"][:,j])]
    choke = choke == 1

    if nozzle_type == "conv":

        plt1 = plt.subplot(1,2,1)  
        
        plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["M9"][:,j][~np.isnan(map["M9"][:,j])],color="b",linestyle="--")
        plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["M9"][:,j][~np.isnan(map["M9"][:,j])][choke],color="b")

        plt2 = plt.subplot(1,2,2)

        plt1.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos1]+shift11,map["M9"][:,j][~np.isnan(map["M9"][:,j])][pos1]+shift12,\
        str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign1,verticalalignment=valign1)

    elif nozzle_type == "conv-div":

        plt1 = plt.subplot(1,2,1) 

        plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["M9"][:,j][~np.isnan(map["M9"][:,j])],color="b",linestyle="--")
        plt1.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["M9"][:,j][~np.isnan(map["M9"][:,j])][choke],color="b")

        plt1.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos1]+shift11,map["M9"][:,j][~np.isnan(map["M9"][:,j])][pos1]+shift12,\
        str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign1,verticalalignment=valign1)

        plt2 = plt.subplot(1,2,2)

        plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])],map["A9A8"][:,j][~np.isnan(map["A9A8"][:,j])],\
        color="b",linestyle="--")
        plt2.plot(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][choke],map["A9A8"][:,j][~np.isnan(map["A9A8"][:,j])][choke],\
        color="b")

        plt1.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos1]+shift11,map["M9"][:,j][~np.isnan(map["M9"][:,j])][pos1]+shift12,\
        str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign1,verticalalignment=valign1)
        plt2.text(map["fuelParameter*"][:,j][~np.isnan(map["fuelParameter*"][:,j])][pos2]+shift21,map["A9A8"][:,j][~np.isnan(map["A9A8"][:,j])][pos2]+shift22,\
        str(M0[j]),color='#FF007F',size=10,horizontalalignment=halign2,verticalalignment=valign2)

if nozzle_type == "conv":

    plt1.set_title(r"$\bf{OUTLET \ MACH \ NUMBER}$", fontsize = 12)
    plt2.set_title(r"$\bf{NOZZLE \ AREA \ RELATION}$", fontsize = 12)

    plt1.set_xlabel(r"$ \frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)
    plt1.set_ylabel(r"$M_{9} \ [-]$",loc='center',fontsize=16)

    plt1.set_xlim(1,7)
    plt1.set_ylim(0.4,1.1)

    plt1.grid(True,linewidth=0.25,color='k',which="major")
    plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
    plt1.minorticks_on()

    ax[0].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[0].transAxes,color='#FF007F',size=14,\
    bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')
    ax[1].text(0.5,0.5,"NOT AVAILABLE",transform=ax[1].transAxes,color='k',size=18,\
    bbox=dict(facecolor='w'),verticalalignment='center',horizontalalignment='center')

elif nozzle_type == "conv-div":

    ax[0].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[0].transAxes,color='#FF007F',size=14,\
    bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')
    ax[1].text(0.075,0.96,r"$M_{0} \ [-]$",transform=ax[1].transAxes,color='#FF007F',size=14,\
    bbox=dict(facecolor='w',edgecolor='k',linewidth=0.5),verticalalignment='center',horizontalalignment='center')

    plt1.set_title(r"$\bf{OUTLET \ MACH \ NUMBER}$", fontsize = 12)
    plt2.set_title(r"$\bf{NOZZLE \ AREA \ RELATION}$", fontsize = 12)

    plt1.set_xlim(1,7)
    plt2.set_xlim(1,7)

    plt1.set_ylim(0,2.5)
    plt2.set_ylim(0.8,2)

    plt1.set_xlabel(r"$ \frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)
    plt1.set_ylabel(r"$M_{9} \ [-]$",loc='center',fontsize=16)
    
    plt2.set_xlabel(r"$ \frac{\eta_{\rm cc}fL}{C_{\rm pc}T_{\rm 0}} \ [-]$",loc='right',fontsize=20)
    plt2.set_ylabel(r"$\frac{A_{9}}{A_{8}} \ [-]$",loc='center',fontsize=18)

    plt1.grid(True,linewidth=0.25,color='k',which="major")
    plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
    plt1.minorticks_on()
    plt2.grid(True,linewidth=0.25,color='k',which="major")
    plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
    plt2.minorticks_on()

fig.suptitle(r"$\rm{CYCLE \ PROPERTIES \ - (2)}$", fontsize = 14)
plt.show()



