from Solvers.OperatingPointSolver import engOperation
from Miscellaneous.AuxilliaryFunctions import componentPlot
from alive_progress import alive_bar
from Components.DesignVariables import p_ref, T_ref, Cp_c, gamma_e, R, A8

from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import numpy as np, csv

## ENGINE CHARACTERIZATION :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Operating Lines and Characteristic Curves

# Choose nozzle type ("conv" / "conv-div") -----------------------------------------------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------

nozzle_type = "conv"

# Choose accuracy and Mach numbers to plot -----------------------------------------------------------------------------------------------------------------

N_min, N_max =  61/108.5, 123/108.5
Num_points = 25

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

N_c = np.linspace(N_min,N_max,Num_points)
M0 = 0

map = {}
var = ["m0", "M0", "mc", "pic", "N", "mt", "pit", \
"T4tT0", "T5tT0", "M9", "A9A8", "p9p0", "T9T0", "E*", "Isp*", "TSFC*", "fuelParameter*", "choke"]

for x in var:
    map[x] = np.empty([Num_points])

with alive_bar(Num_points) as bar:

    for i in range(Num_points):

        m_0, T2t_T0, p2t_p0, eta_d, m_2, T3t_T2t, p3t_p2t, eta_c, \
        m_3, T4t_T3t, p4t_p3t, m_4, T41t_T4t, p41t_p4t, m_41, T5t_T41t, p5t_p41t, eta_t, N_t, m_5, \
        choked, T9_T5t, p9_p5t, eta_n, M9, A9_A8, p9_p0, T9_T0, E, Isp, TSFC, load_param, fuel_param, N   = \
        engOperation(M0, N_c[i], nozzle_type, 15, 0.35)

        # Distinguish between points in choke or non-choke conditions:

        map["m0"][i] = m_0
        map["M0"][i] = M0

        map["mc"][i] = m_2
        map["pic"][i] = p3t_p2t
        map["N"][i] = N

        map["T4tT0"][i] = T4t_T3t*T3t_T2t*T2t_T0
        map["T5tT0"][i] = T5t_T41t*T41t_T4t*T4t_T3t*T3t_T2t*T2t_T0

        map["mt"][i] = m_41
        map["pit"][i] = 1/p5t_p41t

        map["p9p0"][i] = p9_p0
        map["T9T0"][i] = T9_T0
        map["M9"][i] = M9
        map["A9A8"][i] = A9_A8

        # Characteristic curves:

        map["E*"][i] = E
        map["Isp*"][i] = Isp
        map["TSFC*"][i] = TSFC

        map["fuelParameter*"][i] = fuel_param
        map["choke"][i] = choked
        bar()

## Operating Lines ------------------------------------------------------------------------------------------------------------------------------------------
            
# Compressor ------------------------------------------------------------------------------------------------------------------------------------------------

fig = plt.figure(num = 1, figsize=(14,8), edgecolor='k')

plt1 = plt.subplot(1,2,1)
componentPlot("C",True,'viridis',0.5)

file2 = open("Custom/CSV Files/Validation/M-PR_2.csv", 'r')
data2 = np.double(np.array(list(csv.reader(file2))))

plt1.plot(data2[:,0],data2[:,1],'b--^')
plt1.plot(map["mc"],map["pic"],color="k",linewidth=2)

plt1.legend([r"$\rm{Current \ Model}$",r"$\rm{Ref. \ 2}$"], fontsize = 12)
plt1.set_title(r"$\bf{COMPRESSOR}$", fontsize = 12)

# Turbine ---------------------------------------------------------------------------------------------------------------------------------------------------

plt2 = plt.subplot(1,2,2)
componentPlot("T",True,'viridis',0.5)

plt2.plot(map["mt"],map["pit"],color="k",linewidth=2)

plt2.legend([r"$\rm{Current \ Model}$"], fontsize = 12)
plt2.set_title(r"$\bf{TURBINE}$", fontsize = 12)

plt.suptitle(r"$\rm{OPERATING \ LINES}$", fontsize = 14)
plt.show()

## Validation -----------------------------------------------------------------------------------------------------------------------------------------------

# (1) -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 2, figsize = (14,8), edgecolor = 'k')
plt1 = plt.subplot(1,2,1) 

file1 = open("Custom/CSV Files/Validation/T-RPM_1.csv", 'r')
data1 = np.double(np.array(list(csv.reader(file1))))
file2 = open("Custom/CSV Files/Validation/T-RPM_2.csv", 'r')
data2 = np.double(np.array(list(csv.reader(file2))))
file3 = open("Custom/CSV Files/Validation/T-RPM_3.csv", 'r')
data3 = np.double(np.array(list(csv.reader(file3))))

plt1.plot(map["N"],map["E*"],"k")
plt1.plot(data1[:,0],data1[:,1],'r--s')
plt1.plot(data2[:,0],data2[:,1],'b--^')
plt1.plot(data3[:,0],data3[:,1],'g--o')

plt1.set_title(r"$\bf{CORRECTED \ THRUST}$", fontsize = 12)
plt1.set_xlim(60000,110000)
plt1.set_ylim(50,250)

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

plt1.set_xlabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='right',fontsize=20)
plt1.set_ylabel(r"$\frac{E}{p_{\rm 0}/p_{\rm ref}} \ \left[\,\rm _{N}\,\right]$",loc='center',fontsize=20)
plt1.legend([r"$\rm{Current \ Model}$",r"$\rm{Ref. \ 1}$",r"$\rm{Ref. \ 2}$",r"$\rm{Ref. \ 3}$"], fontsize = 12)

plt2 = plt.subplot(1,2,2)
splineT = interp1d(map["N"],map["E*"])

interp_data1, error1 = [], []
interp_data2, error2 = [], []
interp_data3, error3 = [], []

for i in range(np.size(data1[:,0])):

    if data1[i,0] >= 61000 and data1[i,0] <= 108500:
        interp_data1.append(data1[i,0])
        error1.append(np.abs(splineT(data1[i,0])-data1[i,1])/data1[i,1]*100)

for i in range(np.size(data2[:,0])):

    if data2[i,0] >= 61000 and data2[i,0] <= 108500:
        interp_data2.append(data2[i,0])
        error2.append(np.abs(splineT(data2[i,0])-data2[i,1])/data2[i,1]*100)

for i in range(np.size(data3[:,0])):

    if data3[i,0] >= 61000 and data3[i,0] <= 108500:
        interp_data3.append(data3[i,0])
        error3.append(np.abs(splineT(data3[i,0])-data3[i,1])/data3[i,1]*100)

plt2.plot(interp_data1,error1,'r--s')
plt2.plot(interp_data2,error2,'b--^')
plt2.plot(interp_data3,error3,'g--o')

plt2.set_title(r"$\bf{CORRECTED \ THRUST \ RELATIVE \ ERROR}$", fontsize = 12)
plt2.set_xlim(60000,110000)
plt2.set_ylim(0,40)
plt2.legend([r"$\rm{Ref. \ 1}$",r"$\rm{Ref. \ 2}$",r"$\rm{Ref. \ 3}$"], fontsize = 12)

plt2.set_xlabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='right',fontsize=20)
plt2.set_ylabel(r"$\frac{E^* - E^*_{\rm m}}{E^*_{\rm m}} \ \left[ \% \right]$",\
loc='center',fontsize=20)

plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

fig.suptitle(r"$\rm{VALIDATION \ - (1)}$", fontsize = 14)
plt.show()

# (2) -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 3, figsize = (14,8), edgecolor = 'k')
plt1 = plt.subplot(1,2,1) 

file1 = open("Custom/CSV Files/Validation/M-RPM_1.csv", 'r')
data1 = np.double(np.array(list(csv.reader(file1))))
file2 = open("Custom/CSV Files/Validation/M-RPM_2.csv", 'r')
data2 = np.double(np.array(list(csv.reader(file2))))

plt1.plot(map["N"],map["m0"],"k")
plt1.plot(data1[:,0],data1[:,1]*np.sqrt((273.15+26)/288.15)/(102000/101325),'r--s')
plt1.plot(data2[:,0],data2[:,1],'b--^')

plt1.set_title(r"$\bf{CORRECTED \ MASS \ FLOW}$", fontsize = 12)
plt1.set_xlim(60000,110000)
plt1.set_ylim(0.2,0.6)
plt1.legend([r"$\rm{Current \ Model}$",r"$\rm{Ref. \ 1}$",r"$\rm{Ref. \ 2}$"], fontsize = 12)

plt1.set_xlabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='right',fontsize=20)
plt1.set_ylabel(r"$\frac{\it \dot m_{\rm 0} \sqrt{\!T_{\rm 0}/T_{\rm ref}}}{p_{\rm 0}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$",\
loc='center',fontsize=20)

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

plt2 = plt.subplot(1,2,2)
splineT = interp1d(map["N"],map["m0"])

interp_data1, error1 = [], []
interp_data2, error2 = [], []

for i in range(np.size(data1[:,0])):

    if data1[i,0] >= 61000 and data1[i,0] <= 108500:
        interp_data1.append(data1[i,0])
        error1.append(np.abs(splineT(data1[i,0])-data1[i,1])/data1[i,1]*100)

for i in range(np.size(data2[:,0])):

    if data2[i,0] >= 61000 and data2[i,0] <= 108500:
        interp_data2.append(data2[i,0])
        error2.append(np.abs(splineT(data2[i,0])-data2[i,1])/data2[i,1]*100)

plt2.plot(interp_data1,error1,'r--s')
plt2.plot(interp_data2,error2,'b--^')

plt2.set_title(r"$\bf{CORRECTED \ MASS \ FLOW \ RELATIVE \ ERROR}$", fontsize = 12)
plt2.set_xlim(60000,110000)
plt2.set_ylim(0,20)

plt2.set_xlabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='right',fontsize=20)
plt2.set_ylabel(r"$\frac{m^*_0 - m^*_{\rm 0,m}}{m^*_{\rm 0,m}} \ \left[ \% \right]$",\
loc='center',fontsize=20)
plt2.legend([r"$\rm{Ref. \ 1}$",r"$\rm{Ref. \ 2}$"], fontsize = 12)

plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

fig.suptitle(r"$\rm{VALIDATION \ - (2)}$", fontsize = 14)
plt.show()

# (3) -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 4, figsize = (14,8), edgecolor = 'k')
plt1 = plt.subplot(1,2,1) 

file1 = open("Custom/CSV Files/Validation/FF-RPM_1.csv", 'r')
data1 = np.double(np.array(list(csv.reader(file1))))
file2 = open("Custom/CSV Files/Validation/FF-RPM_2.csv", 'r')
data2 = np.double(np.array(list(csv.reader(file2))))
file3 = open("Custom/CSV Files/Validation/FF-RPM_3.csv", 'r')
data3 = np.double(np.array(list(csv.reader(file3))))

eta_cc = 0.95   # [-]
L = 42.8*10**6  # [J/kg]
T0 = 288.15     # [K]

plt1.plot(map["N"],map["fuelParameter*"]*Cp_c*T0/(eta_cc*L)*1000,"k")
plt1.plot(data1[:,0],data1[:,1]*1000,'r--s')
plt1.plot(data2[:,0],data2[:,1],'b--^')
plt1.plot(data3[:,0],data3[:,1],'g--o')

plt1.set_title(r"$\bf{FUEL \ CONSUMPTION}$", fontsize = 12)
plt1.set_xlim(60000,110000)
plt1.set_ylim(0,12)
plt1.legend([r"$\rm{Current \ Model}$",r"$\rm{Ref. \ 1}$",r"$\rm{Ref. \ 2}$",r"$\rm{Ref. \ 3}$"], fontsize = 12)

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

plt1.set_xlabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='right',fontsize = 18)
plt1.set_ylabel(r"$\it \dot m_{\rm f} \ \left[\,\frac{\rm kg}{\rm s}\right]$",loc='center',fontsize = 18)

plt2 = plt.subplot(1,2,2)
splineT = interp1d(map["N"],map["fuelParameter*"]*Cp_c*T0/(eta_cc*L)*1000)

interp_data1, error1 = [], []
interp_data2, error2 = [], []
interp_data3, error3 = [], []

for i in range(np.size(data1[:,0])):

    if data1[i,0] >= 61000 and data1[i,0] <= 108500:
        interp_data1.append(data1[i,0])
        error1.append(np.abs(splineT(data1[i,0])-data1[i,1]*1000)/(data1[i,1]*1000)*100)

for i in range(np.size(data2[:,0])):

    if data2[i,0] >= 61000 and data2[i,0] <= 108500:
        interp_data2.append(data2[i,0])
        error2.append(np.abs(splineT(data2[i,0])-data2[i,1])/data2[i,1]*100)

for i in range(np.size(data3[:,0])):

    if data3[i,0] >= 61000 and data3[i,0] <= 108500:
        interp_data3.append(data3[i,0])
        error3.append(np.abs(splineT(data3[i,0])-data3[i,1])/data3[i,1]*100)

plt2.plot(interp_data1,error1,'r--s')
plt2.plot(interp_data2,error2,'b--^')
plt2.plot(interp_data3,error3,'g--o')

plt2.set_title(r"$\bf{FUEL \ CONSUMPTION \ RELATIVE \ ERROR}$", fontsize = 12)
plt2.set_xlim(60000,110000)
plt2.set_ylim(0,20)

plt2.set_xlabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='right',fontsize=20)
plt2.set_ylabel(r"$\frac{\dot{m_{\rm f}} - \dot{m_{\rm f,m}}}{\dot{m_{\rm f,m}}} \ \left[ \% \right]$",\
loc='center',fontsize=20)
plt2.legend([r"$\rm{Ref. \ 1}$",r"$\rm{Ref. \ 2}$",r"$\rm{Ref. \ 3}$"], fontsize = 12)

plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

fig.suptitle(r"$\rm{VALIDATION \ - (3)}$", fontsize = 14)
plt.show()

# (4) -------------------------------------------------------------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(1,2, num = 5, figsize = (14,8), edgecolor = 'k')

plt1 = plt.subplot(1,2,1)
plt1.set_title(r"$\bf{EXHASUT \ GAS \ TEMPERATURE}$", fontsize = 12)

file1 = open("Custom/CSV Files/Validation/EGT-RPM_1.csv", 'r')
data1 = np.double(np.array(list(csv.reader(file1))))
file2 = open("Custom/CSV Files/Validation/EGT-RPM_2.csv", 'r')
data2 = np.double(np.array(list(csv.reader(file2))))

plt1.plot(map["N"],map["T5tT0"]*T0,"k")
plt1.plot(data1[:,0],data1[:,1],'r--s')
plt1.plot(data2[:,0],data2[:,1],'b--^')

plt1.set_xlim(60000,110000)
plt1.set_ylim(500,1000)

plt1.set_xlabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='right',fontsize=20)
plt1.set_ylabel(r"$\it T_{\rm 5t} \ \left[ K \right]$",\
loc='center',fontsize=18)
plt1.legend([r"$\rm{Current \ Model}$",r"$\rm{Ref. \ 1}$",r"$\rm{Ref. \ 2}$"], fontsize = 12)

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

plt2 = plt.subplot(1,2,2)
splineT = interp1d(map["N"],map["T5tT0"]*T0)

interp_data1, error1 = [], []
interp_data2, error2 = [], []

for i in range(np.size(data1[:,0])):

    if data1[i,0] >= 61000 and data1[i,0] <= 108500:
        interp_data1.append(data1[i,0])
        error1.append(np.abs(splineT(data1[i,0])-data1[i,1])/data1[i,1]*100)

for i in range(np.size(data2[:,0])):

    if data2[i,0] >= 61000 and data2[i,0] <= 108500:
        interp_data2.append(data2[i,0])
        error2.append(np.abs(splineT(data2[i,0])-data2[i,1])/data2[i,1]*100)

plt2.plot(interp_data1,error1,'r--s')
plt2.plot(interp_data2,error2,'b--^')

plt2.set_title(r"$\bf{EXHAUST \ GAS \ TEMPERATURE \ RELATIVE \ ERROR}$", fontsize = 12)
plt2.set_xlim(60000,110000)
plt2.set_ylim(0,20)

plt2.set_xlabel(r"$\frac{N_{\rm}}{\sqrt{T_{\rm 0}/T_{\rm ref}}}$",loc='right',fontsize=20)
plt2.set_ylabel(r"$\frac{T_{\rm 5t} - T^*_{\rm 5t}}{T^*_{\rm 5t}} \ \left[ \% \right]$",\
loc='center',fontsize=20)
plt2.legend([r"$\rm{Ref. \ 1}$",r"$\rm{Ref. \ 2}$"], fontsize = 12)

plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

fig.suptitle(r"$\rm{VALIDATION \ - (4)}$", fontsize = 14)
plt.show()

# (5) -------------------------------------------------------------------------------------------------------------------------------------------------------

# Discharge coefficient correction for bidimensional flow in conical nozzle:

k1, k2 = 4.75, 1.475

def cd(p9t_p9):

    CD_min, CD_max = 0.68, 0.974
    return (CD_max-CD_min)*(1-k2**(-k1*(p9t_p9-1))) + CD_min

def m(p9t_p9):

    M9 = np.sqrt(2/(gamma_e-1)*(p9t_p9**((gamma_e-1)/gamma_e)-1))

    return cd(p9t_p9)*A8*(p_ref*1e5)/np.sqrt(R*T_ref)*np.sqrt(gamma_e)*M9* \
    (1+(gamma_e-1)/2*M9**2)**(-(gamma_e+1)/(2*(gamma_e-1)))

file = open("Custom/CSV Files/Validation/NozzleCalibration.csv", 'r')
csvreader = csv.reader(file)
calib_nozzle = np.double(np.array(list(csvreader)))
   
m_calib_nozzle = calib_nozzle[:,1]/np.sqrt(T_ref)*(p_ref*10**5)
NPR_calib_nozzle = calib_nozzle[:,0]
NPR_max = ((gamma_e+1)/2)**(gamma_e/(gamma_e-1))

fig, ax = plt.subplots(1,2, num = 6,figsize=(14,8),edgecolor='k')
plt1 = plt.subplot(1,2,1) 
plt2 = plt.subplot(1,2,2)

plt1.plot(np.linspace(1,NPR_max,100),m(np.linspace(1,NPR_max,100)),'k-')
plt1.plot(NPR_calib_nozzle,m_calib_nozzle,'r--s')
plt1.legend([r"$\rm{Corrected \ One-Dimensional}$",r"$\rm{Experimental}$"], fontsize = 12)
plt1.set_title(r"$\bf{NOZZLE \ PRESSURE \ RATIO \ - \ CORRECTED \ ISENTROPIC \ MASS \ FLOW}$")

plt1.set_xlabel(r"$\frac{\it p_{\rm 5t}}{\it p_{\rm 9}} \ [-]$", \
loc='right',fontsize=20)
plt1.set_ylabel(r"$\frac{\it \dot m_{\rm 5} \sqrt{\!T_{\rm 5t}/T_{\rm ref}}}{p_{\rm 5t}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$", \
loc='center',fontsize=22)

plt1.set_xlim([1,2])
plt1.set_ylim([0,0.8])

plt1.grid(True,linewidth=0.25,color='k',which="major")
plt1.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt1.minorticks_on()

plt2.plot(np.linspace(1,10,100),cd(np.linspace(1,10,100)),'k-')

plt2.set_xlabel(r"$\frac{\it p_{\rm 5t}}{\it p_{\rm 9}} \ [-]$", \
loc='right',fontsize=22)
plt2.set_ylabel(r"$C_{\rm D} \ [-]$", \
loc='center',fontsize=18)

plt2.set_xlim([1,10])
plt2.set_ylim([0.5,1])

plt2.grid(True,linewidth=0.25,color='k',which="major")
plt2.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
plt2.minorticks_on()

plt2.set_title(r"$\bf{DISCHARGE \ COEFFICIENT - NOZZLE \ PRESSURE \ RATIO}$")
fig.suptitle(r"$\rm{NOZZLE \ CHARACTERISTIC \ CURVES}$", fontsize = 14)

plt.show()




