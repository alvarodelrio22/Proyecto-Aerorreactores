from Components.ComponentMap import compressor, turbine
from scipy.optimize import newton
import numpy as np, csv, warnings, matplotlib.pyplot as plt, matplotlib as mpl

mpl.rcParams['mathtext.fontset'] = 'cm'
warnings.filterwarnings('ignore')

## DESIGN POINT DETERMINATION: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Iterating along component map definition

gamma_e = 1.344
R = 287
T_ref = 288.15
p_ref = 1.01325

## Nozzle Outlet (9) ----------------------------------------------------------------------------------------------------------------------------------------

# Designed to maximum expansion: choke conditions.

eta_n = 1

# As NPR > 5 a variable convergent - divergent segment is considered for the nozzle. A full expansion is considered for the design configuration.
# For NPR greater than 1.845, choke conditions happen for an engine embedded in ambient conditions. 
# No bigger NPR can be achieved without this phenomenon happening. It is assumed that the ambient conditions in the exit are p0, T0.

A8 = 0.0032 #  [m^2]

file = open("Custom/CSV Files/Validation/NozzleCalibration.csv", 'r')
csvreader = csv.reader(file)
calib_nozzle = np.double(np.array(list(csvreader)))

# Discharge coefficient correction for bidimensional flow in conical nozzle:

k1, k2 = 4.75, 1.475

def cd(p9t_p9):

    CD_min, CD_max = 0.68, 0.974
    return (CD_max-CD_min)*(1-k2**(-k1*(p9t_p9-1))) + CD_min

def m(p9t_p9):

    M9 = np.sqrt(2/(gamma_e-1)*(p9t_p9**((gamma_e-1)/gamma_e)-1))

    return cd(p9t_p9)*A8*(p_ref*1e5)/np.sqrt(R*T_ref)*np.sqrt(gamma_e)*M9* \
    (1+(gamma_e-1)/2*M9**2)**(-(gamma_e+1)/(2*(gamma_e-1)))
   
m_calib_nozzle = calib_nozzle[:,1]/np.sqrt(T_ref)*(p_ref*10**5)
NPR_calib_nozzle = calib_nozzle[:,0]
NPR_max = ((gamma_e+1)/2)**(gamma_e/(gamma_e-1))

fig, ax = plt.subplots(1,2,num=1,figsize=(14,8),edgecolor='k')
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





