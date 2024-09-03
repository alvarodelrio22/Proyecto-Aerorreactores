
import csv, matplotlib as mpl, numpy as np, warnings
from matplotlib import pyplot as plt

from DesignVariables import m_c_design, m_t_design, pi_c_design, pi_t_design
from scipy.interpolate import RectBivariateSpline

mpl.rcParams['mathtext.fontset'] = 'cm'
warnings.filterwarnings('ignore')

type = "T"

if type == "C":

    m_scale_factor = 1
    pi_scale_factor = 1

elif type == "T":

    m_scale_factor = 1/75
    pi_scale_factor = 2/3

## MAP PLOT: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
# Choose number of points to plot, beta lines, speed lines,
# and number of contours

## Load Generic data:

map_names = ["pi", "eta", "m"]
map = {}

for x in map_names:

    map[x] = []

    file = open("Custom/CSV Files/"+type+x+".csv", 'r')
    csvreader = csv.reader(file)
    data = list(csvreader)
    map[x] = np.array(data)

if type == "C":

    map["N"] = np.array([61,81,97,111,123])/108.5
    map["beta"] = np.linspace(0,1,13)
    precision = 3

elif type == "T":

    map["N"] = np.linspace(1/2,7/6,8)
    map["beta"] = np.linspace(0,1,11)

    map["N"]  = [0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1]
    map["beta"] = [0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2]
    precision = 3

# ----------------------------------------------------------------------------------------------------------------------------------------------------------
    
refined_map = {}
refined_map["m"] = []
refined_map["pi"] = []
refined_map["eta"] = []
plt.rcParams["figure.figsize"] = [12,7.5]

Num_refinement = 26
Num_contour = 100

beta_refined = np.linspace(np.min(map["beta"]),np.max(map["beta"]),Num_refinement)
N_refined = np.linspace(np.min(map["N"]),np.max(map["N"]),Num_refinement)
pos_N = np.linspace(1,np.size(N_refined),Num_refinement)
pos_beta = np.linspace(1,np.size(beta_refined),Num_refinement)  

beta_mesh, N_mesh = np.meshgrid(beta_refined,N_refined)
interp_m = RectBivariateSpline(map["beta"], map["N"], np.transpose(map["m"]), kx=precision, ky=precision)
interp_pi = RectBivariateSpline(map["beta"], map["N"], np.transpose(map["pi"]), kx=precision, ky=precision)
interp_eta = RectBivariateSpline(map["beta"], map["N"], np.transpose(map["eta"]), kx=precision, ky=precision)

for i in range(Num_refinement):
 
 refined_map["m"].append([])
 refined_map["pi"].append([])
 refined_map["eta"].append([])

 refined_map["m"][i] = m_scale_factor*interp_m.ev(beta_refined,N_refined[i]*np.ones(Num_refinement))
 refined_map["pi"][i] = pi_scale_factor*(interp_pi.ev(beta_refined,N_refined[i]*np.ones(Num_refinement))-1)+1
 refined_map["eta"][i] = interp_eta.ev(beta_refined,N_refined[i]*np.ones(Num_refinement))

if type == "C":

    print(np.max(refined_map["eta"]),np.min(refined_map["eta"]))

    comp_map = plt.contourf(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.6,0.8,Num_contour),cmap = 'jet')
    contour_lines = plt.contour(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.6,0.8,Num_contour),vmax = 0.6,
    cmap = "grey",linestyle=':',linewidths = 0.1)

    for i in range(Num_refinement-1):

        plt.plot(refined_map["m"][int(pos_N[i])],refined_map["pi"][int(pos_N[i])], color = 'k',linewidth = 0.6,linestyle = '-.')
        plt.text(refined_map["m"][int(pos_N[i])][0]-0.025,refined_map["pi"][int(pos_N[i])][0],str(np.round(N_refined[int(pos_N[i])],2)),size = 10)

    plt.grid(True,linewidth=0.25,color='k',which="major")
    plt.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
    plt.minorticks_on()

    plt.xlim([0,0.6])
    plt.ylim([1,5])

    xlabel = plt.xlabel(r"$\frac{\it \dot m_{\rm 2} \sqrt{\!T_{\rm 2t}/T_{\rm ref}}}{p_{\rm 2t}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$",loc='right')
    xlabel.set_fontsize(20)
    
    ylabel = plt.ylabel(r"$\pi_{\rm C} \ [-]$",loc='center',rotation=0)
    ylabel.set_fontsize(18)
    ylabel.set_horizontalalignment('right')

    plt.title(r"$\bf{COMPRESSOR \ MAP}$",fontsize = 16)
    plt.scatter(m_c_design,pi_c_design,100,marker='*',color='yellow',linewidth=1.25)

    Nlabel = plt.text(0.15,4.35,r"$\frac{\left(N_{\rm C}/N_{\rm C,ref}\right)}{\sqrt{\!T_{\rm 2t}/T_{\rm ref}}} \ [-]$")
    Nlabel.set_fontsize(18)

    colbar = plt.colorbar(comp_map)
    colbar.set_label(r"$\ \it η_{\rm C} \ [-]$",fontsize = 18, rotation=0, horizontalalignment='left')
    colbar.set_ticks(np.linspace(0.6,0.8,5))
    colbar.set_ticklabels(["0.60","0.65","0.70","0.75","0.80"])

    plt.plot([0.22, refined_map["m"][int(pos_N[15])][0]],
    [4.20, refined_map["pi"][int(pos_N[15])][0]],color='k',linewidth = 0.75)


elif type == "T":

    comp_map = plt.contourf(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.50,0.95,Num_contour),cmap = 'jet')
    contour_lines = plt.contour(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.50,0.95,Num_contour),vmax = 0.5,
    cmap = "grey",linestyle=':',linewidths = 0.1)

    for i in range(Num_refinement-1):

        plt.plot(refined_map["m"][int(pos_N[i])],refined_map["pi"][int(pos_N[i])], color = 'k',linewidth = 0.35,linestyle = '-.')
        plt.text(refined_map["m"][int(pos_N[i])][int(pos_beta[i])]-0.001,refined_map["pi"][int(pos_N[i])][int(pos_beta[i])],
        str(np.round(N_refined[int(pos_N[i])],2)),size = 8)

    plt.grid(True,linewidth=0.25,color='k',which="major")
    plt.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
    plt.minorticks_on()

    fig = plt.figure(num=1, figsize=(14,8), edgecolor='k')

    plt.xlim([0.14,0.28])
    plt.ylim([1,3])

    xlabel = plt.xlabel(r"$\frac{\it \dot m_{\rm 41} \sqrt{\!T_{\rm 41t}/T_{\rm ref}}}{p_{\rm 41t}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$",loc='right')
    xlabel.set_fontsize(20)
    
    ylabel = plt.ylabel(r"$\pi_{\rm T} \ [-]$",loc='center',rotation=0)
    ylabel.set_fontsize(18)
    ylabel.set_horizontalalignment('right')

    plt.title(r"$\bf{TURBINE \ MAP}$",fontsize = 16)
    plt.scatter(m_t_design,pi_t_design,100,marker='*',color='yellow',linewidth=1.25)

    Nlabel = plt.text(0.15,2.1,r"$\frac{\left(N_{\rm T}/N_{\rm T,ref}\right)}{\sqrt{\!T_{\rm 41t}/T_{\rm ref}}} \ [-]$")
    Nlabel.set_fontsize(18)

    colbar = plt.colorbar(comp_map)
    colbar.set_label(r"$\ \it η_{\rm T} \ [-]$",fontsize = 18, rotation=0, horizontalalignment='left')
    colbar.set_ticks(np.linspace(0.55,0.95,9))
    colbar.set_ticklabels(["0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.9","0.95"])

    plt.plot([0.16, refined_map["m"][-1][int(Num_refinement/4)]],
    [2.0, refined_map["pi"][-1][int(Num_refinement/4)]],color='k',linewidth = 0.75)


plt.show()




