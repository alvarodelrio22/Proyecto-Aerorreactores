
import csv, matplotlib as mpl, numpy as np, warnings
from matplotlib import pyplot as plt

from Components.DesignVariables import m_c_design, m_t_design, pi_c_design, pi_t_design
from scipy.interpolate import RectBivariateSpline

mpl.rcParams['mathtext.fontset'] = 'cm'
warnings.filterwarnings('ignore')

## PLOT COMPONENT MAP (EMPTY/NON-EMPTY) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def componentPlot(type,show,colmap,alphachannel):

    if type == "C":

        m_scale_factor = 1
        pi_scale_factor = 1

    elif type == "T":

        m_scale_factor = 1/64
        pi_scale_factor = 7/8

    # Not as indicated in GSP Manual.
    # Reference spool speed [rpm] of the component design point,
    # not necessarily coincident with the actual cycle design point,
    # referred to as simply "design".

    # Load Generic data:

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

        #map["N"] = np.linspace(1/2,7/6,8)
        #map["beta"] = np.linspace(0,1,11)

        map["N"] = [0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2]
        map["beta"] = [0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1]
        precision = 3

    # ----------------------------------------------------------------------------------------------------------------------------------------------------------
        
    refined_map = {}
    refined_map["m"] = []
    refined_map["pi"] = []
    refined_map["eta"] = []
    plt.rcParams["figure.figsize"] = [14,8]

    Num_refinement = 30

    beta_refined = np.linspace(np.min(map["beta"]),np.max(map["beta"]),Num_refinement)
    N_refined = np.linspace(np.min(map["N"]),np.max(map["N"]),Num_refinement)
    pos_N = np.linspace(1,np.size(N_refined),Num_refinement)
    pos_beta = np.linspace(1,np.size(beta_refined),Num_refinement)  

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

        if show:

            comp_map = plt.contourf(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.6,0.8,100),
            cmap=colmap,linestyle=':',linewidths=1,alpha=alphachannel)
            colbar = plt.colorbar(comp_map,orientation='horizontal',location='bottom',pad=-0.15,shrink=0.8)
    
        else:

            plt.contour(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.6,0.8,100),
            cmap = "grey",linestyle=':',linewidths = 0.1)

        for i in range(Num_refinement-1):

            plt.plot(refined_map["m"][int(pos_N[i])],refined_map["pi"][int(pos_N[i])], color = 'k',linewidth = 0.6,linestyle = '-.')
            plt.text(refined_map["m"][int(pos_N[i])][0]-0.035,refined_map["pi"][int(pos_N[i])][0],str(np.round(N_refined[int(pos_N[i])],2)),size = 10)

        plt.grid(True,linewidth=0.25,color='k',which="major")
        plt.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
        plt.minorticks_on()

        plt.xlim([0,0.6])
        plt.ylim([1,5])

        xlabel = plt.xlabel(r"$\frac{\it \dot m_{\rm 2} \sqrt{\!T_{\rm 2t}/T_{\rm ref}}}{p_{\rm 2t}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$",loc='right')
        xlabel.set_fontsize(20)
        
        if show:
            
            ylabel = plt.ylabel(r"$\pi_{\rm C} \ [-]$")
            ylabel.set_fontsize(18)

            Nlabel = plt.text(0.1,4.3,r"$\frac{\left(N_{\rm C}/N_{\rm C,ref}\right)}{\sqrt{\!T_{\rm 2t}/T_{\rm ref}}} \ [-]$")
            Nlabel.set_fontsize(18)

            colbar.set_label(r"$\eta_{\rm C} \ [-]$",fontsize = 16)
            colbar.set_ticks(np.linspace(0.6,0.8,5))
            colbar.set_ticklabels(["0.60","0.65","0.70","0.75","0.80"])

        else:

            plt.title(r"$\bf{COMPRESSOR \ MAP}$",fontsize = 16)

            ylabel = plt.ylabel(r"$\pi_{\rm C} \ [-]$",loc='center',rotation=0)
            ylabel.set_fontsize(18)
            ylabel.set_horizontalalignment('right')
            
            Nlabel = plt.text(0.1,2.85,r"$\frac{\left(N_{\rm C}/N_{\rm C,ref}\right)}{\sqrt{\!T_{\rm 2t}/T_{\rm ref}}} \ [-]$")
            Nlabel.set_fontsize(18)   

        plt.plot([0.22, refined_map["m"][int(pos_N[18])][0]],
        [4.20, refined_map["pi"][int(pos_N[18])][0]],color='k',linewidth = 0.75)


    elif type == "T":

        if show:

            turb_map = plt.contourf(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.5,0.85,100),\
            cmap=colmap,linestyle=':',linewidths=1,alpha=alphachannel)
            colbar = plt.colorbar(turb_map,orientation='vertical',location='left',pad =-0.2,shrink=0.7)

        else:

            plt.contour(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.5,0.85,100),vmax=0.5,\
            cmap = "grey",linestyle=':',linewidths = 0.1)

        for i in range(Num_refinement-1):

            plt.plot(refined_map["m"][int(pos_N[i])],refined_map["pi"][int(pos_N[i])], color = 'k',linewidth = 0.35,linestyle = '-.')
            plt.text(refined_map["m"][int(pos_N[i])][int(pos_beta[i])]-0.001,refined_map["pi"][int(pos_N[i])][int(pos_beta[i])],
            str(np.round(N_refined[int(pos_N[i])],2)),size = 8)

        plt.grid(True,linewidth=0.25,color='k',which="major")
        plt.grid(True,linewidth=0.15,linestyle=':',color='k',which="minor")
        plt.minorticks_on()

        plt.xlim([0.14,0.32])
        plt.ylim([1,3.5])

        xlabel = plt.xlabel(r"$\frac{\it \dot m_{\rm 41} \sqrt{\!T_{\rm 41t}/T_{\rm ref}}}{p_{\rm 41t}/p_{\rm ref}} \ \left[\,\frac{\rm kg}{\rm s}\right]$",loc='right')
        xlabel.set_fontsize(20)
        
        if show:

            ylabel = plt.ylabel(r"$\pi_{\rm T} \ [-]$",loc='center')
            ylabel.set_fontsize(18)

            Nlabel = plt.text(0.20,2.3,r"$\frac{\left(N_{\rm T}/N_{\rm T,ref}\right)}{\sqrt{\!T_{\rm 41t}/T_{\rm ref}}} \ [-]$")
            Nlabel.set_fontsize(18)

            colbar.set_label(r"$\eta_{\rm T} \ [-]$",fontsize = 16)
            colbar.set_ticks(np.linspace(0.5,0.85,8))
            colbar.set_ticklabels(["0.5","0.55","0.60","0.65","0.70","0.75","0.80","0.85"])

        else:

            plt.title(r"$\bf{TURBINE \ MAP}$",fontsize = 16)

            ylabel = plt.ylabel(r"$\pi_{\rm T} \ [-]$",loc='center',rotation=0)
            ylabel.set_fontsize(18)
            ylabel.set_horizontalalignment('right')

            Nlabel = plt.text(0.22,2.5,r"$\frac{\left(N_{\rm T}/N_{\rm T,ref}\right)}{\sqrt{\!T_{\rm 41t}/T_{\rm ref}}} \ [-]$")
            Nlabel.set_fontsize(18)
            Nlabel.set_fontsize(18)

        plt.plot([0.22, refined_map["m"][-1][int(Num_refinement/4)]],
        [2.2, refined_map["pi"][-1][int(Num_refinement/4)]],color='k',linewidth = 0.75)

    return plt

## RELAXATION FACTOR ASSIGNATION ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Define an auxilliary function for an adaptative relaxation factor:

def relaxationFactor(beta, N, beta_limits, N_limits, values):

    if np.isnan(beta) or np.isnan(N):
        return np.NaN
    elif beta < 0 or beta > 1 or N < 0.45 or N > 1.08:
        return np.NaN

    for i in range(len(beta_limits)-1):

        for j in range(len(N_limits)-1):

            if beta_limits[i] >= beta >= beta_limits[i+1] and N_limits[j] <= N <= N_limits[j+1]:

                relaxation_factor = values[i][j]
                return relaxation_factor