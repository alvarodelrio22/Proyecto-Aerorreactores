
import csv, matplotlib as mpl, math as m, numpy as np
from matplotlib import pyplot as plt
from Components.DesignVariables import m_HPC_design, m_LPC_design, m_HPT_design, m_LPT_design, \
pi_LPC_design, pi_HPC_design, pi_HPT_design, pi_LPT_design
from scipy.interpolate import RectBivariateSpline

def componentPlot(type,show):

    m_c_ref = 19.7890   #[kg/s]
    pi_c_ref = 6.6631   #[bar]
    m_t_ref = 19.8100  #[kg/s]
    pi_t_ref = 2.5000  #[bar]

    if type == "HPC":

        m_scale_factor = m_HPC_design/m_c_ref
        pi_scale_factor = (pi_HPC_design - 1)/(pi_c_ref - 1)
        label = "C"

    elif type == "LPC":

        m_scale_factor = m_LPC_design/m_c_ref
        pi_scale_factor = (pi_LPC_design - 1)/(pi_c_ref - 1)
        label = "C"   

    elif type == "HPT":

        m_scale_factor = m_HPT_design/m_t_ref
        pi_scale_factor = (pi_HPT_design - 1)/(pi_t_ref - 1)
        label = "T"

    elif type == "LPT":

        m_scale_factor = m_LPT_design/m_t_ref
        pi_scale_factor = (pi_LPT_design - 1)/(pi_t_ref - 1)    
        label = "T"

    # Not as indicated in GSP Manual.
    # Reference spool speed [rpm] of the component design point,
    # not necessarily coincident with the actual cycle design point,
    # referred to as simply "design".

    ## Map plot
        
    # Choose number of points to plot, beta lines, speed lines,
    # and number of contours

    ## Load Generic data:

    map_names = ["pi", "eta", "m"]
    map = {}

    for x in map_names:

        map[x] = []

        file = open("CSV Files/"+label+x+"Gen.csv", 'r')
        csvreader = csv.reader(file)
        data = list(csvreader)
        map[x] = np.array(data)

    if label == "C":

        map["N"] = [0.45,0.5,0.6,0.7,0.8,0.85,0.90,0.92,0.94,0.955,0.98,1,1.04,1.08]
        map["beta"] = [0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1]

    elif label == "T":

        map["N"] = [0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2]
        map["beta"] = [0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1]

    # ----------------------------------------------------------------------------------------------------------------------------------------------------------
        
    refined_map = {}
    refined_map["m"] = []
    refined_map["pi"] = []
    refined_map["eta"] = []
    plt.rcParams["figure.figsize"] = [12,7.5]

    Num_refinement = 30
    Num_contour = 100

    beta_refined = np.linspace(np.min(map["beta"]),np.max(map["beta"]),Num_refinement)
    N_refined = np.linspace(np.min(map["N"]),np.max(map["N"]),Num_refinement)
    pos_N = np.linspace(1,np.size(N_refined),Num_refinement)
    pos_beta = np.linspace(1,np.size(beta_refined),Num_refinement)  

    beta_mesh, N_mesh = np.meshgrid(beta_refined,N_refined)
    interp_m = RectBivariateSpline(map["beta"], map["N"], np.transpose(map["m"]), kx=3, ky=3)
    interp_pi = RectBivariateSpline(map["beta"], map["N"], np.transpose(map["pi"]), kx=3, ky=3)
    interp_eta = RectBivariateSpline(map["beta"], map["N"], np.transpose(map["eta"]), kx=3, ky=3)

    for i in range(Num_refinement):
    
        refined_map["m"].append([])
        refined_map["pi"].append([])
        refined_map["eta"].append([])

        refined_map["m"][i] = m_scale_factor*interp_m.ev(beta_refined,N_refined[i]*np.ones(Num_refinement))
        refined_map["pi"][i] = pi_scale_factor*(interp_pi.ev(beta_refined,N_refined[i]*np.ones(Num_refinement))-1)+1
        refined_map["eta"][i] = interp_eta.ev(beta_refined,N_refined[i]*np.ones(Num_refinement))

    if label == "C":

        if show == True:
            comp_map = plt.contourf(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.5,0.9,100),cmap = 'jet')
            colbar = plt.colorbar(comp_map)

        plt.contour(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.5,0.9,100),vmax = 0.5,
        cmap = "grey",linestyle=':',linewidths = 0.2)

        for i in range(Num_refinement-1):

            plt.plot(refined_map["m"][int(pos_N[i])],refined_map["pi"][int(pos_N[i])], color = 'k',linewidth = 0.8,linestyle = '-.')
            plt.text(refined_map["m"][int(pos_N[i])][1],refined_map["pi"][int(pos_N[i])][1]-0.2,str(np.round(N_refined[int(pos_N[i])],2)),size = 10)

        plt.grid(linewidth = 0.1, color = 'k')

    elif label == "T":

        if show == True:
            comp_map = plt.contourf(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.5,0.9,100),cmap = 'jet')
            colbar = plt.colorbar(comp_map)
        
        plt.contour(refined_map["m"],refined_map["pi"],refined_map["eta"],np.linspace(0.5,1,100),vmax = 0.5,
        cmap = "grey",linestyle=':',linewidths = 0.1)

        for i in range(Num_refinement-1):

            plt.plot(refined_map["m"][int(pos_N[i])],refined_map["pi"][int(pos_N[i])], color = 'k',linewidth = 0.35,linestyle = '-.')
            plt.text(refined_map["m"][int(pos_N[i])][int(pos_beta[i])]-0.05,refined_map["pi"][int(pos_N[i])][int(pos_beta[i])]-0.01,
            str(np.round(N_refined[int(pos_N[i])],2)),size = 8)

        plt.grid(linewidth = 0.1, color = 'k')


    if type == "HPC":

        plt.xlim([np.floor(np.min(refined_map["m"])),np.ceil(np.max(refined_map["m"]))])
        plt.ylim([np.floor(np.min(refined_map["pi"])),np.ceil(np.max(refined_map["pi"]))])

        xlabel = plt.xlabel(r"$\frac{\it G_{\rm 25} \sqrt{T_{\rm 25t}/T_{\rm ref}}}{p_{\rm 25t}/p_{\rm ref}} \ \left[\frac{\rm kg}{\rm s}\right]$",loc='center')
        xlabel.set_fontsize(17)
        
        ylabel = plt.ylabel(r"$\pi_{\rm HPC}$",loc='center')
        ylabel.set_fontsize(17)

        plt.title('HIGH PRESSURE COMPRESSOR MAP',fontsize = 14, weight = 'bold')
        plt.scatter(m_HPC_design,pi_HPC_design,100,marker='*',color='yellow',linewidth=1.25)

        Nlabel = plt.text(np.ceil(np.max(refined_map["m"]))/2.5,np.ceil(np.max(refined_map["pi"]))*(4.75/8),
        r"$\frac{(N_{\rm HPC}/N_{\rm HPC,ref})}{\sqrt{T_{\rm 25t}/T_{\rm ref}}}$", fontsize=18)
        Nlabel.set_fontsize(15)

        if show == True:
            colbar.set_label(r"$\it η_{\rm HPC} \ \rm \left[-\right]$",fontsize = 14)
            colbar.set_ticks(np.linspace(0.5,0.9,9))
            colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

        plt.plot([np.ceil(np.max(refined_map["m"]))*(9.5/20), refined_map["m"][int(pos_N[14])][Num_refinement-1]],
        [np.ceil(np.max(refined_map["pi"]))*(4.75/8), refined_map["pi"][int(pos_N[14])][Num_refinement-1]],color='k',linewidth = 0.75)

    elif type == "LPC":

        plt.xlim([np.floor(np.min(refined_map["m"])),np.ceil(np.max(refined_map["m"]))])
        plt.ylim([np.floor(np.min(refined_map["pi"])),np.ceil(np.max(refined_map["pi"]))+1])

        xlabel = plt.xlabel(r"$\frac{\it G_{\rm 2} \sqrt{T_{\rm 2t}/T_{\rm ref}}}{p_{\rm 2t}/p_{\rm ref}} \ \left[\frac{\rm kg}{\rm s}\right]$",loc='center')
        xlabel.set_fontsize(17)
        
        ylabel = plt.ylabel(r"$\pi_{\rm LPC}$",loc='center')
        ylabel.set_fontsize(17)

        plt.title('LOW PRESSURE COMPRESSOR MAP',fontsize = 14, weight = 'bold')
        plt.scatter(m_LPC_design,pi_LPC_design,100,marker='*',color='yellow',linewidth=1.25)

        Nlabel = plt.text(np.ceil(np.max(refined_map["m"]))/2.5,np.ceil(np.max(refined_map["pi"]))*(4.75/8),
        r"$\frac{(N_{\rm HPC}/N_{\rm HPC,ref})}{\sqrt{T_{\rm 25t}/T_{\rm ref}}}$", fontsize=18)
        Nlabel.set_fontsize(15)

        if show == True:
            colbar.set_label(r"$\it η_{\rm LPC} \ \rm \left[-\right]$",fontsize = 14)
            colbar.set_ticks(np.linspace(0.5,0.9,9))
            colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

        plt.plot([np.ceil(np.max(refined_map["m"]))*(9.5/20), refined_map["m"][int(pos_N[14])][Num_refinement-1]],
        [np.ceil(np.max(refined_map["pi"]))*(4.75/8), refined_map["pi"][int(pos_N[14])][Num_refinement-1]],color='k',linewidth = 0.75)

    elif type == "HPT":

        plt.xlim([np.floor(np.min(refined_map["m"])),np.ceil(np.max(refined_map["m"]))])
        plt.ylim([np.floor(np.min(refined_map["pi"])),np.ceil(np.max(refined_map["pi"]))-0.8])

        xlabel = plt.xlabel(r"$\frac{\it G_{\rm 4} \sqrt{T_{\rm 4t}/T_{\rm ref}}}{p_{\rm 4t}/p_{\rm ref}} \ \left[\frac{\rm kg}{\rm s}\right]$",loc='center')
        xlabel.set_fontsize(17)
        
        ylabel = plt.ylabel(r"$\pi_{\rm HPT}$",loc='center')
        ylabel.set_fontsize(17)

        plt.title('HIGH PRESSURE TURBINE MAP',fontsize = 14, weight = 'bold')
        plt.scatter(m_HPT_design,pi_HPT_design,100,marker='*',color='yellow',linewidth=1.25)

        Nlabel = plt.text(np.ceil(np.max(refined_map["m"]))*(3/4),np.ceil(np.max(refined_map["pi"]))*(4.5/8),
        r"$\frac{(N_{\rm HPT}/N_{\rm HPT,ref})}{\sqrt{T_{\rm 4t}/T_{\rm ref}}}$", fontsize=18)
        Nlabel.set_fontsize(15)

        if show == True:   
            colbar.set_label(r"$\it η_{\rm HPT} \ \rm \left[-\right]$",fontsize = 14)
            colbar.set_ticks(np.linspace(0.5,1,11))
            colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95","1"])

        plt.plot([np.ceil(np.max(refined_map["m"]))*(3.2/4), refined_map["m"][-1][int(Num_refinement/4)]],
        [np.ceil(np.max(refined_map["pi"]))*(4.5/8), refined_map["pi"][-1][int(Num_refinement/4)]],color='k',linewidth = 0.75)

    elif type == "LPT":

        plt.xlim([np.floor(np.min(refined_map["m"])),np.ceil(np.max(refined_map["m"]))+1])
        plt.ylim([np.floor(np.min(refined_map["pi"])),np.ceil(np.max(refined_map["pi"]))-0.2])

        xlabel = plt.xlabel(r"$\frac{\it G_{\rm 45} \sqrt{T_{\rm 45t}/T_{\rm ref}}}{p_{\rm 45t}/p_{\rm ref}} \ \left[\frac{\rm kg}{\rm s}\right]$",loc='center')
        xlabel.set_fontsize(17)
        
        ylabel = plt.ylabel(r"$\pi_{\rm LPT}$",loc='center')
        ylabel.set_fontsize(17)

        plt.title('LOW PRESSURE TURBINE MAP',fontsize = 14, weight = 'bold')
        plt.scatter(m_LPT_design,pi_LPT_design,100,marker='*',color='yellow',linewidth=1.25)

        Nlabel = plt.text(np.ceil(np.max(refined_map["m"]))*(3/4),np.ceil(np.max(refined_map["pi"]))*(5.75/8),
        r"$\frac{(N_{\rm LPT}/N_{\rm LPT,ref})}{\sqrt{T_{\rm 45t}/T_{\rm ref}}}$", fontsize=18)
        Nlabel.set_fontsize(15)

        if show == True:
            colbar.set_label(r"$\it η_{\rm LPT} \ \rm \left[-\right]$",fontsize = 14)
            colbar.set_ticks(np.linspace(0.5,1,11))
            colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95","1"])

        plt.plot([np.ceil(np.max(refined_map["m"]))*(3.2/4), refined_map["m"][-1][int(Num_refinement/4)]],
        [np.ceil(np.max(refined_map["pi"]))*(5.75/8), refined_map["pi"][-1][int(Num_refinement/4)]],color='k',linewidth = 0.75)

    return plt