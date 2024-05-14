import numpy as np, matplotlib.tri as tri
import warnings, time

from Solvers.LowPressureCouplingSolver import lpCoupling
from AuxilliaryFunctions.RelaxationFactor import relaxation_factor
from Components.MapPlotFunction import componentPlot
from alive_progress import alive_bar

warnings.filterwarnings("ignore")
print(" ")

# Show contour or scatter ---------------------------------------------------------------------------------------------------------------------------------

contour = False

# Points to plot: -----------------------------------------------------------------------------------------------------------------------------------------

Num_points_beta = 30
Num_points_N = 30

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
 
min_beta, max_beta = 0, 1
min_N, max_N = 0.45, 1.08

beta_LPC = np.linspace(min_beta,max_beta,Num_points_beta)
N_LPC = np.linspace(min_N,max_N,Num_points_N) 

var = ["mLPC", "piLPC", "etaLPC", "mLPT", "piLPT", "etaLPT", \
       "mHPC", "piHPC", "etaHPC", "mHPT", "piHPT", "etaHPT"]
map = {}

# Adapt manually the initial guesses and relaxation factor to the area of the map for better efficiency: ---------------------------------------------------

num_iter0 = 15

relaxation_beta_lim = [1, 0.9, 0.85, 0]
relaxation_N_lim = [0.45, 0.7, 0.95, 1.08]
relaxation_matrix = \
[[0.00,0.75,0.80],
 [0.00,0.25,0.15],
 [0.15,0.00,0.00]]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

for x in var:
    map[x] = np.empty([Num_points_beta, Num_points_N])

with alive_bar(Num_points_beta*Num_points_N) as bar:

    for i in range(Num_points_beta):

        for j in range(Num_points_N):

            m_2, p25t_p2t, eta_LPC, m_25, p3t_p25t, eta_HPC, m_41, \
            p45t_p41t, eta_HPT, m_45, p5t_p45t, eta_LPT = lpCoupling(beta_LPC[i],N_LPC[j],num_iter0,\
            relaxation_factor(beta_LPC[i],N_LPC[j],relaxation_beta_lim,relaxation_N_lim,relaxation_matrix),True)
            
            map["mLPC"][i,j] = m_2
            map["piLPC"][i,j] = p25t_p2t
            map["etaLPC"][i,j] = eta_LPC

            map["mHPC"][i,j] = m_25
            map["piHPC"][i,j] = p3t_p25t
            map["etaHPC"][i,j] = eta_HPC

            map["mHPT"][i,j] = m_41
            map["piHPT"][i,j] = 1/p45t_p41t
            map["etaHPT"][i,j] = eta_HPT

            map["mLPT"][i,j] = m_45
            map["piLPT"][i,j] = 1/p5t_p45t
            map["etaLPT"][i,j] = eta_LPT

            time.sleep(5e-5)
            bar()

# Generation of triangular mesh for data interpolation: ---------------------------------------------------------------------------------------------------

def triangulate(type,limit_skewness,limit_area):

    if type != "LPC" and type != "HPC" and type != "HPT" and type != "LPT":
        print("Invalid selection. Choose between the following components: LPC, HPC, HPT, LPT")
        return np.NaN

    # Masking off of unwanted triangles according to a limit area value: applicable to hide large elements.
    # Number of normalized standard deviations in skewness value define the rest of the mask.

    triangulation = tri.Triangulation(map["m"+type],map["pi"+type])
    triangles = triangulation.triangles

    x_tri = map["m"+type][triangles] - np.roll(map["m"+type][triangles], 1, axis=1)
    y_tri = map["pi"+type][triangles] - np.roll(map["pi"+type][triangles], 1, axis=1)

    sides = np.sort(np.sqrt(x_tri**2 + y_tri**2),axis=1)

    # Skewness is defined as proportional to the cosine of maximum angle of the triangle (always between -1 and 1/2 and it is opposite the hypotenuse):

    skewness = (1/2 - (sides[:,0]**2 + sides[:,1]**2 - sides[:,-1]**2)/(2*sides[:,0]*sides[:,1]))/(3/2)
    normalized_skewness = skewness/np.max(skewness)

    # Heron's formula to determine triangle areas:

    s = (sides[:,0]+ sides[:,1]+ sides[:,2])/2
    area = np.sqrt(s*(s - sides[:,0])*(s - sides[:,1])*(s - sides[:,2]))
    normalized_area = area/np.max(area)

    # Performing the intersection of both masks:

    mask_area = normalized_area > limit_area
    mask_skewness = normalized_skewness > limit_skewness
    mask = np.logical_or(mask_skewness,mask_area)

    triangulation.set_mask(mask)
    
    return triangulation

# Getting rid of NaN information in variables:

for x in var:
    map[x] = map[x][~np.isnan(map[x])]

# LPC -------------------------------------------------------------------------------------------------------------------------------------------------------

map["mLPC"] = map["mLPC"][~np.isnan(map["mLPC"])]
map["piLPC"] = map["piLPC"][~np.isnan(map["piLPC"])]
map["etaLPC"] = map["etaLPC"][~np.isnan(map["etaLPC"])]

pltLPC = componentPlot("LPC",False)

if contour:
    triangulation = triangulate("LPC",1,0.15)
    coupled_mapLPC = pltLPC.tricontourf(triangulation,map["etaLPC"],np.linspace(0.5,0.9,100), cmap='jet')
else:
    coupled_mapLPC = pltLPC.scatter(map["mLPC"],map["piLPC"],10,vmin=0.5,vmax=0.9,c=map["etaLPC"],marker="x",cmap="jet")

pltLPC.title('LOW PRESSURE COUPLING - LPC',fontsize = 14, weight = 'bold')

colbar = pltLPC.colorbar(coupled_mapLPC)
colbar.set_label(r"$\it η_{\rm LPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])
pltLPC.show()

# HPC -------------------------------------------------------------------------------------------------------------------------------------------------------

map["mHPC"] = map["mHPC"][~np.isnan(map["mHPC"])]
map["piHPC"] = map["piHPC"][~np.isnan(map["piHPC"])]
map["etaHPC"] = map["etaHPC"][~np.isnan(map["etaHPC"])]

pltHPC = componentPlot("HPC",False)

if contour:
    triangulation = triangulate("HPC",0.9,0.02)
    coupled_mapHPC = pltHPC.tricontourf(triangulation,map["etaHPC"],np.linspace(0.5,0.9,100), cmap='jet')
else:    
    coupled_mapHPC = pltHPC.scatter(map["mHPC"],map["piHPC"],10,vmin=0.5,vmax=0.9,c=map["etaHPC"],marker="x",cmap="jet")

pltHPC.title('HIGH PRESSURE COUPLING - HPC',fontsize = 14, weight = 'bold')

colbar = pltHPC.colorbar(coupled_mapHPC)
colbar.set_label(r"$\it η_{\rm HPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

pltHPC.show()

# HPT ------------------------------------------------------------------------------------------------------------------------------------------------------

map["mHPT"] = map["mHPT"][~np.isnan(map["mHPT"])]
map["piHPT"] = map["piHPT"][~np.isnan(map["piHPT"])]
map["etaHPT"] = map["etaHPT"][~np.isnan(map["etaHPT"])]

pltHPT = componentPlot("HPT",False)

if contour:
    triangulation = triangulate("HPT",0.8,1)
    coupled_mapHPT = pltHPT.tricontourf(triangulation,map["etaHPT"],np.linspace(0.7,0.95,100), cmap='jet')
else:
    coupled_mapHPT = pltHPT.scatter(map["mHPT"],map["piHPT"],10,vmin=0.7,vmax=0.95,c=map["etaHPT"],marker="x",cmap="jet")

pltHPT.title('HIGH PRESSURE COUPLING - HPT',fontsize = 14, weight = 'bold')

colbar = pltHPT.colorbar(coupled_mapHPT)
colbar.set_label(r"$\it η_{\rm HPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltHPT.show()

# LPT -------------------------------------------------------------------------------------------------------------------------------------------------------

map["mLPT"] = map["mLPT"][~np.isnan(map["mLPT"])]
map["piLPT"] = map["piLPT"][~np.isnan(map["piLPT"])]
map["etaLPT"] = map["etaLPT"][~np.isnan(map["etaLPT"])]

pltLPT = componentPlot("LPT",False)

if contour:
    triangulation = triangulate("LPT",1,0.002)
    coupled_mapLPT = pltLPT.tricontourf(triangulation,map["etaLPT"],np.linspace(0.7,0.95,100), cmap='jet')
else:
    coupled_mapLPT = pltLPT.scatter(map["mLPT"],map["piLPT"],10,vmin=0.7,vmax=0.95,c=map["etaLPT"],marker="x",cmap="jet")

pltLPT.title('LOW PRESSURE COUPLING - LPT',fontsize = 14, weight = 'bold')

colbar = pltLPT.colorbar(coupled_mapLPT)
colbar.set_label(r"$\it η_{\rm LPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltLPT.show()