import numpy as np, matplotlib.tri as tri
import warnings, time

from Solvers.LowPressureCouplingSolver import lpCoupling
from Components.MapPlotFunction import componentPlot
from alive_progress import alive_bar

warnings.filterwarnings("ignore")
print(" ")

Num_points = 25
 
beta_LPC = np.linspace(0,1,Num_points)
N_LPC = np.linspace(0.45,1.08,Num_points) 

var = ["mLPC", "piLPC", "etaLPC", "betaLPC", "NLPC", "mLPT", "piLPT", "etaLPT", "betaLPT", "NLPT", \
       "mHPC", "piHPC", "etaHPC", "betaHPC", "NHPC", "mHPT", "piHPT", "etaHPT", "betaHPT", "NHPT"]
map = {}

for x in var:
    map[x] = np.zeros([Num_points, Num_points])

with alive_bar(len(beta_LPC)*len(N_LPC)) as bar:

    for i in range(len(beta_LPC)):

        for j in range(len(N_LPC)):

            m_2, p25t_p2t, eta_LPC, m_25, p3t_p25t, eta_HPC, beta_HPC, N_HPC, m_41, \
            p45t_p41t, eta_HPT, beta_HPT, N_HPT, m_45, p5t_p45t, eta_LPT, beta_LPT, N_LPT = lpCoupling(beta_LPC[i],N_LPC[j],50,0.65,True)
            
            map["mLPC"][i,j] = m_2
            map["piLPC"][i,j] = p25t_p2t
            map["etaLPC"][i,j] = eta_LPC
            map["betaLPC"][i,j] = beta_LPC[i]
            map["NLPC"][i,j] = N_LPC[j]

            map["mHPC"][i,j] = m_25
            map["piHPC"][i,j] = p3t_p25t
            map["etaHPC"][i,j] = eta_HPC
            map["betaHPC"][i,j] = beta_HPC
            map["NHPC"][i,j] = N_HPC

            map["mHPT"][i,j] = m_41
            map["piHPT"][i,j] = 1/p45t_p41t
            map["etaHPT"][i,j] = eta_HPT
            map["betaHPT"][i,j] = beta_HPT
            map["NHPT"][i,j] = N_HPT

            map["mLPT"][i,j] = m_45
            map["piLPT"][i,j] = 1/p5t_p45t
            map["etaLPT"][i,j] = eta_LPT
            map["betaLPT"][i,j] = beta_LPT
            map["NLPT"][i,j] = N_LPT

            time.sleep(5e-4)
            bar()

# Getting rid of NaN information in variables:

for x in var:
   map[x] = map[x][~np.isnan(map[x])]

def triangulate(type,std_limit_sides,std_limit_area):

    if type != "LPC" and type != "HPC" and type != "HPT" and type != "LPT":
        print("Invalid selection. Choose between the following components: LPC, HPC, HPT, LPT")
        return np.NaN

    # Masking off of unwanted triangles according to a limit area value: applicable to hide large elements.
    # Number of normalized standard deviations in side value define the rest of the mask to avoid skewness.

    triangulation = tri.Triangulation(map["m"+type],map["pi"+type])
    triangles = triangulation.triangles

    x_tri = map["m"+type][triangles] - np.roll(map["m"+type][triangles], 1, axis=1)
    y_tri = map["pi"+type][triangles] - np.roll(map["pi"+type][triangles], 1, axis=1)

    # Heron's formula to determine triangle areas:

    sides = np.sqrt(x_tri**2 + y_tri**2)
    std_sides = np.std(sides.flatten())
    mean_sides = np.mean(sides.flatten())

    s = (sides[:,0]+ sides[:,1]+ sides[:,2])/2
    areas = np.sqrt(s*(s - sides[:,0])*(s - sides[:,1])*(s - sides[:,2]))
    std_areas = np.std(areas)
    mean_areas = np.mean(areas)

    normalized_sides = np.max(np.abs((sides - mean_sides)/std_sides),axis=1)
    normalized_area = np.abs((areas - mean_areas)/std_areas)

    # Perform the intersection of both masks:

    mask_sides = normalized_sides > std_limit_sides
    mask_area = normalized_area > std_limit_area
    mask = np.logical_and(mask_sides,mask_area)

    triangulation.set_mask(mask)
    
    return triangulation

# LPC -------------------------------------------------------------------------------------------------------------------------------------------------------

triangulation = triangulate("LPC",1.5,1.5)

pltLPC = componentPlot("LPC",False)
coupled_mapLPC = pltLPC.tricontourf(triangulation,map["etaLPC"],np.linspace(0.5,0.9,100), cmap='jet')
pltLPC.triplot(triangulation, 'b-', linewidth = 1)
pltLPC.show()

# coupled_mapLPC = pltLPC.contourf(map["mLPC"], map["piLPC"], map["etaLPC"], np.linspace(0.5,0.9,100), cmap='jet')
# coupled_mapLPC = pltLPC.scatter(map["mLPC"],map["piLPC"],5,c=map["etaLPC"],marker="^",cmap="jet")

pltLPC.title('LOW PRESSURE COUPLING - LPC',fontsize = 14, weight = 'bold')

colbar = pltLPC.colorbar()
colbar.set_label(r"$\it η_{\rm LPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

# HPC -------------------------------------------------------------------------------------------------------------------------------------------------------

triangulation = triangulate("HPC",1.5,1.5)

pltHPC = componentPlot("HPC",False)
coupled_mapHPC = pltHPC.tricontourf(triangulation,map["etaHPC"],np.linspace(0.5,0.9,100), cmap='jet')
pltHPC.triplot(triangulation, 'b-', linewidth = 1)
pltHPC.show()

#coupled_mapHPC = pltHPC.contourf(m_grid,pi_grid,np.transpose(eta_grid),np.linspace(0.5,0.9,100),cmap = 'jet')
#coupled_mapHPC = pltHPC.scatter(map["mHPC"],map["piHPC"],5,c=map["etaHPC"],marker="^",cmap="jet")

pltHPC.title('HIGH PRESSURE COUPLING - HPC',fontsize = 14, weight = 'bold')

colbar = pltHPC.colorbar(coupled_mapHPC)
colbar.set_label(r"$\it η_{\rm HPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

pltHPC.show()

# HPT ------------------------------------------------------------------------------------------------------------------------------------------------------

triangulation = triangulate("HPT",1.5,1.5)

pltHPT = componentPlot("HPT",False)
coupled_mapHPT = pltHPT.tricontourf(triangulation,map["etaHPT"],np.linspace(0.5,0.9,100), cmap='jet')
pltHPT.triplot(triangulation, 'b-', linewidth = 1)
pltHPT.show()

#coupled_mapHPT = pltHPT.contourf(m_grid,pi_grid,np.transpose(eta_grid),np.linspace(0.7,0.95,100),cmap = 'jet')
#coupled_mapHPT = pltHPT.scatter(map["mHPT"],map["piHPT"],5,c=map["etaHPT"],marker="^",cmap="jet")

pltHPT.title('HIGH PRESSURE COUPLING - HPT',fontsize = 14, weight = 'bold')

colbar = pltHPT.colorbar(coupled_mapHPT)
colbar.set_label(r"$\it η_{\rm HPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltHPT.show()

# LPT -------------------------------------------------------------------------------------------------------------------------------------------------------

triangulation = triangulate("LPT",1.5,1.5)

pltLPT = componentPlot("LPT",False)
coupled_mapLPT = pltLPT.tricontourf(triangulation,map["etaLPT"],np.linspace(0.5,0.9,100), cmap='jet')
pltLPT.triplot(triangulation, 'b-', linewidth = 1)
pltLPT.show()

#coupled_mapLPT = pltLPT.contourf(m_grid,pi_grid,np.transpose(eta_grid),np.linspace(0.7,0.95,100),cmap = 'jet')
#coupled_mapLPT = pltLPT.scatter(map["mLPT"],map["piLPT"],5,map["etaLPT"],marker="^",cmap="jet")

pltLPT.title('LOW PRESSURE COUPLING - LPT',fontsize = 14, weight = 'bold')

colbar = pltLPT.colorbar(coupled_mapLPT)
colbar.set_label(r"$\it η_{\rm LPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltLPT.show()