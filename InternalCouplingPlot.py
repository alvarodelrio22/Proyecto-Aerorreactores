from Solvers.LowPressureCouplingSolver import lpCoupling
import numpy as np, warnings
from Components.MapPlotFunction import componentPlot
import matplotlib.mlab as mlab
import scipy.interpolate as sc
from Components.ComponentMap import compressor, turbine
import concavity

warnings.filterwarnings("ignore")

Num_points = 15
 
beta_LPC = np.linspace(0,1,Num_points)
N_LPC = np.linspace(0.45,1.08,Num_points) 

var = ["mLPC", "piLPC", "etaLPC", "betaLPC", "NLPC", "mLPT", "piLPT", "etaLPT", "betaLPT", "NLPT", \
       "mHPC", "piHPC", "etaHPC", "betaHPC", "NHPC", "mHPT", "piHPT", "etaHPT", "betaHPT", "NHPT"]
map = {}

for x in var:
    map[x] = np.zeros([Num_points, Num_points])

for i in range(len(beta_LPC)):

    for j in range(len(N_LPC)):

        print("Iteration: " + str(i*len(beta_LPC) + (j+1)) + "/" + str(len(beta_LPC)*len(N_LPC)))

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

for x in var:
   map[x] = map[x][~np.isnan(map[x])]

# LPC -------------------------------------------------------------------------------------------------------------------------------------------------------

#m_enhanced = np.linspace(np.min(map["mLPC"]),np.max(map["mLPC"]),1000)
#pi_enhanced = np.linspace(np.min(map["piLPC"]),np.max(map["piLPC"]),1000)

#m_grid, pi_grid = np.meshgrid(m_enhanced,pi_enhanced)
#eta_grid = sc.griddata((map["mLPC"],map["piLPC"]),map["etaLPC"],(m_grid,pi_grid),method='cubic')

pltLPC = componentPlot("LPC",False)
#coupled_mapLPC = pltLPC.contourf(m_grid,pi_grid,np.transpose(eta_grid),np.linspace(0.5,0.9,100), cmap='jet')

# coupled_mapLPC = pltLPC.contourf(map["mLPC"], map["piLPC"], map["etaLPC"], np.linspace(0.5,0.9,100), cmap='jet')
coupled_mapLPC = pltLPC.scatter(map["mLPC"],map["piLPC"],5,c=map["etaLPC"],marker="^",cmap="jet")

pltLPC.title('LOW PRESSURE COUPLING - LPC',fontsize = 14, weight = 'bold')

colbar = pltLPC.colorbar(coupled_mapLPC)
colbar.set_label(r"$\it η_{\rm LPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

pltLPC.show()

# HPC -------------------------------------------------------------------------------------------------------------------------------------------------------

#m_enhanced = np.linspace(np.min(map["mHPC"]),np.max(map["mHPC"]),1000)
#pi_enhanced = np.linspace(np.min(map["piHPC"]),np.max(map["piHPC"]),1000)

#m_grid, pi_grid = np.meshgrid(m_enhanced,pi_enhanced)
#eta_grid = sc.griddata((map["mHPC"],map["piHPC"]),map["etaHPC"],(m_grid,pi_grid),method='cubic')

pltHPC = componentPlot("HPC",False)
#coupled_mapHPC = pltHPC.contourf(m_grid,pi_grid,np.transpose(eta_grid),np.linspace(0.5,0.9,100),cmap = 'jet')
coupled_mapHPC = pltHPC.scatter(map["mHPC"],map["piHPC"],5,c=map["etaHPC"],marker="^",cmap="jet")

pltHPC.title('HIGH PRESSURE COUPLING - HPC',fontsize = 14, weight = 'bold')

colbar = pltHPC.colorbar(coupled_mapHPC)
colbar.set_label(r"$\it η_{\rm HPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

pltHPC.show()

# HPT ------------------------------------------------------------------------------------------------------------------------------------------------------

#m_enhanced = np.linspace(np.min(map["mHPT"]),np.max(map["mHPT"]),1000)
#pi_enhanced = np.linspace(np.min(map["piHPT"]),np.max(map["piHPT"]),1000)

#m_grid, pi_grid = np.meshgrid(m_enhanced,pi_enhanced)
#eta_grid = sc.griddata((map["mHPT"],map["piHPT"]),map["etaHPT"],(m_grid,pi_grid),method='cubic')

pltHPT = componentPlot("HPT",False)
#coupled_mapHPT = pltHPT.contourf(m_grid,pi_grid,np.transpose(eta_grid),np.linspace(0.7,0.95,100),cmap = 'jet')
coupled_mapHPT = pltHPT.scatter(map["mHPT"],map["piHPT"],5,c=map["etaHPT"],marker="^",cmap="jet")

pltHPT.title('HIGH PRESSURE COUPLING - HPT',fontsize = 14, weight = 'bold')

colbar = pltHPT.colorbar(coupled_mapHPT)
colbar.set_label(r"$\it η_{\rm HPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltHPT.show()

# LPT -------------------------------------------------------------------------------------------------------------------------------------------------------

#m_enhanced = np.linspace(np.min(map["mLPT"]),np.max(map["mLPT"]),1000)
#pi_enhanced = np.linspace(np.min(map["piLPT"]),np.max(map["piLPT"]),1000)

#m_grid, pi_grid = np.meshgrid(m_enhanced,pi_enhanced)
#eta_grid = sc.griddata(map["mLPT"],map["piLPT"],map["etaLPT"],(m_grid,pi_grid),method='cubic')

pltLPT = componentPlot("LPT",False)
#coupled_mapLPT = pltLPT.contourf(m_grid,pi_grid,np.transpose(eta_grid),np.linspace(0.7,0.95,100),cmap = 'jet')
coupled_mapLPT = pltLPT.scatter(map["mLPT"],map["piLPT"],5,map["etaLPT"],marker="^",cmap="jet")

pltLPT.title('LOW PRESSURE COUPLING - LPT',fontsize = 14, weight = 'bold')

colbar = pltLPT.colorbar(coupled_mapLPT)
colbar.set_label(r"$\it η_{\rm LPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltLPT.show()