from Solvers.LowPressureCouplingSolver import lpCoupling
import numpy as np, warnings
from Components.MapPlotFunction import componentPlot

warnings.filterwarnings("ignore")

Num_points = 50
 
beta_LPC = np.linspace(0,1,Num_points)
N_LPC = np.linspace(0.45,1.08,Num_points) 

var = ["mLPC", "piLPC", "etaLPC", "mLPT", "piLPT", "etaLPT", "mHPC", "piHPC", "etaHPC", "mHPT", "piHPT", "etaHPT"]
map = {}

for x in var:
    map[x] = np.zeros([Num_points, Num_points])

for i in range(len(beta_LPC)):

    for j in range(len(N_LPC)):

        print("Iteration: " + str(i*len(beta_LPC) + (j+1)) + "/" + str(len(beta_LPC)*len(N_LPC)))

        m_2, m_25, p25t_p2t, p3t_p25t, eta_LPC, eta_HPC, m_41, m_45, p45t_p41t, p5t_p45t, eta_LPT, eta_HPT = lpCoupling(beta_LPC[i],N_LPC[j],200,0.25,True)

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
      
warnings.filterwarnings("default")

pltLPC = componentPlot("LPC",False)
coupled_mapLPC = pltLPC.contourf(map["mLPC"],map["piLPC"],map["etaLPC"],np.linspace(0.5,0.9,100),cmap = 'jet')

pltLPC.title('LOW PRESSURE COUPLING - LPC',fontsize = 14, weight = 'bold')

colbar = pltLPC.colorbar(coupled_mapLPC)
colbar.set_label(r"$\it η_{\rm LPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

pltLPC.show()

pltLPT = componentPlot("LPT",False)
coupled_mapT = pltLPT.contourf(map["mLPT"],map["piLPT"],map["etaLPT"],np.linspace(0.5,0.95,100),cmap = 'jet')

pltLPT.title('LOW PRESSURE COUPLING - LPT',fontsize = 14, weight = 'bold')

colbar = pltLPT.colorbar(coupled_mapT)
colbar.set_label(r"$\it η_{\rm LPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltLPT.show()

pltHPC = componentPlot("HPC",False)
coupled_mapHPC = pltHPC.contourf(map["mHPC"],map["piHPC"],map["etaHPC"],np.linspace(0.5,0.9,100),cmap = 'jet')

pltHPC.title('HIGH PRESSURE COUPLING - HPC',fontsize = 14, weight = 'bold')

colbar = pltHPC.colorbar(coupled_mapHPC)
colbar.set_label(r"$\it η_{\rm HPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

pltHPC.show()

print(map["mLPT"])
print(map["piLPT"])
print(map["etaLPT"])

pltHPT = componentPlot("HPT",False)
coupled_mapHPT = pltHPT.contourf(map["mHPT"],map["piHPT"],map["etaHPT"],np.linspace(0.5,1,100),cmap = 'jet')

pltHPT.title('HIGH PRESSURE COUPLING - HPT',fontsize = 14, weight = 'bold')

colbar = pltHPT.colorbar(coupled_mapHPT)
colbar.set_label(r"$\it η_{\rm HPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltHPT.show()