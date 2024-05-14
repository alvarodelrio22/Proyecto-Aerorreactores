import numpy as np, warnings, time

from Solvers.HighPressureCouplingSolver import hpCoupling
from Components.MapPlotFunction import componentPlot
from alive_progress import alive_bar

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

contour = True

# Coupled points: -----------------------------------------------------------------------------------------------------------------------------------------

Num_points_beta = 250
Num_points_N = 250

beta_HPC = np.linspace(0,1,Num_points_beta)
N_HPC = np.linspace(0.45,1.08,Num_points_N)

var = ["mHPC", "piHPC", "etaHPC", "mHPT", "piHPT", "etaHPT"]
map = {}

with alive_bar(Num_points_beta*Num_points_N) as bar:

    for x in var:
        map[x] = np.empty([Num_points_beta, Num_points_N])

    for i in range(Num_points_beta):

        for j in range(Num_points_N):

            print("")

            m_25, p3t_p25t, eta_HPC, m_41, p45t_p41t, eta_HPT = hpCoupling(beta_HPC[i],N_HPC[j],150,0.85,True)

            map["mHPC"][i,j] = m_25
            map["piHPC"][i,j] = p3t_p25t
            map["etaHPC"][i,j] = eta_HPC

            map["mHPT"][i,j] = m_41
            map["piHPT"][i,j] = 1/p45t_p41t
            map["etaHPT"][i,j] = eta_HPT

            time.sleep(5e-4)
            bar()

# HPC -----------------------------------------------------------------------------------------------------------------------------------------------------

pltC = componentPlot("HPC",False)

if contour:
    coupled_mapC = pltC.contourf(map["mHPC"],map["piHPC"],map["etaHPC"],np.linspace(0.5,0.9,100),cmap = 'jet')
else:
    coupled_mapC = pltC.scatter(map["mHPC"],map["piHPC"],10,vmin=0.5,vmax=0.9,c=map["etaHPC"],marker="x",cmap="jet")

pltC.title('HIGH PRESSURE COUPLING - HPC',fontsize = 14, weight = 'bold')

colbar = pltC.colorbar(coupled_mapC)
colbar.set_label(r"$\it η_{\rm HPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

pltC.show()

# HPT -----------------------------------------------------------------------------------------------------------------------------------------------------

pltT = componentPlot("HPT",False)

if contour:
    coupled_mapT = pltT.contourf(map["mHPT"],map["piHPT"],map["etaHPT"],np.linspace(0.7,0.95,100),cmap = 'jet')
else:
    coupled_mapT = pltT.scatter(map["mHPT"],map["piHPT"],10,vmin=0.7,vmax=0.95,c=map["etaHPT"],marker="x",cmap="jet")

pltT.title('HIGH PRESSURE COUPLING - HPT',fontsize = 14, weight = 'bold')

colbar = pltT.colorbar(coupled_mapT)
colbar.set_label(r"$\it η_{\rm HPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltT.show()

