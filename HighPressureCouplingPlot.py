import numpy as np, warnings, time

from Solvers.HighPressureCouplingSolver import hpCoupling
from Components.MapPlotFunction import componentPlot
from alive_progress import alive_bar

warnings.filterwarnings("ignore")

Num_points = 25
beta_HPC = np.linspace(0,1,Num_points)
N_HPC = np.linspace(0.45,1.08,Num_points)

var = ["mC", "piC", "etaC", "mT", "piT", "etaT"]
map = {}

with alive_bar(len(beta_HPC)*len(N_HPC)) as bar:

    for x in var:
        map[x] = np.zeros([Num_points, Num_points])

    for i in range(len(beta_HPC)):

        for j in range(len(N_HPC)):

            print("Iteration: " + str(i*len(beta_HPC) + (j+1)) + "/" + str(len(beta_HPC)*len(N_HPC)))

            m_25, p3t_p25t, eta_HPC, m_41, p45t_p41t, eta_HPT = hpCoupling(beta_HPC[i],N_HPC[j],100,0.9,True)

            map["mC"][i,j] = m_25
            map["piC"][i,j] = p3t_p25t
            map["etaC"][i,j] = eta_HPC

            map["mT"][i,j] = m_41
            map["piT"][i,j] = 1/p45t_p41t
            map["etaT"][i,j] = eta_HPT

            time.sleep(5e-4)
            bar()

pltC = componentPlot("HPC",False)
coupled_mapC = pltC.contourf(map["mC"],map["piC"],map["etaC"],np.linspace(0.5,0.9,100),cmap = 'jet')

pltC.title('HIGH PRESSURE COUPLING - HPC',fontsize = 14, weight = 'bold')

colbar = pltC.colorbar(coupled_mapC)
colbar.set_label(r"$\it η_{\rm HPC}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

pltC.show()

pltT = componentPlot("HPT",False)
coupled_mapT = pltT.contourf(map["mT"],map["piT"],map["etaT"],np.linspace(0.7,0.95,100),cmap = 'jet')

pltT.title('HIGH PRESSURE COUPLING - HPT',fontsize = 14, weight = 'bold')

colbar = pltT.colorbar(coupled_mapT)
colbar.set_label(r"$\it η_{\rm HPT}$",fontsize = 14)
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltT.show()

