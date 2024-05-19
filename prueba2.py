from Solvers.ExternalCouplingSolver import nonChokedExternalCoupling
from Miscellaneous.AuxilliaryFunctions import componentPlot
from alive_progress import alive_bar
import numpy as np

Num_refinement = 100
N = np.linspace(0.45,1.08,Num_refinement)
M0 = 0

var = ["mLPC", "piLPC", "etaLPC", "mLPT", "piLPT", "etaLPT", \
"mHPC", "piHPC", "etaHPC", "mHPT", "piHPT", "etaHPT"]
map = {}

for x in var:
    map[x] = np.empty(Num_refinement)

with alive_bar(Num_refinement) as bar:

    for i in range(Num_refinement):

        m_2, p25t_p2t, eta_LPC, m_25, p3t_p25t, eta_HPC, m_41, \
        p45t_p41t, eta_HPT, m_45, p5t_p45t, eta_LPT = nonChokedExternalCoupling(M0, N[i], 25, 0, True)

        map["mLPC"][i] = m_2
        map["piLPC"][i] = p25t_p2t
        map["etaLPC"][i] = eta_LPC

        map["mHPC"][i] = m_25
        map["piHPC"][i] = p3t_p25t
        map["etaHPC"][i] = eta_HPC

        map["mHPT"][i] = m_41
        map["piHPT"][i] = 1/p45t_p41t
        map["etaHPT"][i] = eta_HPT

        map["mLPT"][i] = m_45
        map["piLPT"][i] = 1/p5t_p45t
        map["etaLPT"][i] = eta_LPT
        
        bar()
        
for x in var:
    map[x] = map[x][~np.isnan(map[x])]

# LPC -------------------------------------------------------------------------------------------------------------------------------------------------------

pltLPC = componentPlot("LPC",False)

pltLPC.plot(map["mLPC"],map["piLPC"],color="r",linestyle="--",marker="o")
pltLPC.title('LOW PRESSURE COUPLING - LPC',fontsize = 14, weight = 'bold')
pltLPC.show()

# HPC -------------------------------------------------------------------------------------------------------------------------------------------------------

pltHPC = componentPlot("HPC",False)

pltHPC.plot(map["mHPC"],map["piHPC"],color="r",linestyle="--")
pltHPC.title('HIGH PRESSURE COUPLING - HPC',fontsize = 14, weight = 'bold')
pltHPC.show()

# HPT -------------------------------------------------------------------------------------------------------------------------------------------------------

pltHPT = componentPlot("HPT",False)

pltHPT.plot(map["mHPT"],map["piHPT"],color="r",linestyle="--")
pltHPT.title('HIGH PRESSURE COUPLING - HPT',fontsize = 14, weight = 'bold')
pltHPT.show()

# LPT -------------------------------------------------------------------------------------------------------------------------------------------------------

pltLPT = componentPlot("LPT",False)

pltLPT.plot(map["mLPT"],map["piLPT"],color="r",linestyle="--")
pltLPT.title('LOW PRESSURE COUPLING - LPT',fontsize = 14, weight = 'bold')
pltLPT.show()