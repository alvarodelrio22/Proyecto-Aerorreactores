import numpy as np, warnings, time

from Solvers.InternalCouplingSolver import hpCoupling
from Miscellaneous.AuxilliaryFunctions import componentPlot, relaxationFactor
from alive_progress import alive_bar

warnings.filterwarnings("ignore")
print(" ")

# Show contour or scatter: ----------------------------------------------------------------------------------------------------------------------------------

contour = True

# Points to plot: -------------------------------------------------------------------------------------------------------------------------------------------

Num_points_beta = 10
Num_points_N = 10

# Boundaries: -----------------------------------------------------------------------------------------------------------------------------------------------
 
min_beta, max_beta = 0, 1
min_N, max_N = 0.45, 1.08

# Adapt manually the initial guesses and relaxation factor to the area of the map for better efficiency: ----------------------------------------------------

num_iter0 = 100

relaxation_B_lim = [max_beta, 0.90, 0.75, min_beta]
relaxation_N_lim = [min_N, 0.75, 0.95, max_N]
relaxation_matrix = \
[[0.65,0.75,0.85],  \
[0.65,0.75,0.90],   \
[0.70,0.85,0.90]]

# -----------------------------------------------------------------------------------------------------------------------------------------------------------

beta_HPC = np.linspace(min_beta,max_beta,Num_points_beta)
N_HPC = np.linspace(min_N,max_N,Num_points_N)
var = ["mHPC", "piHPC", "etaHPC", "mHPT", "piHPT", "etaHPT"]
map = {}

for x in var:
    map[x] = np.empty([Num_points_beta, Num_points_N])

with alive_bar(Num_points_beta*Num_points_N) as bar:

    for i in range(Num_points_beta):

        for j in range(Num_points_N):

            m_25, p3t_p25t, eta_HPC, m_41, p45t_p41t, eta_HPT = hpCoupling(beta_HPC[i],N_HPC[j],num_iter0,\
            relaxationFactor(beta_HPC[i],N_HPC[j],relaxation_B_lim,relaxation_N_lim,relaxation_matrix),True)

            map["mHPC"][i,j] = m_25
            map["piHPC"][i,j] = p3t_p25t
            map["etaHPC"][i,j] = eta_HPC

            map["mHPT"][i,j] = m_41
            map["piHPT"][i,j] = 1/p45t_p41t
            map["etaHPT"][i,j] = eta_HPT

            time.sleep(5e-4)
            bar()

# HPC -------------------------------------------------------------------------------------------------------------------------------------------------------

pltC = componentPlot("HPC",False,"",0)

if contour:
    coupled_mapC = pltC.contourf(map["mHPC"],map["piHPC"],map["etaHPC"],np.linspace(0.5,0.9,100),cmap = 'jet')
else:
    coupled_mapC = pltC.scatter(map["mHPC"],map["piHPC"],10,vmin=0.5,vmax=0.9,c=map["etaHPC"],marker="x",cmap="jet")

colbar = pltC.colorbar(coupled_mapC)
colbar.set_label(r"$ \ \it η_{\rm HPC} \ [-]$",rotation=0,fontsize = 18,horizontalalignment='left')
colbar.set_ticks(np.linspace(0.5,0.9,9))
colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

pltC.show()

# HPT -------------------------------------------------------------------------------------------------------------------------------------------------------

pltT = componentPlot("HPT",False,"",0)

if contour:
    coupled_mapT = pltT.contourf(map["mHPT"],map["piHPT"],map["etaHPT"],np.linspace(0.7,0.95,100),cmap = 'jet')
else:
    coupled_mapT = pltT.scatter(map["mHPT"],map["piHPT"],10,vmin=0.7,vmax=0.95,c=map["etaHPT"],marker="x",cmap="jet")

colbar = pltT.colorbar(coupled_mapT)
colbar.set_label(r"$ \ \it η_{\rm HPC} \ [-]$",rotation=0,fontsize = 18,horizontalalignment='left')
colbar.set_ticks(np.linspace(0.7,0.95,6))
colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

pltT.show()

