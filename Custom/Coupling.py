import numpy as np, matplotlib.tri as tri, warnings

from Solvers.CouplingSolver import coupling
from Miscellaneous.AuxilliaryFunctions import componentPlot, relaxationFactor
from alive_progress import alive_bar

warnings.filterwarnings("ignore")
print(" ")

## ENGINE INTERNAL COUPLING :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Show contour or scatter: ----------------------------------------------------------------------------------------------------------------------------------

contour = False

# Points to plot: -------------------------------------------------------------------------------------------------------------------------------------------

Num_points_beta = 25
Num_points_N = 25

# Boundaries: -----------------------------------------------------------------------------------------------------------------------------------------------
 
min_beta, max_beta = 0, 1
min_N, max_N = 1/9, 10/9

# Adapt manually the initial guesses and relaxation factor to the area of the map for better efficiency: ----------------------------------------------------
# The solution may depend on initialization conditions due to potential non-injectivity of the mapping.

num_iter0 = 150

relaxation_B_lim = [max_beta, 0.90, 0.75, min_beta]
relaxation_N_lim = [min_N, 0.75, 0.95, max_N]
relaxation_matrix = \
[[0.95,0.95,0.95],  \
[0.95,0.95,0.95],   \
[0.95,0.95,0.95]]

## COUPLING -------------------------------------------------------------------------------------------------------------------------------------------------

beta_c = np.linspace(min_beta,max_beta,Num_points_beta)
N_c = np.linspace(min_N,max_N,Num_points_N)
var = ["mc", "pic", "etac", "mt", "pit", "etat"]
map = {}

for x in var:
    map[x] = np.empty([Num_points_beta, Num_points_N])

with alive_bar(Num_points_beta*Num_points_N) as bar:

    for i in range(Num_points_beta):

        for j in range(Num_points_N):

            m_2, p3t_p2t, eta_c, m_41, p5t_p41t, eta_t = coupling(beta_c[i],N_c[j],num_iter0,\
            relaxationFactor(beta_c[i],N_c[j],relaxation_B_lim,relaxation_N_lim,relaxation_matrix),True)

            map["mc"][i,j] = m_2
            map["pic"][i,j] = p3t_p2t
            map["etac"][i,j] = eta_c

            map["mt"][i,j] = m_41
            map["pit"][i,j] = 1/p5t_p41t
            map["etat"][i,j] = eta_t

            bar()

# Compressor ------------------------------------------------------------------------------------------------------------------------------------------------

    pltC = componentPlot("C",False,"",0)

    if contour:
        coupled_mapC = pltC.contour(map["mc"],map["pic"],map["etac"],np.linspace(0.5,0.9,100),cmap = 'jet')
    else:
        coupled_mapC = pltC.scatter(map["mc"],map["pic"],10,vmin=0.5,vmax=0.9,c=map["etac"],marker="x",cmap="jet")

    colbar = pltC.colorbar(coupled_mapC)
    colbar.set_label(r"$ \ \it η_{\rm C} \ [-]$",rotation = 0,fontsize = 18,horizontalalignment='left')
    colbar.set_ticks(np.linspace(0.5,0.9,9))
    colbar.set_ticklabels(["0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90"])

    pltC.show()

# Turbine ---------------------------------------------------------------------------------------------------------------------------------------------------

    pltT = componentPlot("T",False,"",0)

    if contour:
        coupled_mapT = pltT.contour(map["mt"],map["pit"],map["etat"],np.linspace(0.7,0.95,100),cmap = 'jet')
    else:
        coupled_mapT = pltT.scatter(map["mt"],map["pit"],10,vmin=0.7,vmax=0.95,c=map["etat"],marker="x",cmap="jet")

    colbar = pltT.colorbar(coupled_mapT)
    colbar.set_label(r"$ \ \it η_{\rm T} \ [-]$",rotation = 0,fontsize = 18, horizontalalignment='left')
    colbar.set_ticks(np.linspace(0.7,0.95,6))
    colbar.set_ticklabels(["0.70","0.75","0.80","0.85","0.90","0.95"])

    pltT.show()