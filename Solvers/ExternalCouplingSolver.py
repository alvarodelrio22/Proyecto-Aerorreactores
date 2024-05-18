
from Solvers.InternalCouplingSolver import lpCoupling

from AuxilliaryFunctions.RelaxationFactor import relaxation_factor as relaxation_factor_HPC
from Components.DesignVariables import gamma_c, gamma_e, Cp_c, Cp_e, N_ref_LPC, N_ref_LPT, b_25, b_3, eta_mLP, f_assumed
import numpy as np