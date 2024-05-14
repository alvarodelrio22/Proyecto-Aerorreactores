
#Define an auxilliary function for an adaptative relaxation factor:

def relaxation_factor(beta, N, beta_limits, N_limits, values):

    for i in range(len(beta_limits)-1):

        for j in range(len(N_limits)-1):

            if beta_limits[i] >= beta >= beta_limits[i+1] and N_limits[j] <= N <= N_limits[j+1]:

                relaxation_factor = values[i][j]
                return relaxation_factor