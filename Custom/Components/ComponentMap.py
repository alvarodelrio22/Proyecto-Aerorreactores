import csv, numpy as np
from scipy.interpolate import RectBivariateSpline
from Components.DesignVariables import m_c_design, m_t_design, pi_c_design, pi_t_design
from scipy.optimize import newton

## COMPRESSOR: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Definition of map interpolation:

beta_data = np.linspace(0,1,13)
N_data = np.array([61,81,97,111,123])/108.5
splineC, generic_mapC = {}, {}
    
Num_refinement = 50
precision = 3

beta_refined = np.linspace(np.min(beta_data), np.max(beta_data),Num_refinement)
N_refined = np.linspace(np.min(N_data), np.max(N_data),Num_refinement)

for name in {"m", "pi", "eta"}:
    file = open("Custom/CSV Files/"+"C"+name+".csv", 'r')
    csvreader = csv.reader(file)
    data = list(csvreader)
    generic_mapC[name] = np.double(np.array(data))
    splineC[name] = RectBivariateSpline(beta_data, N_data, np.transpose(generic_mapC[name]), kx = precision, ky = precision)

# Function definition:       

def compressor(input1,input2,inputname1,inputname2,outputname):

    m_scale_factor = 1
    pi_scale_factor = 1

    # Second input name:

    if inputname2 != "beta" and inputname2 != "N":
        print("Not a valid entry")
        exit()
    elif inputname2 == "N":
        if input2 < np.min(N_data) or input2 > np.max(N_data):
            return np.NaN
    elif inputname2 == "beta":
        if input2 < np.min(beta_data) or input2 > np.max(beta_data):
            return np.NaN
        
    # Output name:

    if outputname == "m" or outputname == "eta" or outputname == "pi" or outputname == "beta" or outputname == "N":

        if outputname != "beta" and outputname != "N":

            interp = splineC[outputname]

        # Input names are beta and N:
        
        if inputname1 == "beta" and inputname2 == "N":

            if input1 < np.min(beta_data) or input1 > np.max(beta_data):
                input1 = np.NaN

            if outputname == "m":
                value = m_scale_factor*interp.ev(input1,input2)
            elif outputname == "pi":
                value = pi_scale_factor*(interp.ev(input1,input2) - 1) + 1
            elif outputname == "eta":
                value = interp.ev(input1,input2)

            return value
        
        interp2 = splineC[inputname1]

        # Case 1 for the first input name when N or beta:
        
        if inputname1 != "N" and inputname2 == "beta":

            if inputname1 == "m":

                if input1 < np.min(m_scale_factor*generic_mapC[inputname1]) or input1 > np.max(m_scale_factor*generic_mapC[inputname1]):         
                    return np.NaN

                def f(x):
                    return np.abs(m_scale_factor*interp2.ev(input2,x) - input1)
                
            elif inputname1 == "pi":

                if input1 < np.min(pi_scale_factor*(generic_mapC[inputname1] - 1) + 1) or input1 > \
                np.max(pi_scale_factor*(generic_mapC[inputname1] - 1) + 1):        
                    return np.NaN

                def f(x):
                    return np.abs(pi_scale_factor*(interp2.ev(input2,x) - 1) + 1 - input1)
                
            elif inputname1 == "eta":

                if input1 < np.min(generic_mapC[inputname1]) or input1 > np.max(generic_mapC[inputname1]):       
                    return np.NaN

                def f(x):
                    return np.abs(interp2.ev(input2,x) - input1)
            
            N_values = newton(f,N_refined,tol = 1e-6, maxiter = 100)
            N = np.NaN
            tolerance = 1e-5

            for k in range(len(N_values)):

                # Seek a value from above (k -> -k):

                if abs(N_values[-k-1]-N_values[-k]) < tolerance:
                    N = N_values[-k-1]
                    break

            if N < np.min(N_refined) or N > np.max(N_refined) or np.isnan(N):
                return np.NaN

            if outputname == "m":
                value = m_scale_factor*interp.ev(input2,N)
            elif outputname == "pi":
                value = pi_scale_factor*(interp.ev(input2,N) - 1) + 1
            elif outputname == "eta":
                value = interp.ev(input2,N)
            elif outputname == "N":
                value = N

            return value
        
        # Case 2 for the first input name when N or beta:
        
        elif inputname1 != "beta" and inputname2 == "N":

            if inputname1 == "m":

                if input1 < np.min(m_scale_factor*generic_mapC[inputname1]) or input1 > np.max(m_scale_factor*generic_mapC[inputname1]):    
                    return np.NaN
        
                def f(x):
                    return np.abs(m_scale_factor*interp2.ev(x,input2) - input1)
                
            elif inputname1 == "pi":

                if input1 < np.min(pi_scale_factor*(generic_mapC[inputname1] - 1) + 1) or input1 > \
                np.max(pi_scale_factor*(generic_mapC[inputname1] - 1) + 1):        
                    return np.NaN

                def f(x):
                    return np.abs(pi_scale_factor*(interp2.ev(x,input2) - 1) + 1 - input1)
                
            elif inputname1 == "eta":

                if input1 < np.min(generic_mapC[inputname1]) or input1 > np.max(generic_mapC[inputname1]):       
                    return np.NaN

                def f(x):
                    return np.abs(interp2.ev(x,input2) - input1)
            
            beta_values = newton(f,beta_refined,tol = 1e-6, maxiter = 100)
            beta = np.NaN
            tolerance = 1e-5

            for k in range(len(beta_values)):

                # Seek a value from above (k -> -k):

                if abs(beta_values[-k-1]-beta_values[-k]) < tolerance:
                    beta = beta_values[-k-1]
                    break

            if beta < np.min(beta_refined) or beta > np.max(beta_refined) or np.isnan(beta):
                return np.NaN                  

            if outputname == "m":
                value = m_scale_factor*interp.ev(beta,input2)
            elif outputname == "pi":
                value = pi_scale_factor*(interp.ev(beta,input2) - 1) + 1
            elif outputname == "eta":
                value = interp.ev(beta,input2)
            elif outputname == "beta":
                value = beta

            return value
        
    else:

        print("Unrecognized variable to output")
        exit()

## TURBINE: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Definition of map interpolation:

#beta_data = np.linspace(0,1,11)
#N_data = np.linspace(1/2,7/6,8)
beta_data = [0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1]
N_data = [0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2]
splineT, generic_mapT = {}, {}

Num_refinement = 50
precision = 3

beta_refined = np.linspace(np.min(beta_data), np.max(beta_data),Num_refinement)
N_refined = np.linspace(np.min(N_data), np.max(N_data),Num_refinement)

for name in {"m", "pi", "eta"}:
    file = open("Custom/CSV Files/"+"T"+name+".csv", 'r')
    csvreader = csv.reader(file)
    data = list(csvreader)
    generic_mapT[name] = np.double(np.array(data))
    splineT[name] = RectBivariateSpline(beta_data, N_data, np.transpose(generic_mapT[name]), kx = precision, ky = precision)

# Function definition:

def turbine(input1,input2,inputname1,inputname2,outputname):

    m_scale_factor = 1/64
    pi_scale_factor = 7/8

    # Second input name:

    if inputname2 != "beta" and inputname2 != "N":
        print("Not a valid entry")
        exit()
    elif inputname2 == "N":
        if input2 < np.min(N_data) or input2 > np.max(N_data):
            return np.NaN
    elif inputname2 == "beta":
        if input2 < np.min(beta_data) or input2 > np.max(beta_data):
            return np.NaN
        
    # Output name:

    if outputname == "m" or outputname == "eta" or outputname == "pi" or outputname == "beta" or outputname == "N":

        if outputname != "beta" and outputname != "N":
        
            interp = splineT[outputname]

        # Input names are beta and N:
        
        if inputname1 == "beta" and inputname2 == "N":

            if input1 < np.min(beta_data) or input1 > np.max(beta_data):
                return np.NaN

            if outputname == "m":
                value = m_scale_factor*interp.ev(input1,input2)
            elif outputname == "pi":
                value = pi_scale_factor*(interp.ev(input1,input2) - 1) + 1
            elif outputname == "eta":
                value = interp.ev(input1,input2)

            return value
        
        interp2 = splineT[inputname1]

        # Case 1 for the first input name when N or beta:
        
        if inputname1 != "N" and inputname2 == "beta":

            if inputname1 == "m":

                if input1 < np.min(m_scale_factor*generic_mapT[inputname1]) or input1 > np.max(m_scale_factor*generic_mapT[inputname1]):        
                    return np.NaN

                def f(x):
                    return np.abs(m_scale_factor*interp2.ev(input2,x) - input1)
                
            elif inputname1 == "pi":

                if input1 < np.min(pi_scale_factor*(generic_mapT[inputname1] - 1) + 1) or input1 > np.max(pi_scale_factor*(generic_mapT[inputname1] - 1) + 1):        
                    return np.NaN

                def f(x):
                    return np.abs(pi_scale_factor*(interp2.ev(input2,x)-1)+1 - input1)
                
            elif inputname1 == "eta":

                if input1 < np.min(generic_mapT[inputname1]) or input1 > np.max(generic_mapT[inputname1]):        
                    return np.NaN

                def f(x):
                    return np.abs(interp2.ev(input2,x) - input1)
            
            N_values = newton(f,N_refined,tol = 1e-6, maxiter = 100)
            N = np.NaN
            tolerance = 1e-5

            for k in range(len(N_values)):

                if abs(N_values[k-1]-N_values[k]) < tolerance:
                    N = N_values[k-1]
                    break

            if N < np.min(N_refined) or N > np.max(N_refined) or np.isnan(N):
                return np.NaN

            if outputname == "m":
                value = m_scale_factor*interp.ev(input2,N)
            elif outputname == "pi":
                value = pi_scale_factor*(interp.ev(input2,N) - 1) + 1
            elif outputname == "eta":
                value = interp.ev(input2,N)
            elif outputname == "N":
                value = N

            return value
        
        # Case 2 for the first input name when N or beta:
        
        elif inputname1 != "beta" and inputname2 == "N":

            if inputname1 == "m":

                if input1 < np.min(m_scale_factor*generic_mapT[inputname1]) or input1 > np.max(m_scale_factor*generic_mapT[inputname1]):        
                    return np.NaN

                def f(x):
                    return np.abs(m_scale_factor*interp2.ev(x,input2) - input1)
                
            elif inputname1 == "pi":

                if input1 < np.min(pi_scale_factor*(generic_mapT[inputname1] - 1) + 1) or input1 > np.max(pi_scale_factor*(generic_mapT[inputname1] - 1) + 1):        
                    return np.NaN

                def f(x):
                    return np.abs(pi_scale_factor*(interp2.ev(x,input2)-1) + 1 - input1)
                
            elif inputname1 == "eta":

                if input1 < np.min(generic_mapT[inputname1]) or input1 > np.max(generic_mapT[inputname1]):        
                    return np.NaN

                def f(x):
                    return np.abs(interp2.ev(x,input2) - input1)
            
            beta_values = newton(f,beta_refined,tol = 1e-6, maxiter = 100)
            beta = np.NaN
            tolerance = 1e-5

            for k in range(len(beta_values)):

                if abs(beta_values[k-1]-beta_values[k]) < tolerance:
                    beta = beta_values[k-1]
                    break

            if beta < np.min(beta_data) or beta > np.max(beta_data) or np.isnan(beta):
                return np.NaN                  

            if outputname == "m":
                value = m_scale_factor*interp.ev(beta,input2)
            elif outputname == "pi":
                value = pi_scale_factor*(interp.ev(beta,input2) - 1) + 1
            elif outputname == "eta":
                value = interp.ev(beta,input2)
            elif outputname == "beta":
                value = beta

            return value
        
    else:

        print("Unrecognized variable to output")
        exit()
