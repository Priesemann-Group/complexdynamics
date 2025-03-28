from argparse import ArgumentParser
import numpy as np
from scipy.integrate import odeint
import functools as functool
import general_methods as gm
import time as t

import matplotlib.pyplot as plt

def _get_arguments():
    """Return script arguments as dictionary
    Args:
        - None
    Returns:
        - dictionary containing script arguments
    """
    parser = ArgumentParser()
    parser.add_argument(
        "ExperimentOutputFolder", help="Path to Output Folder"
    )
    parser.add_argument(
        "ExperimentInputFolder", help="Path to Input-Folder whith '.yml' File"
    )
    parser.add_argument(
        "Itterator", help="Number of Job"
    )
    """ Potentially more arguments """
    return parser.parse_args()


def dSIR_ODE(p, u, t):
    omega, s, beta_0, M_max, h_thres, epsilon, nu, gamma, tau = p 
   
    S, I, H1, H = u
    R = 1 - S - I

    M = M_max - (M_max / h_thres) * epsilon * np.log(1+np.exp((h_thres - H)/epsilon))   # mitigation
    beta = beta_0 * (1-M)           # effective spreading rate
    Gamma_t = 1 + s*np.cos(omega * t)  # seasonal forcing
   
    dS = -beta * Gamma_t * I*S + nu * R
    dI =  beta * Gamma_t * I*S - gamma * I
    dH1 = (I - H1) / tau
    dH =  (H1 - H) / tau
    return [dS, dI, dH1, dH]

def winding_number(p,y_0,ntrans,nattr):
    omega = p[0]
    period = 2*np.pi / omega
    tspan = np.arange(0, ntrans*period, period/100)

    dSIR_ODE_partial = functool.partial(dSIR_ODE, p)
    #wait transient time to attractor
    y = odeint(dSIR_ODE_partial, y_0, tspan)
    # Generate TrainingData
    y_0 = y[-1]
    tspan = np.arange(0, nattr*period, period/100)
    sol = odeint(dSIR_ODE_partial, y_0, tspan)

    inf_vec = sol[:,1]
    sus_vec = sol[:,0]
    npts = sol.shape[0]

    #determine number off peaks
    npeak = 0
    for j in range(2, npts-1):
        if (inf_vec[j-1] < inf_vec[j]) and (inf_vec[j] > inf_vec[j+1]):
            npeak = npeak + 1 
    wind_no = npeak / nattr
    av_inf_no = np.mean(inf_vec)

    #Tests of individual applications
    #print(npeak, nattr, wind_no,av_inf_no)
    #plt.plot(tspan/360, sol[:,0]-1)
    #plt.plot(tspan/360, sol[:,1])
    #plt.show()

    # determine minima and maxima of S and I
    Smin = np.min(sus_vec)
    Smax = np.max(sus_vec)
    Imin = np.min(inf_vec)
    Imax = np.max(inf_vec)
    return [wind_no, av_inf_no, Smin, Smax, Imin, Imax]

def main():
    
    args = _get_arguments()
    
    # Path and file names
    #------------------------------
    PATH_IN = args.ExperimentInputFolder + '/'
    PATH_Out = args.ExperimentOutputFolder + '/'
    
    FileName_In = 'Parameter.yml'
    FileName_Out = 'InfectionData'

    Parameter = gm.yml.read(PATH_IN+FileName_In)
    
    # Integration parameter
    # -----------------------------
    Int = Parameter['Integration']
    ntrans = Int['ntrans']   # no. of transient periods (in years)
    nattr  = Int['nattr']    # no. of points on the attractor in the Poincare section (in years)
    ninit  = Int['ninit']    # no. of initial conditions (per variable)

    # Model parameters
    # -----------------------------
    Model = Parameter['Model']

    beta_0 = Model['beta_0']    # days^{-1}   base spreading rate
    gamma = Model['gamma']      # days^{-1}   recovery rate
    nu = Model['nu']            # days^{-1}   waning rate 
    tau = Model['tau']          # mitigation delay  
    M_max = Model['M_max']      # maximal mitigation
    s = Model['s']              # seasonal amplitude
    omega = Model['omega']      # days^{-1} frequency of yearly seasonal variation
    kappa = Model['kappa']      # weight factor of costs
    epsilon = Model['epsilon']  # feedback curvature
    h_thres = Model['h_thres']

    # control parameters 
    # -----------------------------
    Control = Parameter['Control']
    PatchNr = Parameter['PatchNumber']
    it = int(args.Itterator)
    #print(it)
    n_xpatch = it%PatchNr # Number of patches in x (tau) direction: itterator and n_xpatch are zero-based
    n_ypatch = it//PatchNr # Number of patches in y (M_max) direction: itterator and n_xpatch are zero-based

    tau_min = Control['tau']['min'] 
    tau_max = Control['tau']['max']
    tau_res = Control['tau']['res']
    total_tau_vec = np.arange(tau_min, tau_max, tau_res)
    
    tot_tau_len = len(total_tau_vec)

    if n_xpatch != PatchNr-1:
        tau_vec = total_tau_vec[(n_xpatch)*(tot_tau_len//PatchNr):(n_xpatch+1)*(tot_tau_len//PatchNr)]
    else:
        tau_vec = total_tau_vec[(n_xpatch)*(tot_tau_len//PatchNr):]
    n_tau = len(tau_vec)

    M_max_min = Control['M_max']['min'] 
    M_max_max = Control['M_max']['max']
    M_max_res = Control['M_max']['res']
    total_M_max_vec = np.arange(M_max_min, M_max_max, M_max_res)

    tot_M_max_len = len(total_M_max_vec)

    if n_ypatch != PatchNr-1:
        M_max_vec = total_M_max_vec[(n_ypatch)*(tot_M_max_len//PatchNr):(n_ypatch+1)*(tot_M_max_len//PatchNr)]
    else:
        M_max_vec = total_M_max_vec[(n_ypatch)*(tot_M_max_len//PatchNr):]
    n_M_max = len(M_max_vec)

    #print(tau_vec, '\n\n')
    #print(M_max_vec)

    #plt.pcolor(X,Y,patch)
    #plt.show()
    
    # Initialization of matrices
    # -------------------------------------
    wind_no_mat = np.zeros((n_M_max, n_tau, ninit*ninit))
    av_inf_no_mat = np.zeros((n_M_max, n_tau, ninit*ninit))   # average infection numbers
    Smin_mat = np.zeros((n_M_max, n_tau, ninit*ninit))
    Smax_mat = np.zeros((n_M_max, n_tau, ninit*ninit))
    Imin_mat = np.zeros((n_M_max, n_tau, ninit*ninit))
    Imax_mat = np.zeros((n_M_max, n_tau, ninit*ninit))

    total_integrations = n_tau*n_M_max*ninit*ninit
    # Calculation
    # -------------------------------------
    for i_tau,tau in enumerate(tau_vec):
        for i_M_max,M_max in enumerate(M_max_vec):
            p = [omega, s, beta_0, M_max, h_thres, epsilon, nu, gamma, tau]
            S_0_vec = np.linspace(0, 0.999, ninit)
            count_initial = 0
            for S_0 in S_0_vec:
                I_0_vec = np.linspace(0.001, 1-S_0, ninit)
                for I_0 in I_0_vec:
                    y_0 = [S_0, I_0, 0, 0]
                    wind_no, av_inf_no, Smin, Smax, Imin, Imax = winding_number(p,y_0,ntrans,nattr)
                    wind_no_mat[i_M_max, i_tau, count_initial] = wind_no
                    av_inf_no_mat[i_M_max, i_tau, count_initial] = av_inf_no
                    Smin_mat[i_M_max, i_tau, count_initial] = Smin
                    Smax_mat[i_M_max, i_tau, count_initial] = Smax
                    Imin_mat[i_M_max, i_tau, count_initial] = Imin
                    Imax_mat[i_M_max, i_tau, count_initial] = Imax
                    gm.progress_featback.printProgressBar(iteration=(i_tau*n_M_max+i_M_max)*ninit*ninit+count_initial, total=n_tau*n_M_max*ninit*ninit)
                    count_initial +=1
    
    np.savez_compressed(PATH_Out+FileName_Out+'_Patch_{}'.format(it), 
                        tau_vec=tau_vec,
                        M_max_vec=M_max_vec,
                        wind_no=wind_no_mat, 
                        av_inf_no=av_inf_no_mat, 
                        Smin=Smin_mat, Smax=Smax_mat, 
                        Imin=Imin_mat, Imax=Imax_mat)
    return total_integrations

if __name__ =='__main__':
    t_start = t.time()
    total_integrations = main()
    print('The Program run {} Minutes for (in total) {} time-integrations of the Infection.'.format((t.time()-t_start)/60, total_integrations))