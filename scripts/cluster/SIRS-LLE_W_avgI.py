from argparse import ArgumentParser
import numpy as np
import functools as functool
import general_methods as gm
import time as t

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

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
    parser.add_argument(
        "ProgressFeatback", help="bool (0,1) whether featback of progess is returned."
    )
    """ Potentially more arguments """
    return parser.parse_args()


def dSIR_ODE(p, t, u):
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

def dSIR_withJacoby_ODE(p, t, u):
    omega, s, beta_0, M_max, h_thres, epsilon, nu, gamma, tau = p 
   
    S, I, H1, H = u[:4]
    R = 1 - S - I

    M = M_max - (M_max / h_thres) * epsilon * np.log(1+np.exp((h_thres - H)/epsilon))   # mitigation
    beta = beta_0 * (1-M)           # effective spreading rate
    
    dbeta = - beta_0 * (M_max / h_thres * np.exp(1 / epsilon * (h_thres - H))) / (1 + np.exp(1/epsilon*(h_thres-H)))
    
    Gamma_t = 1 + s*np.cos(omega * t)  # seasonal forcing
   
    dS = -beta * Gamma_t * I * S + nu * R
    dI =  beta * Gamma_t * I * S - gamma * I
    dH1 = (I - H1) / tau
    dH =  (H1 - H) / tau

    ddS = (-beta*Gamma_t*I-nu) * u[4] + (-beta*Gamma_t*S-nu) * u[5] + (-dbeta*S*I) * Gamma_t * u[7]
    ddI = (beta*Gamma_t*I) * u[4] + (beta*Gamma_t*S-gamma) * u[5] + (dbeta*S*I) * Gamma_t * u[7]
    ddH1 = 1 / tau*u[5] - 1 / tau*u[6]
    ddH = 1 / tau*u[6] - 1 / tau*u[7]
    
    return [dS, dI, dH1, dH, ddS, ddI, ddH1, ddH]

def run_model(p, y_0, ntrans, nattr, d_0, nitter):
    #definitions of system (w.r.t p)
    dSIR_ODE_partial = functool.partial(dSIR_ODE, p)  # creates the function with p pluged in
    dSIR_withJacoby_ODE_partial = functool.partial(dSIR_withJacoby_ODE, p)
    lamb=0

    period = 2*np.pi / p[0]
    tend_trans = ntrans*period
    tend_span = nattr*period
    ttrans = np.arange(0, ntrans*period, period/100)
    tspan = np.arange(0, nattr*period, period/100)
    
    #wait transient time to attractor
    y = solve_ivp(dSIR_ODE_partial, (0, tend_trans), y_0, t_eval=ttrans, method='DOP853', rtol=1e-7, atol=1e-10)['y'].T

    #y = odeint(dSIR_ODE_partial, y_0, ttrans)
    y_01 = y[-1]

    #adjusting perturbation direction
    dir = np.array([-1, 1, 0, 0])
    delta_0 = d_0 * dir / np.linalg.norm(dir)
    u_0 = np.append(y_01, delta_0)
    #u = odeint(dSIR_withJacoby_ODE_partial, u_0, ttrans)
    u = solve_ivp(dSIR_withJacoby_ODE_partial, (0, tend_trans), u_0, t_eval=ttrans, method='DOP853', rtol=1e-7, atol=1e-10)['y'].T
    d_end = np.linalg.norm(u[-1, 4:])
    if d_end == 0:
        delta_0 = d_0 * dir / np.linalg.norm(dir)
    else:
        delta_0 = d_0/d_end*u[-1, 4:]
    u_0 = np.append(u[-1, :4], delta_0)

    #calculation of Lyapunovexponent
    npeak = 0
    navgI = 0
    for i in range(nitter):
        #integration
        sol = solve_ivp(dSIR_withJacoby_ODE_partial, (0, tend_span), u_0, t_eval=tspan, method='DOP853', rtol=1e-7, atol=1e-10)['y'].T
        #need to use DOP853 method in solve_ivp for accurate computation 
        
        # Winding Number, average infections:
        inf_vec = sol[:,1]
        npts = sol.shape[0]
        
        for j in range(2, npts-1):
            if (inf_vec[j-1] < inf_vec[j]) and (inf_vec[j] > inf_vec[j+1]):
                npeak = npeak + 1 
        
        navgI = navgI + np.mean(inf_vec)

        #calculation of norm and rescaling of initial perturbation
        d_end = np.linalg.norm(sol[-1, 4:])
        delta_0 = d_0/d_end*sol[-1, 4:]
        u_0 = np.append(sol[-1, :4], delta_0)
        
        #adding to lyapunovexponent (needs subsequent rescaling)
        lamb += np.log(d_end/d_0)

    wind_no = npeak / (nitter*tspan[-1])
    mean_infections = navgI / nitter 
    LLE = lamb/(nitter*tspan[-1])
    return LLE, wind_no, mean_infections 

def main():
    
    args = _get_arguments()
    progress_featback = int(args.ProgressFeatback)
    
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
    l_itterations = Int['nitter']
    d_0 = Int['delta']
    
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

    # Initialization of matrices
    # -------------------------------------
    LLE = np.zeros((n_M_max, n_tau))
    W = np.zeros((n_M_max, n_tau))
    Imean = np.zeros((n_M_max, n_tau))
    total_sets = n_tau*n_M_max

    # Calculation
    # -------------------------------------
    for i_tau,tau in enumerate(tau_vec):
        for i_M_max,M_max in enumerate(M_max_vec):
            p = [omega, s, beta_0, M_max, h_thres, epsilon, nu, gamma, tau]
            
            y_0 = [0.99, 0.01, 0, 0] #S, I, H, H1

            lamb, winding, infections = run_model(p, y_0, ntrans, nattr, d_0, l_itterations)
            LLE[i_M_max, i_tau] = lamb
            W[i_M_max, i_tau] = winding
            Imean[i_M_max, i_tau] = infections
            
            if progress_featback:
                gm.progress_featback.printProgressBar(iteration=(i_tau*n_M_max+i_M_max), total=total_sets)
    
    np.savez_compressed(PATH_Out+FileName_Out+'_Patch_{}'.format(it), 
                        tau_vec=tau_vec,
                        M_max_vec=M_max_vec,
                        LLE=LLE,
                        Winding=W,
                        Imean=Imean)
    # Here need to save the others as well                     
    return total_sets

if __name__ =='__main__':
    t_start = t.time()
    total_integrations = main()
    print('The Program run {} Minutes for (in total) {} time-integrations of the Infection.'.format((t.time()-t_start)/60, total_integrations))