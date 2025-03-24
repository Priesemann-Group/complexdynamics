from scipy.linalg import eigvals
from scipy.optimize import root
import numpy as np
import csv




def softplus(H,mmax,hthres,epsilon):
    slope = mmax/hthres
    return slope*epsilon*np.log(np.exp(1/epsilon*(hthres-H))+1)+(1-mmax)

def softplus_deriv(H,mmax,hthres,epsilon):
    slope = mmax/hthres
    return -slope*np.exp(1/epsilon*(hthres-H))/(1+np.exp(1/epsilon*(hthres-H)))

def FP_I(I, *args):
    b = args[0]
    o = args[1]
    g = args[2]
    hthres = args[3]
    mmax = args[4]
    epsilon = args[5]
    return b*I+o/softplus(I,mmax,hthres,epsilon)-b*o/g*(1-I)

def FPparams(p,iguess=0.001):
    hthres = p['hthres']
    mmax = p['mmax']
    epsilon = p['epsilon']
    b = p['beta']
    o = p['nu']
    g = p['gamma']
    
    sol = root(FP_I,iguess, args=(b,o,g,hthres,mmax,epsilon))

    ISTAR = sol.x[0]
    RSTAR = g/o*ISTAR
    SSTAR = g/b*1/softplus(ISTAR,mmax,hthres,epsilon)
    
    return SSTAR,ISTAR,RSTAR

def jacobian(p):
    hthres = p['hthres']
    mmax = p['mmax']
    epsilon = p['epsilon']
    b = p['beta']
    o = p['nu']
    g = p['gamma']
    tau = p['tau']
    SSTAR,ISTAR,RSTAR = FPparams(p)
    spstar = softplus(ISTAR,mmax,hthres,epsilon)
    spderivstar = softplus_deriv(ISTAR,mmax,hthres,epsilon)

    row1 = [-b*spstar*ISTAR-o, -b*spstar*SSTAR-o,0,-b*spderivstar*SSTAR*ISTAR]
    row2 = [b*spstar*ISTAR, b*spstar*SSTAR-g, 0, b*spderivstar*SSTAR*ISTAR]
    row3 = [0,1/tau,-1/tau,0]
    row4 = [0,0,1/tau,-1/tau]
    J = [row1,row2,row3,row4]
    return J

def largestEW(p):
    J = jacobian(p)
    EW = eigvals(J)
    return np.max(EW.real)


def calc_hopfcurve(path, p, beta, mmaxline=np.linspace(0.5,1,100), tauline=np.linspace(1,35,100)):
    p['beta'] = beta 
    
    with open(path, 'w+', newline='') as file:

        writer = csv.writer(file)
        writer.writerow(('tau','mmax'))

        hopfpoint = 0
        for tau in tauline:
            p['tau'] = tau
            EW_FIX=[]
            for mmax in mmaxline:
                p['mmax'] = mmax
                EW = largestEW(p)
                EW_FIX.append(EW)
            if np.min(EW_FIX) >=0:
                hopfpoint = 0
            elif np.max(EW_FIX) <=0:
                hopfpoint = 1.02
            else:
                hopfpoint = mmaxline[np.argmin(np.absolute(EW_FIX))]
            writer.writerow((tau,hopfpoint))
    return None


def calc_stabilitymatrix(p,beta,mmaxline=np.linspace(0,1,100), tauline=np.linspace(1,50,100)):
    p['beta']=beta

    M = np.zeros([len(tauline),len(mmaxline)])

    for i,mmax in enumerate(mmaxline):
        p['mmax'] = mmax
        for j,tau in enumerate(tauline):
            p['tau'] = tau
            EW = largestEW(p)
            M[j,i] = np.sign(EW)

    return M

