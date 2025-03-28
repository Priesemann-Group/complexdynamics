import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors
import pandas as pd
import numpy as np
from params import cols
from scipy.stats import multivariate_normal
from scipy import signal
from pylab import cm


def set_rcParams():
    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["axes.spines.top"] = False
    mpl.rcParams["axes.titlesize"]= 8
    mpl.rcParams["xtick.labelsize"] = 8
    mpl.rcParams["ytick.labelsize"] = 8
    mpl.rcParams["axes.labelsize"] = 8
    mpl.rcParams["legend.fontsize"] = 8
    mpl.rcParams["legend.title_fontsize"] = 8
    return None


def setabc(ax,lab,x=-0.25,y=1.1):
    ax.text(x,y,f'({lab})', size=8, color='black', transform=ax.transAxes)
    return None

def get_arnold_cmap():
    windingnumbers = [1,5/4,4/3,3/2,5/3,2,5/2,3]
    windingnumberlabels = ['$1$',r'$\frac{5}{4}$',r'$\frac{4}{3}$',r'$\frac{3}{2}$',r'$\frac{5}{3}$','$2$',r'$\frac{5}{2}$','$3$','']

    cmap = cm.get_cmap('YlGnBu',len(windingnumbers)+1)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[-1] = (0, 0, 0, 1.0)
    return windingnumbers, windingnumberlabels, mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)



def plot_stability(ax,hopfcurve,mmaxvals,tauvals,alpha=1):
    
    stability = pd.read_csv(hopfcurve, sep=',', header=0)
    
    x = stability['tau']
    y = np.array(stability['mmax'])
    
   
    ax.fill_between(x,y,color=cols['cstable'],alpha=alpha)
    ax.fill_between(x,y,mmaxvals[-1]*1.1*np.ones(len(x)),color=cols['cunstable'],alpha=alpha)
    
    ax.plot(x,y, lw=1, color='darkgray')
    
    ax.set_xlim(tauvals[0],tauvals[-1])
    ax.set_ylim(mmaxvals[0],mmaxvals[-1])

    ax.set_xlabel(r'mitigation delay $\tau$')
    ax.set_ylabel('maximal mitigation $m_{max}$')
    return None
    
    
def plot_metastability(ax,summer,winter,mmaxvals,tauvals,alpha=1):
    stability_s = pd.read_csv(summer, sep=',', header=0)
    stability_w = pd.read_csv(winter, sep=',', header=0)
    
    x = np.array(stability_s['tau'])
    hopfcurve_low = np.minimum(stability_s['mmax'],stability_w['mmax'])
    hopfcurve_high = np.maximum(stability_s['mmax'],stability_w['mmax'])
      
    ax.fill_between(x,hopfcurve_low, color=cols['cstable'],alpha=alpha)
    ax.fill_between(x,hopfcurve_low,hopfcurve_high, color=cols['cmetastable'],alpha=alpha)
    ax.fill_between(x,hopfcurve_high,tauvals[-1]*np.ones(len(x)), color=cols['cunstable'],alpha=alpha)
    
    ax.plot(x,stability_s['mmax'], lw=1, color='darkgray')
    ax.plot(x,stability_w['mmax'], lw=1, color='darkgray')
    
    ax.set_xlim(tauvals[0],tauvals[-1])
    ax.set_ylim(mmaxvals[0],mmaxvals[-1])

    ax.set_xlabel(r'mitigation delay $\tau$')
    ax.set_ylabel('maximal mitigation $m_{max}$')
    return None



def plot_arnold(ax,res,cmap,param='a',tol=0.025,windingnumbers=[1,3/2,2]):

    tauline = np.sort(list(set(res['tau'])))
    sline = np.flip(np.sort(list(set(res[param]))))

    M = np.zeros([len(sline),len(tauline)])

    for i,tau in enumerate(tauline):
        subres = res[res['tau']==tau]
        for j,s in enumerate(sline):
            subsubres = subres[subres[param]==s]

            W = np.array(subsubres['W'])[0]
            for t in windingnumbers:
                if abs(W-t) <=tol:
                    M[j,i] = t
                    
    Mrev = M.copy()
    for j,t in enumerate(windingnumbers):
        Mrev[M==t] = j+1
    Mrev[M==0] = len(windingnumbers)+2
    
    show = ax.imshow(Mrev,aspect='auto', extent=(tauline.min(), tauline.max(), sline.min(), sline.max()), \
          cmap=cmap,interpolation='None')
               
    return show, Mrev, tauline, sline



def plot_bifurcation(res,ax,param,s=0.01,lw=0.4,shift=False,color='royalblue',marker='.'):
    
    line = np.sort(list(set(res[param])))
    
    for p in line:
        subres = res[res[param]==p]
        peaks = subres['i']
        
        
        ax.scatter(p*np.ones(len(peaks)),100*peaks,c=color,s=s,linewidths=lw,marker=marker)
        ax.set_ylabel('$I_k$ (%)')
        ax.set_xlabel(param)
        
    return None



def plot_perturbation(ax,res,variable='W',param='Mmax',sigma_tau=4, sigma_Mmax=0.03, relative=False, absolute=False):
    
    tauline = np.sort(list(set(res['tau'])))
    Mmaxline = np.sort(list(set(res[param])))

    W = res[variable].values.reshape((len(Mmaxline), len(tauline)))
                    
    weights = multivariate_normal.pdf(res[['tau','Mmax']], mean=[tauline.mean(),Mmaxline.mean()], cov=np.diag([sigma_tau**2, sigma_Mmax**2]))
    weights = weights.reshape((len(Mmaxline), len(tauline)))/weights.sum()
    
    N_tau_padding = max(int(len(tauline)/5), 50)
    N_Mmax_padding = max(int(len(Mmaxline)/5), 50)
    W_padded = np.pad(W, ((N_tau_padding,N_tau_padding),(N_Mmax_padding,N_Mmax_padding)), mode='edge')
    
    averages = signal.fftconvolve(W_padded, weights, mode='same')
    averages = averages[N_tau_padding: -N_tau_padding, N_Mmax_padding: -N_Mmax_padding]
    if not relative:
        difference = averages-W
    else:
        difference = (averages-W)/W
    
    if absolute:
        difference=np.abs(difference)
    
    ### Plot
    to_plot = difference
    lim = vmin=np.abs(to_plot).max()
    
    divnorm=colors.TwoSlopeNorm(vmin=np.min(to_plot[::-1,:]), vcenter=0, vmax=np.max(to_plot[::-1,:]))
    show_pert = ax.imshow(to_plot[::-1,:],aspect='auto',extent=(tauline.min(), tauline.max(), Mmaxline.min(), Mmaxline.max()),\
                      cmap='RdBu_r',norm=divnorm,interpolation='None')
    
    return show_pert

def parameter_space_histogram(ax, tau, Mmax, data, column, weight={'type':'gaussian', 'sigma_tau':4, 'sigma_Mmax':0.015}, global_range=True, bins=30, normalized=False, cumulative=True, density=True, histtype='step', color='black', lw=1, ls='-', label=None):
    data = data.dropna(subset = [column])
    if weight['type'] == 'gaussian':
        weights = multivariate_normal.pdf(data[['tau','Mmax']], mean=[tau,Mmax], cov=np.diag([weight['sigma_tau']**2, weight['sigma_Mmax']**2]))
    if weight['type'] == 'uniform':
        weights = None
        loc_tau = (data['tau'] >= tau-.5*weight['delta_tau']) & (data['tau'] < tau+.5*weight['delta_tau'])
        loc_Mmax = (data['Mmax'] >= Mmax-.5*weight['delta_Mmax']) & (data['Mmax'] < Mmax+.5*weight['delta_Mmax'])
        data = data.loc[loc_tau & loc_Mmax]

    data_to_plot = data[column]
    mean = np.average(data_to_plot, weights=weights/weights.sum())
    if normalized:
        #print(np.average(data_to_plot, weights=weights/weights.sum()))
        data_to_plot /= mean
        #print(data_to_plot.min(), data_to_plot.max())

    limits = None
    if global_range:
        limits = [data_to_plot.min(), data_to_plot.max()]

    ax.hist(data_to_plot, weights=weights, density=density, cumulative=cumulative, histtype=histtype, bins=bins, range=limits, color=color, linewidth=lw, linestyle=ls, label=label)

    return mean