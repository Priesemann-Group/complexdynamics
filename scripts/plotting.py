import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from params import cols
from scipy.stats import multivariate_normal
from scipy import signal

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


# Sets labels to the axis 
def setabc(ax,lab,x=-0.25,y=1.1):
    ax.text(x,y,f'({lab})', size=8, color='black', transform=ax.transAxes)
    return None


def plot_stability(ax,hopfcurve,mmaxvals,tauvals,alpha=1,exponent=False):
    
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
