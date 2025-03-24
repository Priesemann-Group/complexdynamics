import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

class model:
    def __init__(
        self,
        y0, beta, gamma, nu, tau, mmax, hthres, epsilon, a, w,
        tmin, tmax, stepsize, maxstep
    ):
        
        self.y0 = y0
        self.beta = beta
        self.gamma = gamma
        self.nu = nu
        self.tau = tau 
        
        self.mmax = mmax
        self.hthres = hthres
        self.epsilon = epsilon
        self.a = a
        self.w = w
        
        self.tmin = tmin
        self.tmax = tmax
        self.stepsize = stepsize
        self.maxstep = maxstep
        
    
    def FOI_softplus(self,H,I):
        return self.softplus(H)*I*self.beta
    
    def softplus(self,H):
        slope = self.mmax/self.hthres
        return slope*self.epsilon*np.log(np.exp(1/self.epsilon*(self.hthres-H))+1)+(1-self.mmax)
    
    def seasonalforcing(self,t):
        return (1+self.a*np.sin(self.w*t))
    
    def i_peaks(self,t,y):
        FOI = self.FOI_softplus(y[4],y[1])
        return FOI*y[0]*self.seasonalforcing(t) - self.gamma*y[1]
       
    def fun(self,t,y):
        S,I,R,H1,H = y

        FOI = self.FOI_softplus(H,I)
        
        dS = -FOI*S*self.seasonalforcing(t) + self.nu*R
        dI = FOI*S*self.seasonalforcing(t) - self.gamma*I
        dR = self.gamma*I - self.nu*R
        dH1 = 1/self.tau*(I-H1)
        dH = 1/self.tau*(H1-H)
        return [dS,dI,dR,dH1,dH]

    def run(self):
        event1 = lambda t,x: self.i_peaks(t,x)
        event1.direction = -1
            
       
        toutput = np.arange(self.tmin, self.tmax, self.stepsize)
        res = solve_ivp(self.fun, (self.tmin,self.tmax), self.y0, t_eval=toutput, max_step = self.maxstep,events=[event1])
        self.times = res['t']
        self.data = res['y']
        self.events = res['y_events']
        self.t_events = res['t_events']
        return self.times, self.data, self.events, self.t_events
    
    
