import numpy as np

p={
    'y0':[1-0.001, 0.001, 0, 0.001, 0.001],
    'beta': 0.5,
    'gamma':1/10,
    'nu':1/100,
    'tau': 25,
    'mmax': 0.83,
    'hthres':1/1000,
    'epsilon':0.00025,
    'a': 0.25,
    'w':2*np.pi/360,
    'tmin': 0,
    'tmax': 3000,
    'stepsize': 0.1,
    'maxstep':10,
}


cols = {
    'cstable':'#5ebb78',
    'cunstable':'#658ec1',
    'cmetastable':'#76ddc8',
    1:'#658ec1',
    2:'#76ddc8',
    3:'#76ddc8',
    4:'#5ebb78',
    'feedback_I':'#658ec1',
    'seasonality_I':'#5ebb78',
    'sirs':'royalblue',
    'hopfbifurcation':'#658ec1'
}


scenarios = {
    'mmax': {1:0.9,2:0.850,3:0.835,4:0.75},
    'tau': {1:30,2:30,3:30,4:30},
    'marker':{1:'o',2:'^',3:'*',4:'o'}
}
