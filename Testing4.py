#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 08:16:07 2021

@author: yisinuo
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as plp
from matplotlib.patches import ConnectionPatch
from scipy.optimize import fmin_tnc 

from raytracer import *

#%%
# Set variables
ar = np.sqrt(50 * 50 - 45 * 45)
convex_1 = ReverseSphericalRefraction(25, 1.5168, 1, ar, 0.02)
plane_1 = ReverseSphericalRefraction(20, 1, 1.5168, ar, 0)
f1 = focus(plane_1, convex_1) 
screen_1 = OutputPlane(f1)
convex = SphericalRefraction(20, 1, 1.5168, ar, 0.02)
plane = SphericalRefraction(25, 1.5168, 1, ar, 0)
f2 = focus(convex, plane) 
screen = OutputPlane(f2) 
#------------------------------------------------------------------------
# LENS OPTIMIZATION
def RMS2(c1, c2):
    beam = bundle(6, 10)
    centralray = Ray([0, 0, 0],[0, 0, 1])
    xf, yf = beam.xy(centralray)
    lens1 = OptimizeLens(20, 1, 1.5168, 1e100, c1)
    lens2 = ReverseOptimizeLens(25, 1.5168, 1, 1e100, c2)
    screen = OutputPlane(f1)
    dlist = []
    for p in np.arange(0, len(xf)): 
        r = Ray([xf[p], yf[p], 0], [0, 0, 1])
        lens1.propagate_ray(r)
        lens2.propagate_ray(r)
        screen.propagate_ray(r)  
        x = r.p()[0]
        y = r.p()[1]
        dist_sq = x * x + y * y
        dlist.append(dist_sq)
    R = np.sqrt(sum(dlist) / len(dlist))
    return R 


import scipy.optimize
result = scipy.optimize.minimize(RMS2, 0.1, 0, bounds=[(0,1)])   
print(result)


#%%
foc = f1
def RMS_optimize(params):
    c1, c2= params
    beam = bundle(6, 10)
    centralray = Ray([0, 0, 0],[0, 0, 1])
    xf, yf = beam.xy(centralray)
    lens1 = OptimizeLens(20, 1, 1.5168, 1e100, c1)
    lens2 = ReverseOptimizeLens(25, 1.5168, 1, 1e100, c2)
    screen = OutputPlane(foc)
    dlist = []
    for p in np.arange(0, len(xf)): 
        r = Ray([xf[p], yf[p], 0], [0, 0, 1])
        lens1.propagate_ray(r)
        lens2.propagate_ray(r)
        screen.propagate_ray(r)  
        x = r.p()[0]
        y = r.p()[1]
        dist_sq = x * x + y * y
        dlist.append(dist_sq)
    R = np.sqrt(sum(dlist) / len(dlist))
    return R 


# In the case of fixed output plane at f1 = 121.7 mm, plane first orientation
guess = np.array([0.016, 0.0026])

size1 = scipy.optimize.minimize(RMS_optimize, 
                                x0 = guess, bounds=[(0,1),(0,1)]) 
print(size1)
#%%
foc = f2
def RMS_optimize(params):
    c1, c2= params
    beam = bundle(6, 10)
    centralray = Ray([0, 0, 0],[0, 0, 1])
    xf, yf = beam.xy(centralray)
    lens1 = OptimizeLens(20, 1, 1.5168, 1e100, c1)
    lens2 = ReverseOptimizeLens(25, 1.5168, 1, 1e100, c2)
    screen = OutputPlane(foc)
    dlist = []
    for p in np.arange(0, len(xf)): 
        r = Ray([xf[p], yf[p], 0], [0, 0, 1])
        lens1.propagate_ray(r)
        lens2.propagate_ray(r)
        screen.propagate_ray(r)  
        x = r.p()[0]
        y = r.p()[1]
        dist_sq = x * x + y * y
        dlist.append(dist_sq)
    R = np.sqrt(sum(dlist) / len(dlist))
    return R 

#%%
# In the case of fixed output plane at f2 = 118.5 mm, plane first orientation
size2 = scipy.optimize.minimize(RMS_optimize, 
                                x0 = guess, bounds=[(0,1),(0,1)]) 



#%%
# BEWARE OF RUNTIME!

"""
Print:
1. For Plane at f1 = 121.7mm 

RMS_min = 0.003 mm at curvature combination of (0.018, 0.002) mm 

2. For Plane at f1 = 118.5mm 

RMS_min = 0.003 mm at curvature combination of (0.018, 1.600e-03) mm 
"""

P = []
index = []
c2 = []
opt1 = np.arange(0, 0.02, step = 2e-4) #100 points
opt2 = np.arange(0, 0.02, step = 2e-4)
for p in opt1:
    Q = []
    for q in opt2:
        Q.append(RMS_iteration(p, q, f1))
    P.append(min(Q))
    index.append(Q.index(min(Q)))
    c2.append(opt2[Q.index(min(Q))])
#%%
RMS_min1 = min(P)
c1_min1 = opt2[P.index(min(P))]
c2_min1 = c2[P.index(min(P))]

print('1. For Plane at f1 = 121.7mm \n')
print('RMS_min = %.3f mm at curvature combination of (%.3f, %.3f) mm \n'
      %(RMS_min1, c1_min1, c2_min1))

#%%
P = []
index = []
c2 = []
opt1 = np.arange(0, 0.02, step = 2e-4) #100 points
opt2 = np.arange(0, 0.02, step = 2e-4)
for p in opt1:
    Q = []
    for q in opt2:
        Q.append(RMS_iteration(p, q, f1))
    P.append(min(Q))
    index.append(Q.index(min(Q)))
    c2.append(opt2[Q.index(min(Q))])

RMS_min2 = min(P)
c1_min2 = opt2[P.index(min(P))]
c2_min2 = c2[P.index(min(P))]

print('2. For Plane at f1 = 118.5mm \n')
print('RMS_min = %.3f mm at curvature combination of (%.3f, %.3e) mm \n'
      %(RMS_min2, c1_min2, c2_min2))
    

    

    
    
    
    