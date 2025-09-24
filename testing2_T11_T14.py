#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 17:54:04 2021

@author: yisinuo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as plp
from matplotlib.patches import ConnectionPatch
from scipy.optimize import fmin_tnc 

from raytracer import *

#%%
# ------------------------------------------------------------------------
# Task 12
""" 
Thoughts on creating the bundle: 
Equal density means equal spacing between the points, so geometrically
the cross section of xy plane of the ray bundle would have its points forming
layers of circles at equal spacing to each other, and on each circle at layer 
n (at origin n = 0) there would be 6n points forming a regular polygon.

"""
opp1 = OutputPlane(200)
sp_sf = SphericalRefraction(100, 1, 1.5, 1/0.03, 0.03)        
opp = OutputPlane(250)
r0 = Ray([0,0,0],[0,0,1])

# Test Trace 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting the ray bundles along zy plane ~~~
fig, ax = plt.subplots(figsize = (6, 4), dpi = 100)
plt.plot(*opp1.draw_plane(), color = 'grey', label = 'Output Plane')
ax.add_patch(plp.Arc(*SphericalRefraction.draw_lens(sp_sf), 
                     color = 'turquoise', linewidth=2.5, label = 'Lens', 
                     alpha=0.4))
plt.xlabel('z / mm', fontsize = 12)
plt.ylabel('y / mm', fontsize = 12)
plt.title('Singlet f=100 mm Large Bundle Separation', 
          fontsize = 12, pad = 12)

# Creating bundles with diameter of 5mm 
# blue from the origin and red from (0, -7, 0)
beam = bundle(5, 5)
rb = Ray([0, 0, 0], [0, 0, 1])
rr = Ray([0, -7, 0], [0, 0.07, 1])

# Plot bundle using the bundle.plot() method added in class bundle()
beam.plot(rb, sp_sf, opp1, 'blue', 0.3)
beam.plot(rr, sp_sf, opp1, 'red', 0.3)
plt.ylim(-10, 6)
plt.legend(fontsize = 8, framealpha = 0.5)
plt.grid()
plt.show()

#%%
# Focus Calculation with graphical method -------------------------

# Method 1. defined focus() function in raytracer.py
# Which implement the graphical method with only one paraxial ray

# Method 2. use loop and mean 
beam = bundle(5, 5)
blue = beam.vertices(rb, sp_sf, opp1)
red = beam.vertices(rr, sp_sf, opp1)


fb = []
for i in range(0, len(blue)):
    p1_b = blue[i][-2] # y 
    p2_b = blue[i][-1]
    k1 = (p2_b - p1_b) 
    if k1[1] != 0.:
        a1 = -p1_b[1] / k1[1] # Factor by which ray travelled
        f1 = p1_b[-1] + a1 * k1[-1]
        fb.append(f1)
    else: 
        pass

focus_b = np.mean(fb)
std_b = np.std(fb)

# Comparison
print('---- \n') 

print('Focus calculation with graphical method \n') 

print('The focal point of Blue Large Bundle is at z = %.3f ± %.3f mm '
      'in the Figure: '
      'Singlet f=100 mm Large Bundle Separation\n'
      'In comparison, the expected value is 200 mm and the close paraxial ray\n'
      '(focus() function) gives a of %.3f mm.\n'
      %(focus_b, std_b, focus(sp_sf)))        
        
print('---- \n') 
#%%
# ------------------------------------------------------------------------
# Task 13
# Plotting the spot diagram
plt.figure(figsize = (7,7))
plt.suptitle('Corresponding Spot Diagram of Large Bundle Separation',
             fontsize = 16)
plt.subplot(221, aspect = 'equal')
plt.scatter(*beam.xy(rb), color = 'blue', s = 4)
plt.title('\n\nz=0', fontsize = 14, pad = 10)

plt.subplot(222, aspect = 'equal')
plt.scatter(*beam.xy(rr), color = 'red', s = 4)
plt.title('z=0', fontsize = 14, pad = 10)

# At the paraxial focus
# Blue bundle
xf1 = beam.vertices(rb, sp_sf, opp1)[:, -1, 0]
yf1 = beam.vertices(rb, sp_sf, opp1)[:, -1, 1]
plt.subplot(223, aspect = 'equal')
plt.scatter(xf1, yf1, color = 'blue', s = 4)
plt.title('focus', fontsize = 14, pad = 10)

# Red bundle
xf2 = beam.vertices(rr, sp_sf, opp1)[:, -1, 0]
yf2 = beam.vertices(rr, sp_sf, opp1)[:, -1, 1]
plt.subplot(224)
plt.scatter(xf2, yf2, color = 'red', s = 4)
plt.title('focus', fontsize = 14, pad = 10)
plt.tight_layout()

# RMS spot radius calculation~~~~~~
R1 = []
R2 = []
for i in np.arange(0, len(xf1)):
    R1.append(xf1[i] * xf1[i] + yf1[i] * yf1[i])
    R2.append(xf2[i] * xf2[i] + yf2[i] * yf2[i])

# RMS deviation from optical axis at the focus    
RMS1 = np.sqrt(sum(R1)/len(xf1))
RMS2 = np.sqrt(sum(R2)/len(xf2)) 

# Checkingadius of the geometrical focus
r_focus1 = (abs(max(yf1)-RMS1)+abs(min(yf1)-RMS1))/2
r_focus2 = (abs(max(yf2)-RMS2)+abs(min(yf2)-RMS2))/2


print('---- \n') 
print('Task 13: RMS Spot Radius \n')  
print('At Large Foci Separation')  
print('The RMS deviation from optical axis '
      'at the paraxial focus for Blue Bundle is %.3f mm, and the radius of \n'
      'the geometrical focus is %.3f mm.'
      %(RMS1, r_focus1))    
print('The RMS deviation from optical axis'
      ' at the paraxial focus for Red Bundle is %.3f mm, and the radius of \n'
      'the geometrical focus is %.3f mm.'
      %(RMS2, r_focus2))
print('The calculation of the radius of the geometrical focus agrees with \n'
      'the graph. \n') 

#%%
# ------------------------------------------------------------------------
# Task 14
print('---- \n') 
print('Task 14: Diffraction Limit\n') 

f = 100 # Focus length of 100 mm
D = 5 # Diameter of the bundle 5 mm

# For Bundle 1 of blue light with ~450nm
lambda_b = 450e-6 # Wavelength in mm
diffraction_limitb = lambda_b * f / D
print('The critical diffraction limit of blue ray bundle is %.3f mm, \n'
      'so the radius of the focus is %.0f times smaller than \n'
      'the diffraction limit.\n'
      %(diffraction_limitb, diffraction_limitb / r_focus1))

# For Bundle 2 of red light with ~700nm
lambda_r = 700e-6 # Wavelength in mm
diffraction_limitr = lambda_r * f / D
print('The critical diffraction limit of red ray bundle is %.3f mm, \n'
      'so the radius of the focus is around the same as \n'
      'the diffraction limit (ratio %.2f).\n'
      %(diffraction_limitr, diffraction_limitr / r_focus2))

print('Whether the two foci would be resoluted (distinguished) therefore \n'
      'depends on the diffraction limit of the higher wavelength (red) light\n'
      'which is %.3f mm.\n'%(diffraction_limitr))

# Separation of the two focus
d = abs(RMS1-RMS2)
print('The two rays would be resoluted at the paraxial focus given their \n'
      'separation d is %.2f mm and it is %s that d > diffraction limit.\n'
      %(d,d > diffraction_limitr))
print('Individually -----'
      'The blue ray bundle would not form distinguishable image as the RMS\n'
      'spot size is smaller than the critical diffraction limit.')

print('---- \n')

    


fig = plt.figure(figsize=(7, 7),dpi = 100)
plt.suptitle('Singlet f=100 mm \n'
          'Foci Separation at Diffraction Limit 0.014mm', fontsize = 15)

beam = bundle(5, 5)
rb = Ray([0, 0, 0], [0, 0, 1])
rr = Ray([0, -0.00239, 0], [0, 0.000239, 1]) # Reset the rays
oppb = OutputPlane(focus_b)
oppr = OutputPlane(focus_b) # Set up foci planes

# First zoom-in plot
start = fig.add_subplot(3,2,1) 
beam.plot(rb, sp_sf, oppb, 'blue', 0.6)
beam.plot(rr, sp_sf, oppr, 'red', 0.6)   
start.set_xlim(0, 25)
start.set_ylim(-0.005, 0.005)
plt.grid()

# Second zoom-in plot
foc = fig.add_subplot(3,2,2) 
beam.plot(rb, sp_sf, oppb, 'blue', 0.6)
beam.plot(rr, sp_sf, oppr, 'red', 0.3)
foc.plot(*oppb.draw_plane(), color = 'grey', label = 'Output Plane', 
           alpha = 1, linewidth = 2)
foc.set_ylim(-0.01, 0.026)
foc.set_xlim(199.5, 200.05)
plt.grid()

# The main bundle plot
main = fig.add_subplot(3,2,(3,6))
beam.plot(rb, sp_sf, oppb, 'blue', 0.6)
beam.plot(rr, sp_sf, oppr, 'red', 0.3)
main.plot(*oppb.draw_plane(), color = 'grey', label = 'Output Planes')
main.add_patch(plp.Arc(*SphericalRefraction.draw_lens(sp_sf), 
                     color = 'turquoise', linewidth=2.5, label = 'Lens', 
                     alpha=0.8, zorder = 2))

# Adding connection paths
con1 = ConnectionPatch(xyA=(0, -0.005), coordsA=start.transData, 
                       xyB=(0, 0), coordsB=main.transData, color = 'black')
con2 = ConnectionPatch(xyA=(25, -0.005), coordsA=start.transData, 
                       xyB=(25, 0), coordsB=main.transData, color = 'black')
con3 = ConnectionPatch(xyA=(199.5, -0.01), coordsA=foc.transData, 
                       xyB=(180, 0), coordsB=main.transData, color = 'black')
con4 = ConnectionPatch(xyA=(200.05, -0.01), coordsA=foc.transData, 
                       xyB=(205, 0), coordsB=main.transData, color = 'black')
fig.add_artist(con1)
fig.add_artist(con2)
fig.add_artist(con3)
fig.add_artist(con4)

# Making the graph look nice 
main.add_patch(plp.Rectangle((0, -0.3), 25, 0.6, fill = False, alpha=2, 
                           edgecolor = 'black', linewidth=1.5, zorder=10))
main.add_patch(plp.Rectangle((180, -0.5), 25, 1, fill = False, alpha=2, 
                           edgecolor = 'black', linewidth=1.5, zorder=10))
main.set_ylim(-4, 4)
main.set_xlabel('z / mm', fontsize = 12)
main.set_ylabel('y / mm', fontsize = 12)
plt.grid()
plt.legend(fontsize = 8, framealpha = 0.5)