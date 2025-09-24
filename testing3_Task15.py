#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 07:01:58 2021

@author: yisinuo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as plp
from matplotlib.patches import ConnectionPatch
from scipy.optimize import fmin_tnc 

from raytracer import *



#------------------------------------------------------------------------
# Task 15
# Separation 5mm, radius of conven lens 50mm means aperture radius ar:
ar = np.sqrt(50 * 50 - 45 * 45)

# 1. Plane first
convex_1 = ReverseSphericalRefraction(25, 1.5168, 1, ar, 0.02)
plane_1 = ReverseSphericalRefraction(20, 1, 1.5168, ar, 0)
f1 = focus(plane_1, convex_1) 
screen_1 = OutputPlane(f1) # Set output plane at focus

plt.figure(figsize = (7, 6), dpi = 100)
ax1 = plt.subplot(211)
plt.plot(*screen_1.draw_plane(), color = 'grey', linewidth = 3, 
         label = "Output Plane at Focus z=121.7mm")
plt.plot(*SphericalRefraction.draw_lens(plane_1), color = 'navy', 
         linewidth = 2.5, alpha = 0.5)
ax1.add_patch(plp.Arc(*SphericalRefraction.draw_lens(convex_1, concave = True), 
                     color = 'navy', linewidth = 2.5, 
                     label = 'Plano-Convex Lens Orientation 1', alpha = 0.5))
plt.ylim(-7.5, 7.5)
plt.xlim(10, 125)
#plt.xlim(118, 121)
#plt.ylim(-0.1, 0.1)

# Plotting the ray bundles along zy plane ~~~
plt.ylabel('y / mm', fontsize = 12)
plt.title('Bundle of Diameter 10mm through a Plano-Convex Lens\n'
          'L1 = Plane, L2 = Convex Lens', 
          fontsize = 14, pad = 10)

# Creating bundles with diameter of 10mm 
beam_10 = bundle(6, 10)
r_10 = Ray([0, 0, 0], [0, 0, 1])

# Plot bundle using the bundle.plot() method added in class bundle()
beam_10.plot(r_10, plane_1, screen_1, 'orange', 0.2, convex_1)
plt.legend(fontsize = 8, framealpha = 0.5)


# 2. Convex first
ax2 = plt.subplot(212, sharey = ax1, sharex = ax1)
convex = SphericalRefraction(20, 1, 1.5168, ar, 0.02)
plane = SphericalRefraction(25, 1.5168, 1, ar, 0)
f2 = focus(convex, plane) 
screen = OutputPlane(f2) 

plt.plot(*screen.draw_plane(), color = 'grey', linewidth = 3, 
         label = "Output Plane at Focus z=118.5mm")
plt.plot(*SphericalRefraction.draw_lens(plane), color = 'navy', 
         linewidth = 2.5, alpha = 0.5)
ax2.add_patch(plp.Arc(*SphericalRefraction.draw_lens(convex), 
                     color = 'navy', linewidth = 2.5, 
                     label = 'Plano-Convex Lens Orientation 2', alpha = 0.5))

# Plotting the ray bundles along zy plane ~~~
plt.xlabel('z / mm', fontsize = 12)
plt.ylabel('y / mm', fontsize = 12)
plt.title('L1 = Convex Lens, L2 = Plane', fontsize = 14, pad = 10)

# Plot bundle using the bundle.plot() method added in class bundle()
beam_10.plot(r_10, convex, screen, 'orange', 0.2, plane)
plt.legend(fontsize = 8, framealpha = 0.5)
plt.tight_layout()
plt.show()

#%%
plt.figure(figsize = (7, 6), dpi = 100)
ax1 = plt.subplot(211)
plt.plot(*screen_1.draw_plane(), color = 'grey', linewidth = 3, 
         label = "Output Plane at Focus z=121.7mm")
plt.plot(*SphericalRefraction.draw_lens(plane_1), color = 'turquoise', 
         linewidth = 2.5, alpha = 0.5)
ax1.add_patch(plp.Arc(*SphericalRefraction.draw_lens(convex_1, concave = True), 
                     color = 'turquoise', linewidth = 2.5, 
                     label = 'Plano-Convex Lens Orientation 1', alpha = 0.5))
plt.ylim(-0.0001, 0.0001)
plt.xlim(121.6, 122)
#plt.xlim(118, 121)
#plt.ylim(-0.1, 0.1)

# Plotting the ray bundles along zy plane ~~~
plt.ylabel('y / mm', fontsize = 12)
plt.title('Bundle of Diameter 0.2mm through a Plano-Convex Lens\n'
          'L1 = Plane, L2 = Convex Lens', 
          fontsize = 14, pad = 10)

# Creating bundles with diameter of 10mm 
beam_2 = bundle(6, 0.2)
r_10 = Ray([0, 0, 0], [0, 0, 1])

# Plot bundle using the bundle.plot() method added in class bundle()
beam_2.plot(r_10, plane_1, screen_1, 'orange', 0.2, convex_1)
plt.legend(fontsize = 8, framealpha = 0.5)


# 2. Convex first
ax2 = plt.subplot(212, sharey = ax1)
convex = SphericalRefraction(20, 1, 1.5168, ar, 0.02)
plane = SphericalRefraction(25, 1.5168, 1, ar, 0)
screen = OutputPlane(focus(convex, plane)) # Focus value calculated below

plt.plot(*screen.draw_plane(), color = 'grey', linewidth = 3, 
         label = "Output Plane at Focus z=118.5mm")
plt.plot(*SphericalRefraction.draw_lens(plane), color = 'turquoise', 
         linewidth = 2.5, alpha = 0.5)
ax2.add_patch(plp.Arc(*SphericalRefraction.draw_lens(convex), 
                     color = 'turquoise', linewidth = 2.5, 
                     label = 'Plano-Convex Lens Orientation 2', alpha = 0.5))
plt.xlim(118.4, 119)

# Plotting the ray bundles along zy plane ~~~
plt.xlabel('z / mm', fontsize = 12)
plt.ylabel('y / mm', fontsize = 12)
plt.title('L1 = Convex Lens, L2 = Plane', fontsize = 14, pad = 10)

# Plot bundle using the bundle.plot() method added in class bundle()
beam_2.plot(r_10, convex, screen, 'orange', 0.2, plane)
plt.legend(fontsize = 8, framealpha = 0.5)
plt.tight_layout()
plt.show()


#%%
# Resetting all objects
convex_1 = ReverseSphericalRefraction(25, 1.5168, 1, ar, 0.02)
plane_1 = SphericalRefraction(20, 1, 1.5168, ar, 0)
f1 = focus(plane_1, convex_1) 
screen_1 = OutputPlane(f1) # Set output plane at focus
convex = SphericalRefraction(20, 1, 1.5168, ar, 0.02)
plane = SphericalRefraction(25, 1.5168, 1, ar, 0)
f2 = focus(convex, plane) 
screen = OutputPlane(f2) 
beam_0 = bundle(6, 0.2)
beam_10 = bundle(6, 10)
r_10 = Ray([0, 0, 0], [0, 0, 1])
r0 = Ray([0, 0, 0], [0, 0, 1])

# Spot Diagram 2 ------------------------------------------------------------
plt.figure(figsize = (6, 7))
plt.suptitle('Corresponding Spot Diagram of Different Bundle Diameter D \n'
             'at Different Orientation of the Plano-Convex Lens',
             fontsize = 15)

# At the paraxial focus
# 1. Plane first
# D = 0.2mm
xf1 = beam_0.vertices(r_10, plane_1, screen_1, convex_1)[:, -1, 0] 
yf1 = beam_0.vertices(r_10, plane_1, screen_1, convex_1)[:, -1, 1]
plt.subplot(221, aspect = 'equal')
plt.scatter(xf1, yf1, c = '#ffcf33', s = 4)
plt.title('\n\n\n D=0.2mm ''\n O1=Plane First'
          '\nFocus at z=121.7mm', c = 'SaddleBrown', fontsize = 12, pad = 13)
plt.ylim(-50e-8, 50e-8)
plt.xlim(-50e-8, 50e-8)

# D = 10mm
xf2 = beam_10.vertices(r_10, plane_1, screen_1, convex_1)[:, -1, 0] 
yf2 = beam_10.vertices(r_10, plane_1, screen_1, convex_1)[:, -1, 1]
plt.subplot(222, aspect = 'equal')
plt.scatter(xf2, yf2, c = 'orange', s = 4)
plt.title('\n\n D=10mm ''\n O1=Plane First'
          '\nFocus at z=121.7mm', c = 'SaddleBrown', fontsize = 12, pad = 13)
plt.ylim(-0.02, 0.02)
plt.xlim(-0.02, 0.02)

# 2. Spherical lens First
# D = 0.2mm
xf3 = beam_0.vertices(r_10, convex, screen, plane)[:, -1, 0] 
yf3 = beam_0.vertices(r_10, convex, screen, plane)[:, -1, 1]
plt.subplot(223, aspect = 'equal')
plt.scatter(xf3, yf3, c = '#ffcf33', s = 4)
plt.title('\n O2=Convex First \n'
          'Focus at z=118.5mm ', c = 'SaddleBrown', fontsize = 12, pad = 13)
plt.ylim(-50e-8, 50e-8)
plt.xlim(-50e-8, 50e-8)

# D = 10mm
xf4 = beam_10.vertices(r_10, convex, screen, plane)[:, -1, 0] 
yf4 = beam_10.vertices(r_10, convex, screen, plane)[:, -1, 1]
plt.subplot(224, aspect = 'equal')
plt.scatter(xf4, yf4, color = 'orange', s = 4)
plt.title('\n O2=Convex First \n'
          'Focus at z=118.5mm ', c = 'SaddleBrown', fontsize = 12, pad = 13)
plt.ylim(-0.02, 0.02)
plt.xlim(-0.02, 0.02)
plt.tight_layout()


#%%
# RMS spot radius calculation -------------------------------------------
R1, R2, R3, R4 = [], [], [], []
for i in np.arange(0, len(xf1)):
    R1.append(xf1[i] * xf1[i] + yf1[i] * yf1[i])
    R2.append(xf2[i] * xf2[i] + yf2[i] * yf2[i])
    R3.append(xf3[i] * xf3[i] + yf3[i] * yf3[i])
    R4.append(xf3[i] * xf4[i] + yf4[i] * yf4[i])

# RMS deviation from optical axis at the focus    
RMS1 = np.sqrt(sum(R1)/len(xf1))
RMS2 = np.sqrt(sum(R2)/len(xf2)) 
RMS3 = np.sqrt(sum(R3)/len(xf3)) 
RMS4 = np.sqrt(sum(R4)/len(xf4)) 

#%%
f1_list = []
b1 = beam_0.vertices(r0, plane_1, screen_1, convex_1)    
for i in range(0, len(b1)):
    p1 = b1[i][-2] 
    p2 = b1[i][-1]
    k1 = (p2 - p1)                       
    if k1[1] != 0.:   
        a1 = -p1[1] / k1[1] # Factor by which ray travelled
        f = p1[-1] + a1 * k1[-1]
        f1_list.append(f)
    else: 
        pass
    
    
std1 = np.std(f1_list)
f1 = np.mean(f1_list)    
    
print('---- \n') 
print('Focus Calculation \n')
print('The focus of the plano-convex lens of curvature 0.02 mm at '
      'ORIENTATION 1 (PLANE FIRST) is %.1f mm \n'%(f1))

# 2. Spherical lens First
f2_list = []
b1 = beam_0.vertices(r0, convex, screen, plane)
for i in range(0, len(b1)):
    p1_1 = b1[i][-2] 
    p2_1 = b1[i][-1]
    k1 = (p2_1 - p1_1) 
    if k1[1] != 0.:
        a1 = -p1_1[1] / k1[1] # Factor by which ray travelled
        f = p1_1[-1] + a1 * k1[-1]
        f2_list.append(f)
    else: 
        pass
   
std2 = np.std(f2_list)
f2 = np.mean(f2_list)  

print('The focus of the plano-convex lens of curvature 0.02 mm at '
      'ORIENTATION 2 (CONVEX FIRST) is %.1f mm\n'%(f2))
print('---- \n') 

screen_1 = OutputPlane(f1)
screen = OutputPlane(f2)

#%%
# RMS spot radius calculation -------------------------------------------

def RMS_2lens(bundle, centralray = Ray([0, 0, 0], [0, 0, 1]), 
        lens1 = SphericalRefraction(1, 1, 1, 1000, 0), 
        screen = OutputPlane(100), lens2 = None):
    
    xf = bundle.vertices(centralray, lens1, screen, lens2)[:, -1, 0] 
    yf = bundle.vertices(centralray, lens1, screen, lens2)[:, -1, 1]
    RMS_list = []
    for i in np.arange(0, len(xf)):
        dist = xf[i] * xf[i] + yf[i] * yf[i]
        RMS_list.append(dist)
        
    return np.sqrt(sum(RMS_list)/(len(xf)-0))

# RMS of D = 10 and orientation O1
RMS_10_O1 = RMS_2lens(beam_10, r_10, plane_1, screen_1, convex_1)
# RMS of D = 10 and orientation O2
RMS_10_O2 = RMS_2lens(beam_10, r_10, convex, screen, plane)

# RMS of D = 0.2 and orientation O1
RMS_2_O1 = RMS_2lens(beam_0, r_10, plane_1, screen_1, convex_1)
# RMS of D = 10 and orientation O2
RMS_2_O2 = RMS_2lens(beam_0, r_10, convex, screen, plane)

# Diffraction Limit
D_15 = 10 # Diameter of the bundle 10 mm
lambda_15 = 588e-6 # Wavelength in mm
offset = (25 + 20) / 2 # z-offset of the lens

# 1. ORIENTATION 1 (PLANE FIRST)
f1_length = focus(plane_1, convex_1)  - offset
diffraction_limit1 = lambda_15 * f1_length / D_15

# 2. ORIENTATION 2 (CONVEX FIRST)
f2_length = focus(convex, plane) - offset
diffraction_limit2 = lambda_15 * f2_length / D_15


print('---- \n') 
print('RMS Spot Radius \n')  
print('1. ORIENTATION 1 (PLANE FIRST)------\n')  

print('RMS_2_O1 = %.3e mm '
      'at the paraxial focus for SMALL Bundle D = 0.2 mm.\n'
      %(RMS_2_O1,))   


print('RMS_10_O1 = %.4f mm '
      'at the paraxial focus for BIG Bundle D = 10 mm.\n'
      %(RMS_10_O1))

print('critical diffraction limit = %.4f mm for 588nm wavelength, \n'
      'so comparing to RMS of BIG Bundle %.3f mm it is smaller. \n'
      'The imaging would be very blurry for ORIENTATION 1.\n'
      %(diffraction_limit1, RMS_10_O1))

print('2. ORIENTATION 2 (CONVEX FIRST)------\n')  

print('RMS_2_O2 = %.3e mm '
      'at the paraxial focus for SMALL Bundle D = 0.2 mm.\n'
      %(RMS_2_O2,))  

print('RMS_10_O2 = %.4f mm \n'
      'at the paraxial focus for BIG Bundle D = 10 mm.\n'
      %(RMS_10_O2))


print('critical diffraction limit = %.4f mm for 588nm wavelength, \n'
      'so comparing to RMS of BIG Bundle %.3fmm it is around the same.\n'
      'The imaging would be less blurry for ORIENTATION 2.'
      %(diffraction_limit2, RMS_10_O2))

print('---- \n') 