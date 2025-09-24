#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 18:01:11 2021

@author: yisinuo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as plp
from matplotlib.patches import ConnectionPatch
from scipy.optimize import fmin_tnc 

def snell(incident_dir, normal, n1, n2):
    """
    Function implementing Snell's law to return the direction of the 
    refracted ray
    
    Parameters
    ----------
    incident_dir: direction of the incident ray 
    
    normal: surface normal of the lens surface
        
    n1: refractive index of the medium 1 (medium containing the incident ray)
        
    n2: refractive index of the medium 2
    
        
    Returns
    -------
    refracted_dir_norm 
        Direction of refracted ray
    
    """
    
    n_hat = -(normal / np.linalg.norm(normal))
    k = incident_dir
    ratio = n1 / n2
    dot = -np.dot(n_hat, k)
     
    # The function should return None in case of total internal reflection
    if np.sin(np.arccos(dot)) > (1/ratio):
        return None
    
    # Calculating the direction of the refracted ray using Snell's law in vector form
    else:
        refracted_dir = ratio * k + (ratio * dot - np.sqrt (1 - ratio * ratio * (1 - dot * dot))) * n_hat
        refracted_dir_norm = refracted_dir / np.linalg.norm(refracted_dir)
        return refracted_dir_norm
    
class Ray:
    
    """
    Class representing an optical ray.
    It is by convention that the optical axis is along the final argument 
    of the position array.
    
    Parameters
    ----------
    numpy.ndarray
       p: initial position of the ray
       k: initial direction of the ray 

    
    ----
    Methods defined here:
    p(): returns the current position of the ray
    
    k(): returns the current direction of the ray
    
    append(): append a new point and direction to the ray. p and k are 
    3-element arrays
    
    vertices(): returns all the points along the ray
    
    list(): returns all the points with their corresponding direction 
        in the ray
    
    draw_ray(): gives ray positions such that they can be put into 
            matplotlib.pyplot for plotting the ray along the optical axis
 
    """
    
    def __init__(self, p = [0.0,0.0,0.0], k = [0.0,0.0,0.0]):
        self._point = np.array(p, dtype = float)
        self._dir = np.array(k, dtype = float)
        self._plist = [[np.array(p), np.array(k)]]
        
    def p(self):
        """
        Method of returning most recent position of the ray

        Returns
        -------
        numpy.ndarray
            Most recent position of the ray.

        """
        return np.asarray(self._plist[-1][0], dtype = float)
    
    def k(self):
        """
        Method of returning most recent normalised direction of the ray
        
        Returns
        -------
        numpy.ndarray
            Current direction of the ray.

        """
        
        # For direction (0, 0, 0) (ray stops at the output plane)
        # or else would encounter ZeroDivisionError if try to normalise (0, 0, 0)
        if list(self._plist[-1][1]) == [0, 0, 0]:
            pass
        
        else:
            norm = np.linalg.norm(self._plist[-1][1]) #normalisation factor
            self._plist[-1][1] = self._plist[-1][1] / norm
        
        return np.asarray(self._plist[-1][1], dtype = float)
     
    def append(self, p, k):
        """
        Method of appending a new position and direction to the ray

        Parameters
        ----------
            p: 3D list of the position being added
            k: 3D list of the direction being added 

        Returns
        -------
        list
            A list of points and directions of the ray

        """
        self._plist.append([np.array(p, dtype = float), 
                            np.array(k, dtype = float)])
        return self._plist
    
    def vertices(self):
        pos = []
        for i in range(0,len(self._plist)):
            pos.append(self._plist[i][0])
        return np.array(pos)
    
    def list(self):
        """
        Method for checking the complete list of points and their 
        corresponding diretions. 
        
        Returns
        -------
        numpy.ndarray
            self._plist
        """
        return self._plist
    
    def draw_ray(self):
        
        """
        Method for returning the y, z positions of the propagated ray in two 
        separate list to enable ray-plotting using matplotlib.
        
        Raises
        ------
        Exception
            Prevents inputting to matplotlib.pyplot.plot without any 
            successful ray propagation.
        
        Returns
        -------
        numpy.ndarray
            self.vertices()[:,2] : z coordinates
            self.vertices()[:,1] : y coordinates
        
        """
        try:
            return self.vertices()[:,2], self.vertices()[:,1]
        except:
            raise Exception('ray propagation failed')
            
class OpticalElement:
    """
    Base class representing an optical element responsible for propagating 
    the ray through. 
    
    ----
    Methods defined here:
        propagate_ray(): Base method for ray propagation. 
            Each optical element would have a different implementation.
            
    """
    def propagate_ray(self, ray):
        """
        Parameters
        ----------
        ray : __main__.Ray
            Ray() class object being drawn.

        Raises
        ------
        NotImplementedError

        Returns
        -------
        None.

        """
        
        "propagate a ray through the optical element"
        raise NotImplementedError()
        
        
            
class SphericalRefraction(OpticalElement):
    """
    Class representing a spherical refracting surface.
    
    Parameters
    ----------
    Int or float
        z_intercept: the position of interception of the optical axis
        
        n1: refractive index of the medium 1 
            (medium containing the incident ray)
        
        n2: refractive index of the medium 2 (of the lens)
        
        aperture_radius: radius of aperture
        
        curvature: curvature of the lens
    
    ----
    Methods defined here:
        radius(): returns the radius of the spherical lens and organise the 
            special case of flat plane
        
        intercept(): returns the first valid intercept of the ray to the 
            spherical surface
    
        propagate_ray(): returns the new refracted dircection of the ray
        
        draw_ray(): gives ray positions after propagation such that they can 
            be put into matplotlib.pyplot for plotting the ray 
            along the optical axis
      
    
    """
    def __init__(self, z_intercept = 0, n1 = 1, n2 = 1, aperture_radius = 0, curvature = 0):
        self._z0 = np.array([0, 0, z_intercept], dtype = float)
        self._cur = curvature
        self._n1 = n1
        self._n2 = n2
        self._ar = aperture_radius
        
        # n1 and n2 need to be more than 1
        if n1 < 1 or n2 < 1:
            raise ValueError('The refractive index n1 or n2 '
                             'cannot be less than 1')
            
        # The aperture radius would be trivial if greater than the radius of curvature
        if self._cur != 0:
            if self._ar > 1/abs(self._cur):
                raise ValueError('The aperture radius is expected to be '
                                 'smaller than the radius of curvature')
    
    def radius(self):
        """
        Method of calculating the radius of the spherical lens for the 
        convienience of in-class methods.

        Returns
        -------
        float
            Radius of spherical lens

        """
        if self._cur == 0:
            return None
        else:
            return 1 / self._cur
                
    def intercept(self, ray):
        """.
        Method of calculating the first valid intercept of the ray to the lens
        surface

        Parameters
        ----------
        ray : __main__.Ray
            Ray() class object being drawn.

        Returns
        -------
        itc : np.ndarray
            Position of interception.

        """
        k = ray.k()
        # Intercept for a plane surface
        if self.radius() == None: 
            a = (self._z0[-1] - ray.p()[-1]) / k[-1] # The factor by which the ray propagates before interception
            itc = ray.p() + a * k
            
            return itc
        
        # Intercept for a curved surface
        else: 
            R = self.radius()
            centre = self._z0 + np.array([0, 0, R], dtype = float)
            r = ray.p() - centre
            root = (np.dot(r,ray.k()) * np.dot(r,ray.k())) - (np.dot(r,r) - R * R)
            
            # Return None for invalid intercept
            if root < 0:
                return None
            
            # Picking out the first intercept
            l1 = -np.dot(r, ray.k()) + np.sqrt(root)
            l2 = -np.dot(r, ray.k()) - np.sqrt(root)
            #l = min(abs(l1), abs(l2))
            l = min(abs(l1), abs(l2))
            itc = ray.p() + l * ray.k()

            return itc
        
        # Checking validity of the interception with aperture radius on y-axis
        if abs(itc[1]) > self._ar: 
            return None
           
    def draw_ray(self, ray):
        """
        Method for returning the y, z positions of the propagated ray in two 
        separate list to enable ray-plotting using matplotlib.

        Parameters
        ----------
        ray : __main__.Ray
            Ray() class object being drawn.

        Returns
        -------
        list
            ray_z : z-coordinates of the ray positions
            ray_y : y-coordinates of the ray positions
        """
        
        
        # A straight line in the case of no valid interception
        itc = self.intercept(ray)
        if itc is None:
            ray_z = [ray.p()[-1], ray.k()[-1] * 1000]
            ray_y = [ray.p()[1], ray.k()[1] * 1000]
        
        else:
            pos = ray.vertices()
            ray_z = []
            ray_y = []
            for i in range(0, len(pos)):
                 ray_z.append(pos[i][-1])
                 ray_y.append(pos[i][1])
        return ray_z, ray_y

        
        
    def propagate_ray(self, ray):
        
        """
        Method of propagating a ray through a lens, updating a new direction 
        of the ray
        
        Parameters
        ----------
        ray : __main__.Ray
            Ray() class object being propagated.
        
        Returns
        -------
        numpy.ndarray 
            New direction of the ray after refraction
        
        """
        
        k = ray.k()
        itc = self.intercept(ray)
        
        # When there is no valid intercept, the method returns None
        if itc is None:
            return None
        
        # The surface normal is special for a plane surface
        if self.radius() == None:
            normal = np.array([0., 0., 1.])
            
            
        else:
            centre = self._z0 + np.array([0, 0, self.radius()], dtype = float)
            normal = centre - itc

        
        newk = snell(k, normal, self._n1, self._n2)
        
        # When there is total internal refraction, the method returns None
        if newk is None:
            return None, normal
        
        else:
            ray.append(itc, newk)
            return newk
        
    def draw_lens(self, concave = False):
        """
        Method of converting the parameters of spherical lens of 
        SphericalRefraction into parameters used in matplotlib.patches.Arc so 
        as to enable easier plotting of the lens.
        
        Parameters
        ----------
        concave : bool, optional
            Set orientation of the spherical lens
        
        Returns
        -------
        float 
            Key parameters for matplotlib.patches.Arc
        
        """
        if self.radius() != None:
            R = self.radius() 
            # R would be the "width" and "height" parameters 
            # in matplotlib.patches.Arc
          
            if concave == False:
                x = self._z0[-1] + R
                y = 0
                xy = (x, y)
                angle = (np.arcsin(self._ar / R)) * 180 / np.pi # in deg
                theta1 = 180 - angle
                theta2 = theta1 + 2 * angle
                return xy, 2*R, 2*R, 0, theta1, theta2
            else:
                x = self._z0[-1] - R
                y = 0
                xy = (x, y)
                angle = (np.arcsin(self._ar / R)) * 180 / np.pi # in deg
                theta1 = 180 - angle
                theta2 = theta1 + 2 * angle
                return xy, 2*R, 2*R, 0, theta2, theta1
                
        else: # Plane lens
            x = np.empty(2)
            x.fill(self._z0[-1])
            y = [-self._ar, self._ar]
            return x, y 
            
    def focal(self):
        """
        Method for finding the position of focal point of one spherical lens

        Returns
        -------
        numpy.ndarray
            Position of focal point

        """
        f = 1 / ((self._n2 - 1)*(self._cur))
        return np.array([0., 0., self._z0[-1] + self.radius() + f])
    

class ReverseSphericalRefraction(SphericalRefraction):
    """
    Class representing a spherical refracting surface but in opposite 
    orientation to the normal Spherical Refraction class.
    
    Parameters
    ----------
    Int or float
        z_intercept: the position of interception of the optical axis
        
        n1: refractive index of the medium 1 
            (medium containing the incident ray)
        
        n2: refractive index of the medium 2 (of the lens)
        
        aperture_radius: radius of aperture
        
        curvature: curvature of the lens
    
    ----
    Methods defined here:
        intercept(): returns the second valid intercept of the ray to the 
        spherical surface
    
        propagate_ray(): returns the new refracted dircection of the ray
            
    
    """
                
    def intercept(self, ray):
        """
        Method of calculating the second valid intercept of the ray to the lens
        surface

        Parameters
        ----------
        ray : __main__.Ray
            Ray() class object being drawn.

        Returns
        -------
        itc : np.ndarray
            Position of interception.

        """
        k = ray.k()
        
        # Intercept for a plane surface
        if self.radius() == None: 
            dist = self._z0[-1] - ray.p()[-1]
            if dist < 0 and abs(self._z0[-1]) < abs(ray.p()[-1]): 
                return None #No interception if ray does not go into the plane
            else:
                # The factor by which the ray propagates before interception
                a = (self._z0[-1] - ray.p()[-1]) / k[-1] 
                itc = ray.p() + a * k
        
        # Intercept for a curved surface
        else: 
            R = self.radius()
            centre = self._z0 - np.array([0, 0, R], dtype = float)
            r = ray.p() - centre
            root = (np.dot(r,ray.k()) * np.dot(r,ray.k())) - (np.dot(r,r) 
                                                              - R * R)
            # Return None for invalid intercept
            if root < 0:
                return None

            l1 = -np.dot(r, ray.k()) + np.sqrt(root)
            l2 = -np.dot(r, ray.k()) - np.sqrt(root)
        
            # Picking out the second intercept
            l = max(l1, l2)
            itc = ray.p() + l * ray.k()
        return itc
        
        # Checking validity of the interception with aperture radius on y-axis
        if abs(itc[1]) > self._ar: 
            return None
        
        
    def propagate_ray(self, ray, reverse = None):
        
        """
        Method of propagating a ray through a lens, updating a new direction 
        of the ray
        
        Parameters
        ----------
        ray : __main__.Ray
            Ray() class object being propagated.
        
        Returns
        -------
        numpy.ndarray 
            New direction of the ray after refraction
        
        """
        
        k = ray.k()
        itc = self.intercept(ray)
        
        # When there is no valid intercept, the method returns None
        if itc is None:
            return None
        
        # The surface normal is special for a plane surface
        if self.radius() == None:
                normal = np.array([0., 0., 1.])

        else:
            centre = self._z0 - np.array([0, 0, self.radius()], dtype = float)
            normal = itc - centre
        newk = snell(k, normal, self._n1, self._n2)
        
        # When there is total internal refraction, the method returns None
        if newk is None:
            return None
        
        else:
            ray.append(itc, newk)
            return newk
        
    def draw_lens(self):
        """
        Method of converting the parameters of spherical lens of 
        ReverseSphericalRefraction into parameters used in
        matplotlib.patches.Arc so as to enable easier plotting of the lens.

        Returns
        -------
        float 
            Key parameters for matplotlib.patches.Arc
        
        """
        if self.radius() != None:
            R = self.radius() 
            # R would be the "width" and "height" parameters 
            # in matplotlib.patches.Arc
            x = self._z0[-1] - R
            y = 0
            xy = (x, y)
            angle = (np.arcsin(self._ar / R)) * 180 / np.pi # in deg
            theta1 = 180 - angle
            theta2 = theta1 + 2 * angle
            return xy, 2*R, 2*R, 0, theta2, theta1
                
        else: # Plane lens
            x = np.empty(2)
            x.fill(self._z0[-1])
            y = [-self._ar, self._ar]
            return x, y 
            
            
class OutputPlane(OpticalElement):
    """
        Class representing a output plane where the ray terminates.

        Parameters
        ----------
        Int or float
            z_intercept: the position of interception of the optical axis

         ----
        Methods defined here:
        intercept(): returns the first valid intercept of the ray to the 
            spherical surface
    
        propagate_ray(): appends to the ray the termination point and 
            direction

        """
    def __init__(self, z_intercept):
        self._z0 = np.array([0., 0., z_intercept], dtype = float)
    
    def intercept(self, ray):
        """
        Method of calculating the valid intercept of the ray to the plane

        Parameters
        ----------
        ray : __main__.Ray
            Ray() class object being drawn.

        Returns
        -------
        itc : np.ndarray
            Position of interception.

        """
        k = ray.k()
        dist = self._z0[-1] - ray.p()[-1]
        if dist < 0 and abs(self._z0[-1]) < abs(ray.p()[-1]): 
            return None # No interception if ray does not go into the plane
        else:
            a = dist / k[-1] 
            # The factor by which the ray propagates before interception
            itc = ray.p() + a * ray.k()
            return itc
    
    def propagate_ray(self, ray):
        """
        Method of terminating a ray at the plane. Append the terminating
        direction and the intercept point to the ray.
        
        Parameters
        ----------
        ray : __main__.Ray
            Ray() class object being propagated.
                
        """
    
        ray.append(self.intercept(ray), [0., 0., 0])
        
    def draw_plane(self):
        """
        Method of converting the parameters of output plane into parameters 
        used in matplotlib.pyplot.plot so as to enable easier plotting 
        of the lens.
        
        Parameters
        ----------
        concave : bool, optional
            Set orientation of the spherical lens
        
        Returns
        -------
        float 
            Key parameters for matplotlib.patches.Arc
        
        """
        x = np.empty(2)
        x.fill(self._z0[-1])
        y = [-1000, 1000]
        return x, y
    
    
    
class OptimizeLens(SphericalRefraction):
    """
    Class representing a spherical refracting surface used for lens 
    optimization only, for it does not raise ValueError when the aperture 
    radius is bigger than the radius of curvature.
    
    Parameters
    ----------
    Int or float
        z_intercept: the position of interception of the optical axis
        
        n1: refractive index of the medium 1 
            (medium containing the incident ray)
        
        n2: refractive index of the medium 2 (of the lens)
        
        aperture_radius: radius of aperture
        
        curvature: curvature of the lens
            
    
    """
    def __init__(self, z_intercept = 0, n1 = 1, n2 = 1, 
                 aperture_radius = 0, curvature = 0):
        self._z0 = np.array([0, 0, z_intercept], dtype = float)
        self._cur = curvature
        self._n1 = n1
        self._n2 = n2
        self._ar = aperture_radius
        
        # n1 and n2 need to be more than 1
        if n1 < 1 or n2 < 1:
            raise ValueError('The refractive index n1 or n2 cannot be less than 1')
            
class ReverseOptimizeLens(ReverseSphericalRefraction):
    """
    Class representing a spherical refracting surface of opposite orientation
    used for lens optimization only, for it does not raise ValueError when 
    the aperture radius is bigger than the radius of curvature.
    
    Parameters
    ----------
    Int or float
        z_intercept: the position of interception of the optical axis
        
        n1: refractive index of the medium 1 
            (medium containing the incident ray)
        
        n2: refractive index of the medium 2 (of the lens)
        
        aperture_radius: radius of aperture
        
        curvature: curvature of the lens
            
    
    """
    def __init__(self, z_intercept = 0, n1 = 1, n2 = 1, aperture_radius = 0, curvature = 0):
        self._z0 = np.array([0, 0, z_intercept], dtype = float)
        self._cur = curvature
        self._n1 = n1
        self._n2 = n2
        self._ar = aperture_radius
        
        # n1 and n2 need to be more than 1
        if n1 < 1 or n2 < 1:
            raise ValueError('The refractive index n1 or n2 cannot be less than 1')
            
            
class bundle():
    def __init__(self, layer_num, diameter):
        self._n = layer_num
        self._D = diameter
        
    def xy(self, centralray = Ray([0, 0, 0], [0, 0, 1])):
        """
        A function for generating the Cartesian coordinates in the xy plane of 
        each point in a uniform collimated beam.
    
        Parameters
        ----------
        centralray : __main__.Ray, optional
            Ray object that is at the centre of the bundle. 
            The default is Ray([0, 0, 0], [0, 0, 1]).
        
        Returns
        -------
        numpy.ndarray
            x: x-coordinates of all points in the beam
            y: y-coordinates of all points in the beam
    
        """
        r = self._D / (2 * self._n) # Spacing between the points
        xc = centralray.vertices()[0][0]
        yc = centralray.vertices()[0][1]
        x = [0.+xc]
        y = [0.+yc]
        
        for i in np.arange(0, self._n+1):
            l = i * r # Distance of the point from the central point
            N = 6 * i # Number of sides and points in the regular polygon
            for k in np.arange(1, N+1):
                theta = 2 * (k-1) * np.pi / N
                x.append(l * np.cos(theta)+xc)
                y.append(l * np.sin(theta)+yc)
                
        return np.array(x), np.array(y)
    
    def plot(self, centralray = Ray([0, 0, 0], [0, 0, 1]), 
                    lens1 = SphericalRefraction(1, 1, 1, 1000, 0), 
                    screen = OutputPlane(100), c = 'blue', l=1,
                    lens2 = None):

        """
        Method for plotting the bundle along the optical axis through a lens 
        and to a output plane.
        
        Parameters
        ----------
        centralray : __main__.Ray, optional
            Ray object that is at the centre of the bundle. 
            The default is Ray([0, 0, 0], [0, 0, 1]).
        
        lens1 : __main__.SphericalRefraction, optional
            The lens through which the bundle propagates.
            The default creates an SphericalRefraction object - plane lens at 
            z = 0 of height 1000 and refractive index 1.
            
        lens2 : __main__.SphericalRefraction, optional
            The second lens through which the bundle propagates.
            The default is None as the second lens is optional.
        
        screen :__main__.OutputPlane, optional
            The default is OutputPlane(100).
            
        c : str, optional
            Input the color of the ray. The default is 'blue'.
            
        l : int, optional
            Input the linewidth of the ray. The default is 1.
        
        """
        k = centralray.list()[0][1]
        #k = centralray.k()
        z = centralray.list()[0][0][-1]
        #z = k[-1]
        xy_list = self.xy(centralray)
        for i in np.arange(0, len(xy_list[0])):
            r = Ray([xy_list[0][i], xy_list[1][i], z], k)
            if lens2 == None: 
                #if reverse == None:
                    lens1.propagate_ray(r) 
                #else:
                    #lens1.propagate_ray(r, reverse = 1)
                #screen.propagate_ray(r) 
                #plt.plot(*r.draw_ray(), color = str(c), linewidth = l)
            else: 
                #if reverse == None:
                    lens1.propagate_ray(r) 
                    lens2.propagate_ray(r) 
                #else:
                    #lens1.propagate_ray(r, reverse = 1) 
                    #lens2.propagate_ray(r, reverse = 1)
            screen.propagate_ray(r)
            plt.plot(*r.draw_ray(), color = str(c), linewidth = l) 
                
    
    def vertices(self, centralray = Ray([0, 0, 0], [0, 0, 1]), 
                    lens1 = SphericalRefraction(1, 1, 1, 1000, 0), 
                    screen = OutputPlane(100), lens2 = None):
        """
        Method for retuning the positions of all the points in the bundle.
        
        Parameters
        ----------
        centralray : __main__.Ray, optional
            Ray object that is at the centre of the bundle. 
            The default is Ray([0, 0, 0], [0, 0, 1]).
        
        lens1 : __main__.SphericalRefraction, optional
            The lens through which the bundle propagates.
            The default creates an SphericalRefraction object - plane lens at 
            z = 0 of height 1000 and refractive index 1.
            
        lens2 : __main__.SphericalRefraction, optional
            The second lens through which the bundle propagates.
            The default is None as the second lens is optional.
        
        screen :__main__.OutputPlane, optional
            The default is OutputPlane(100).
        
        Returns
        -------
        numpy.ndarray
            Vertices of all the rays in the bundle.
        
        """
        k = centralray.list()[0][1]
        z = centralray.list()[0][0][-1]
        xy_list = self.xy(centralray)
        vertices = []
        
        for i in np.arange(0, len(xy_list[0])):
            r = Ray([xy_list[0][i], xy_list[1][i], z], k)
    
            if lens2 == None: 
                #if reverse == None:
                    
                lens1.propagate_ray(r)
                #else: 
                    #lens1.propagate_ray(r, reverse = 1)
                screen.propagate_ray(r)  
                vertices.append(r.vertices())
            else: 
                #if reverse == None:
                lens1.propagate_ray(r) 
                lens2.propagate_ray(r)
                #else: 
                    #lens1.propagate_ray(r, reverse = 1) 
                    #lens2.propagate_ray(r, reverse = 1)
                screen.propagate_ray(r) 
                vertices.append(r.vertices())
                
        return np.array(vertices)
    
        
def focus(lens1, lens2 = None):
    """
    

    Parameters
    ----------
    lens1 : __main__.SphericalRefraction
        Primary lens surface.
    lens2 : __main__.SphericalRefraction, optional
        Secondary lens surface. The default is None.

    Returns
    -------
    float
        f
        focus position along the optical axis.

    """
    parax = Ray([0, 0.01, 0],[0, 0, 1])
    lens1.propagate_ray(parax)
    if lens2 == None:
        pass
    else:
        lens2.propagate_ray(parax)
    p = parax.list()[-1][0]
    k = parax.list()[-1][-1]                      
    a = abs(p[1]) / abs(k[1]) # Factor by which ray travelled
    f = p[-1] + a * k[-1]
    return f

def RMS_2lens(bundle, centralray = Ray([0, 0, 0], [0, 0, 1]), 
        lens1 = SphericalRefraction(1, 1, 1, 1000, 0), 
        screen = OutputPlane(100), lens2 = None):
    """
    Function that returns the RMS spot size for any bundle of ray, 
    at any given output plane, through any primary or secondary lens.

    Parameters
    ----------
    bundle : __main__.bundle
        The bundle of ray which is being propagated

    centralray : __main__.Ray, optional
        The centre of the bundle.
        The default is Ray([0, 0, 0], [0, 0, 1]).
        
    lens1 : __main__.SphericalRefraction, optional
        Primary optical element.
        The default is SphericalRefraction(1, 1, 1, 1000, 0).
        
    screen : __main__.OutputPlane, optional
        The default is OutputPlane(100).
        
    lens2 : __main__.SphericalRefraction, optional
        Secondary optical element.
        The default is None.

    Returns
    -------
    float
        RMS spot size.

    """
    
    xf = bundle.vertices(centralray, lens1, screen, lens2)[:, -1, 0] 
    yf = bundle.vertices(centralray, lens1, screen, lens2)[:, -1, 1]
    RMS_list = []
    for i in np.arange(0, len(xf)):
        dist = xf[i] * xf[i] + yf[i] * yf[i]
        RMS_list.append(dist)
        
    return np.sqrt(sum(RMS_list)/(len(xf)-0))

def RMS_iteration(c1, c2, foc):
    """
    Function that returns the RMS spot size for specific ray bundle of 10mm
    and through specific lens of fixed refractive index and position. This is 
    to aid optimization of the special case in the script.

    Parameters
    ----------
    c1 : float or int.
        Curvature of the primary lens.
    c2 : float or int.
        Curvature of the secondary lens.
    foc : float or int.
        Fixed position of the screen.

    Returns
    -------
    R : float
        RMS spot size.

    """
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