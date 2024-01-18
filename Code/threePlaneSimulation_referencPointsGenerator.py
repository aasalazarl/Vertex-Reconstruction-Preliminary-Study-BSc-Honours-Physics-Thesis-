# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 16:46:34 2020

@author: alejandrosalazar
"""

# Import libraries.
import numpy as np
import math as math
#from mpl_toolkits.mplot3d import Axes3D

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Functions.

def average(list):
    return sum(list)/len(list)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Set up traces.

#--------------------------------------------------------------------------
# Dictionaries and variables.
vertices  = {}         # To store real vertex.
vertexReconstruction = 0   # To store vertex reconstruction.

planesIndices = ['x', 'y', 'z']
planeA = {'x': [], 'y': [], 'z': []}
planeB = {'x': [], 'y': [], 'z': []}
planeC = {'x': [], 'y': [], 'z': []}

planeA_recon = {'x': [], 'y': [], 'z': []}
planeB_recon = {'x': [], 'y': [], 'z': []}
planeC_recon = {'x': [], 'y': [], 'z': []}

gridsIndices = ['x', 'y']
gridA = {'x': [], 'y': []}
gridB = {'x': [], 'y': []}
gridC = {'x': [], 'y': []}


#--------------------------------------------------------------------------
# Real vertices coordinates.
# Use polar coordinates. Place traces plane in the y = 1.50 meters plane. 
# All units in meters.
x0 = list(np.linspace(1.50 - 0.02, 1.50 + 0.02, 1))
y0 = list(np.linspace(1.50, 1.50, len(x0)))
z0 = list(np.linspace(-2.50, -2.50, len(x0)))

vertices = {'x0': x0, 'y0': y0, 'z0': z0}

#--------------------------------------------------------------------------
# Intersection finder.

l_lower = -2; l_upper = 0

class intersectionFinder:
    def __init__(self, i, l_lower, l_upper, plane, grid, plane_recon,\
                 zPlane, z1, z2):

        for j in range(0, len(plane)):       

            plane_recon_element = str(planesIndices[j])        
                
            # Find intersections with xy-grid and points at the center 
            # of the squares for vertex reconstruction. 
            if z1 < zPlane + 0.000001 or z2 < zPlane + 0.000001:
                
                plane[str(planesIndices[j])].append(\
                      trace1[str(trace1Indices[j])][i])
                plane[str(planesIndices[j])].append(\
                      trace2[str(trace2Indices[j])][i])
            
                if j < len(planesIndices) - 1 and\
                (z1 < zPlane + 0.0001 or z2 < zPlane + 0.0001):    
                
                    for l in range(l_lower, l_upper):
                    
                        dummy = plane[str(planesIndices[j])][l]
                        grid_element = str(gridsIndices[j])

                        if  dummy > 0 and dummy < 0.08: 
                            segment = 1
                            grid[grid_element].append(segment)
                            plane_recon[plane_recon_element].append(\
                                       average([0, 0.08]))
                    
                        elif dummy >= 2.92 and dummy < 3:           
                            segment = 73
                            grid[grid_element].append(segment)
                            plane_recon[plane_recon_element].append(\
                                       average([2.92, 3]))
                    
                        elif dummy <= 0 or dummy >= 3:
                            grid[grid_element].append('No detection')
                            plane_recon[plane_recon_element].append(\
                                       'No detection')
                                
                        else:
                            segment = math.floor(((dummy - 0.08) * 100)\
                                                 /4 + 2)
                            grid[grid_element].append(segment)
                            plane_recon[plane_recon_element].append(\
                                       average([((segment - 2) * 4 + 8)\
                                                / 100, ((segment - 1) *\
                                                        4 + 8) / 100]))

                elif j == len(planesIndices) - 1:
                    plane_recon[plane_recon_element].append(z1)
                    plane_recon[plane_recon_element].append(z2)


#--------------------------------------------------------------------------
for m in range(0, len(x0)):
    
    l_lower += 2; l_upper += 2
    
    # Traces.
    r = np.linspace(0, 9, 1000000)
    
    # Trace 1.
    angle1 = 160 * (np.pi/180)      # Angle from z-axis.
    trace1Indices = ['x1', 'y1', 'z1']
    x1 = r * np.sin(angle1) + x0[m]
    y1 = np.linspace(1.50, 1.50, 1000000)
    z1 = r * np.cos(angle1) + z0[m]
    trace1 = {'x1': x1, 'y1': y1, 'z1': z1}

    # Trace 2.
    angle2 = 200 * (np.pi/180)
    trace2Indices = ['x2', 'y2', 'z2']
    x2 = r * np.sin(angle2) + x0[m]
    y2 = np.linspace(1.50, 1.50, 1000000)
    z2 = r * np.cos(angle2) + z0[m]
    trace2 = {'x2': x2, 'y2': y2, 'z2': z2}

    #----------------------------------------------------------------------
    # Find intersections of traces with planes.
    # Planes; A: zA = -5.00 m, B: zB = -13/2 m, C: zC = -8 m. 
    zA = -5.00; zB = -5.30; zC = -5.60 

    for i in range(0, len(r)):
        #------------------------------------------------------------------
        # Intersection with plane A.
        if abs(z1[i] + 5.00) < 0.00001 or abs(z2[i] + 5.00) < 0.00001: 
            intersectionFinder(i, l_lower, l_upper, planeA, gridA,\
                               planeA_recon, zA, z1[i], z2[i])
                
        #------------------------------------------------------------------
        # Intersection with plane B.
        if abs(z1[i] + 5.30) < 0.00001 or abs(z2[i] + 5.30) < 0.00001:
            intersectionFinder(i, l_lower, l_upper, planeB, gridB,\
                               planeB_recon, zB, z1[i], z2[i])
    
        #------------------------------------------------------------------    
        # Intersection with plane C.        
        if abs(z2[i] + 5.60) < 0.00001 or abs(z2[i] + 5.60) < 0.00001:
            intersectionFinder(i, l_lower, l_upper, planeC, gridC,\
                               planeC_recon, zC, z1[i], z2[i])

''' How to read the results:
    Example: i-coordinate: [a, b, c, d] (the size of i-coordinate will
    always be even).
    [a, b, ...] corresponds to decay # 1; a and b correspond to traces 
    1 and 2, respectively.
    [..., c, d] corresponds to decay # 2; c and d correspond to traces 
    1 and 2, respectively.'''
    
print('planeA:', planeA)
print('planeB:', planeB)
print('planeC:', planeC)
