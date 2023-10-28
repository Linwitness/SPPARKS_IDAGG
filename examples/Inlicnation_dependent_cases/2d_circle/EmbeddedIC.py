#!/Users/lin.yang/miniconda3/bin/python python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  17 16:51:03 2020

@author: Lin.Yang
"""

import numpy as np
import math


savename = 'IC/circleIC_2a00_3a00.init' #The output filename
size = 512  #Number of sites in the x and y directions, so the size of the domain
img = np.zeros((size,size)) #Figure of all sites with GrainID

grains = 2 #Number of grains
Center = [size/2, size/2] #The center coordinate of the circle grain
Radius = 200 #The radius of the circle grain

#Randomly generate Euler Angles for each grain
EulerAngles = np.zeros((grains,3))
# for i in range(0,grains):
#    EulerAngles[i,0] = 2*math.pi*np.random.uniform(0,1)
#    EulerAngles[i,1] = 0.5*math.pi*np.random.uniform(0,1)
#    EulerAngles[i,2] = 2*math.pi*np.random.uniform(0,1)

#generate Euler Angles specifically for each grain
EulerAngles = np.zeros((grains,3))
EulerAngles[0,0] = 2*math.pi/2
EulerAngles[0,1] = 0
EulerAngles[0,2] = 0
EulerAngles[1,0] = 3*math.pi/2
EulerAngles[1,1] = 0
EulerAngles[1,2] = 0



#Paint each domain site according to distance from Center
for i in range(0,size):
    for j in range(0,size):
        GrainID = 0
        Dist = math.sqrt((Center[0]-j)**2+(Center[1]-i)**2)
        if Dist < Radius:
            GrainID = 2
        else:
            GrainID = 1
        img[j,i] = GrainID

#Write the information in the SPPARKS format and save the file
IC = [0]*(size*size+3)
IC[0] = '# This line is ignored\n'
IC[1] = 'Values\n'
IC[2] = '\n'
k=0
for i in range(0,size):
    for j in range(0,size):
        GrainID = int(img[j,i])
        IC[k+3] = str(k+1) + ' ' + str(int(GrainID)) + ' ' + str(EulerAngles[GrainID-1,0]) + ' ' + str(EulerAngles[GrainID-1,1]) + ' ' + str(EulerAngles[GrainID-1,2]) + '\n'
        k = k + 1

with open(savename, 'w') as file:
    file.writelines( IC )
