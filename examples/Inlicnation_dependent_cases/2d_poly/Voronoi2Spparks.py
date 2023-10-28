#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:10:03 2020

@author: fhilty
"""

import numpy as np
import math
from tqdm import tqdm

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False

savename = 'IC/VoronoiIC_4x3_elong_xy.init' #The output filename
size = [300, 400]  #Number of sites in the x and y directions, so the size of the domain
grains = 200  #Number of grains to create

#read grain centers from other file
# GCords = np.zeros((grains,2))
# img = np.zeros((size[0],size[1]))
# f = open("Case3.txt")
# line = f.readline()
# check = 0
# while line:
#     eachline = line.split()
#     if is_number(eachline[0]):
#         GCords[check,0],GCords[check,1] = float(eachline[0]), float(eachline[1])
#         check += 1
#     line = f.readline()
# f.close()
# print(GCords)
#Randomly generate grain centers
GCords = np.zeros((grains,2))
img = np.zeros((size[1],size[0]))
for i in range(0,grains):
    GCords[i,0],GCords[i,1]= np.random.randint(size[1]),np.random.randint(size[0])

#Paint each domain site according to which grain center is closest
for i in tqdm(range(0,size[0])):
    for j in range(0,size[1]):
        SiteID = 0
        if (GCords[0,0]-j) > size[1]/2:
            j_distance = size[1] - (GCords[0,0]-j)
        else:
            j_distance = (GCords[0,0]-j)
        if (GCords[0,1]-i) > size[0]/2:
            i_distance = size[0] - (GCords[0,1]-i)
        else:
            i_distance = (GCords[0,1]-i)
        MinDist = math.sqrt(j_distance**2+i_distance**2)
        # print(MinDist)
        for k in range(1,grains):
            if (GCords[k,0]-j) > size[1]/2:
                kj_distance = size[1] - (GCords[k,0]-j)
            else:
                kj_distance = (GCords[k,0]-j)
            if (GCords[k,1]-i) > size[0]/2:
                ki_distance = size[0] - (GCords[k,1]-i)
            else:
                ki_distance = (GCords[k,1]-i)
            dist = math.sqrt(kj_distance**2+ki_distance**2)
            if dist < MinDist:
                SiteID = k
                MinDist = dist
        img[j,i] = SiteID

#Randomly generate Euler Angles for each grain
EulerAngles = np.zeros((grains,3))
for i in range(0,grains):
   EulerAngles[i,0] = 0.5*math.pi*(np.random.rand()//0.5)
   EulerAngles[i,1] = 0
   EulerAngles[i,2] = 0

#Write the information in the SPPARKS format and save the file
IC = [0]*(size[0]*size[1]+3)
IC[0] = '# This line is ignored\n'
IC[1] = 'Values\n'
IC[2] = '\n'
k=0
for i in tqdm(range(0,size[0])):
    for j in range(0,size[1]):
        SiteID = int(img[j,i])
        IC[k+3] = str(k+1) + ' ' + str(int(SiteID+1)) + ' ' + str(EulerAngles[SiteID,0]) + ' ' + str(EulerAngles[SiteID,1]) + ' ' + str(EulerAngles[SiteID,2]) + '\n'
        k = k + 1
    # print(f"IC string updated: {i/size*100}%")

with open(savename, 'w') as file:
    file.writelines( IC )
