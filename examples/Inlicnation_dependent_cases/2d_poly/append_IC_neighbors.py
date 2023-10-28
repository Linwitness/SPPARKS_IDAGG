#!/Users/lin.yang/miniconda3/bin/python python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 18:51:03 2022

@author: Lin.Yang
"""

import numpy as np
from shutil import copyfile
import math

# The saved neighbor value is ID (1 ~ N) instead of i (0 ~ N-1))

# IC_filename = "IC/circleIC_000_000.init"
# IC_savename = "IC/circleIC_000_000_neighbors3.init"
IC_filename = "IC/VoronoiIC_20000.init"
IC_savename = "IC/VoronoiIC_20000_neighbor5.init"
size_x,size_y = 2400, 2400
interval = 5
nei_num = (2*interval+3)**2-1
img = np.zeros((size_y,size_x)) #Figure of all sites with GrainID

# copyfile(IC_filename, IC_savename)
with open(IC_filename, 'r') as f_read:
  IC_value = f_read.readlines()

for i in range(size_y): # y-axis
  for j in range(size_x): # x-axis
    img[i,j] = int(i*size_x + j)

IC_nei = []
IC_nei.append("# This line is ignored\n")
IC_nei.append("2 dimension\n")
IC_nei.append(f"{nei_num} max neighbors\n")
IC_nei.append(f"{size_x*size_y} sites\n")
IC_nei.append(f"0 {size_x} xlo xhi\n")
IC_nei.append(f"0 {size_y} ylo yhi\n")
IC_nei.append("0 1 zlo zhi\n")
IC_nei.append("\n")
IC_nei.append("Sites\n")
IC_nei.append("\n")
for i in range(size_y): # y-axis
  for j in range(size_x): # x-axis
    site = i*size_x + j
    IC_nei.append(f"{site + 1} {float(j)} {float(i)} 0.5\n")

IC_nei.append("\n")
IC_nei.append("Neighbors\n")
IC_nei.append("\n")

for i in range(size_y): # y-axis
  for j in range(size_x): # x-axis
    site = i*size_x + j

    # if nei_num == 8:
    #   im, ip = (i-1)%size_y, (i+1)%size_y
    #   jm, jp = (j-1)%size_x, (j+1)%size_x
    #   IC_nei.append(f"{site + 1} {int(img[im,jm]+1)} {int(img[im,j]+1)} {int(img[im,jp]+1)} "
    #                 f"{int(img[i,jm]+1)} {int(img[i,jp]+1)} "
    #                 f"{int(img[ip,jm]+1)} {int(img[ip,j]+1)} {int(img[ip,jp]+1)}\n")

    if nei_num > 0:
      tmp_nei = f"{site + 1}"
      for m in range(-(interval+1),interval+2):
        for n in range(-(interval+1),interval+2):
          if m==0 and n==0: continue
          tmp_i = (i+m)%size_y
          tmp_j = (j+n)%size_x
          tmp_nei += f" {int(img[tmp_i, tmp_j]+1)}"

      IC_nei.append(tmp_nei+"\n")

IC_nei.append("\n")

with open(IC_savename, 'w') as file:
    file.writelines( IC_nei )
    file.writelines( IC_value[1:] )
