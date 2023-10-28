#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri MAY 28 15:59:27 2021

@author: lin.yang
"""
import os
import numpy as np
import dump2energy as de
import math

IWannaOnlyPlot = False

base_name = ["hex"]
tri_type = "ori_ave_aveE"

# inc_delta = [0.8265]
inc_delta = [0.6]
# inc_delta = [0.2]
# inc_delta = [0.0]

# J_scale = 1.844
# J_scale = 1.611
# J_scale = 1.203
J_scale = 1

reference_axis = [[1,0,0]]
# reference_axis = [[0.5*math.sqrt(3),0.5,0]]
# reference_axis = [[0.5*math.sqrt(2),0.5*math.sqrt(2),0]]
# reference_axis = [[0.5,0.5*math.sqrt(3),0]]
# reference_axis = [[0,1,0]]

# random_seed = [589943]
random_seed = [56689]
# random_seed = [2023]

kT = 0.66
core = 5

path  = os.getcwd() + "/"
if math.isnan(core):
    shell_file =  "agg_hex.in"
    core_tag = "_singleCore"
    IC_tag = ""
else:
    shell_file = "agg_hex_multiCore.in"
    core_tag = "_multiCore"+str(core)
    IC_tag = "_neighbor5"

path  = "/Users/lin.yang/projects/SPPARKS-AGG/examples/Test_SimplifyIncE/2d_hex/"

for i in range(len(base_name)):
    for k in range(len(inc_delta)):
        for m in range(len(reference_axis)):
            for n in range(len(random_seed)):
                # for o in range(len(triple_energy)):
                fileBase="h_"+tri_type+"_"+base_name[i]+core_tag+\
                    "_delta"+str(inc_delta[k])+"_m2_J"+str(J_scale)+\
                    "_refer_"+str(round(reference_axis[m][0],2))+"_"+\
                              str(round(reference_axis[m][1],2))+"_"+\
                              str(round(reference_axis[m][2],2))+\
                    "_seed"+str(random_seed[n])+"_kt066"

                ICfile =base_name[i]+"IC"+IC_tag

                figure_name1 = fileBase+"_micro"
                figure_name2 = fileBase+"_enrgy"

                if os.path.exists(path+fileBase+".dump") and IWannaOnlyPlot:
                    figure_path="video/"+figure_name2
                    timestep, energy_figure = de.dump2energy(path+fileBase, 81)
                    de.plot_energy_video(timestep, energy_figure, figure_path)
                else:
                    os.system(f"mpirun -np {1 if math.isnan(core) else core} ~/projects/SPPARKS-AGG/src/spk_agg \
                                -var fileBase {fileBase} \
                                -var ICfile {ICfile} \
                                -var inc_delta {inc_delta[k]} \
                                -var reference_axis0 {reference_axis[m][0]} \
                                -var reference_axis1 {reference_axis[m][1]} \
                                -var random_seed {random_seed[n]} \
                                -var kT {kT} \
                                -var J_scale {J_scale} \
                                    < {shell_file}")

                    os.system(f'ffmpeg -framerate 30 -i Images/{fileBase}.%04d.jpeg -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p video/{figure_name1}.mp4')

                    figure_path="video/"+figure_name2
                    timestep, energy_figure = de.dump2energy(path+fileBase, 81)

                    # plot_energy_figure(0, energy_figure)
                    de.plot_energy_video(timestep, energy_figure, figure_path)
