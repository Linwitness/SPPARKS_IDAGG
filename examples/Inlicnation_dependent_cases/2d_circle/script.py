#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri MAY 28 15:59:27 2021

@author: lin.yang
"""
from cmath import nan
import os
from random import random
import numpy as np
import dump2energy as de
import math

IWannaOnlyPlot = False

base_name = ["000_000"]
# base_name = ["embeddedIC"]


# inc_delta = [0.8265]
inc_delta = [0.95]
# inc_delta = [0.20]
# inc_delta = [0.0]

# energy_scale = [1.825]
# energy_scale = [1.6]
# energy_scale = [1.203]
energy_scale = [1]

reference_axis = [[1,0,0]]
# reference_axis = [[0.5*math.sqrt(3),0.5,0]]
# reference_axis = [[0.5,0.5*math.sqrt(3),0]]
# reference_axis = [[0,1,0]]
# reference_axis = [[0.5*math.sqrt(2),0.5*math.sqrt(2),0]]

random_seed = [56689]
# random_seed = [30]
# random_seed = [827304]
# random_seed = [2022]
# random_seed = [7374]

core = 16

path  = os.getcwd() + "/"
if math.isnan(core):
    shell_file =  "agg_embedded_1Core.in"
    core_tag = "_singleCore"
    IC_tag = ""
else:
    shell_file = "agg_embedded_multiCore.in"
    core_tag = "_multiCore"+str(core)
    IC_tag = "_neighbors5"



for i in range(len(base_name)):
        for k in range(len(inc_delta)):
            for m in range(len(reference_axis)):
                for n in range(len(energy_scale)):
                    for o in range(len(random_seed)):
                        fileBase="cT_ori_aveE_"+base_name[i]+core_tag+"_kt066"+"_seed"+str(random_seed[o])+\
                                "_scale"+str(energy_scale[n])+"_delta"+str(inc_delta[k])+"_m6"+\
                                "_refer_"+str(round(reference_axis[m][0],2))+"_"+\
                                        str(round(reference_axis[m][1],2))+"_"+\
                                        str(round(reference_axis[m][2],2))

                        ICfile ="circleIC_"+base_name[i]+IC_tag

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
                                        -var energy_scale {energy_scale[n]} \
                                        -var random_seed {random_seed[o]} \
                                            < {shell_file}")

                            os.system(f'ffmpeg -framerate 30 -i Images/{fileBase}.%04d.jpeg \
                                        -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p \
                                        video/{figure_name1}.mp4')

                            figure_path="video/"+figure_name2
                            timestep, energy_figure = de.dump2energy(path+fileBase, 81)

                            de.plot_energy_video(timestep, energy_figure, figure_path,0.1)
                            np.save("results/"+fileBase+"_energy",energy_figure)

                            if not os.path.exists("results/"+fileBase+".npy"):
                                timestep, grain_structure_figure = de.dump2img(path+fileBase, 81)
                                np.save("results/"+fileBase,grain_structure_figure)
                            else:
                                grain_structure_figure = np.load(path + "results/"+fileBase+".npy")
                                timestep = 30 * np.array(range(len(grain_structure_figure)))
                            de.plot_structure_video(timestep, grain_structure_figure, "video/"+fileBase+"_structure")
