#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri MAY 28 15:59:27 2021

@author: lin.yang
"""
import os
import numpy as np
import dump2energy as de
#
# k =list([0, 0.1, 0.5, 1, 2, 10])

fileBase="spparks"

base_name = ["000_a00", "0.2a00_1.2a00", "0.5a00_1.5a00", "0.7a00_1.7a00", "a00_2a00"]
# base_name1 = ["000_000", "0.5a00_0.5a00", "a00_a00", "1.5a00_1.5a00"]
# base_name = ["k0", "k0.1", "k0.5", "k1", "k2", "k10"]

path  = "/Users/lin/projects/SPPARKS-AGG/examples/Test_SimplifyIncE/2d_circle/"
file_name = "spparks"


for i in range(len(base_name)):


    ICfile ="circleIC_"+base_name[i]#, "circleIC_0.5a00_0.5a00", "circleIC_a00_a00", "circleIC_1.5a00_1.5a00"]

    figure_name1 = "c_"+base_name[i]+"_micro_2low_bp"
    figure_name2 = "c_"+base_name[i]+"_enrgy_2low_bp"

    os.system(f"mpirun -np 1 ~/projects/SPPARKS-AGG/src/spk_agg -var fileBase {fileBase} -var ICfile {ICfile} < agg_embedded.in")

    os.system(f'ffmpeg -framerate 30 -i Images/spparks.%04d.jpeg -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p video/{figure_name1}.mp4')

    figure_path="video/"+figure_name2
    timestep, energy_figure = de.dump2energy(path+file_name, 81)

    # plot_energy_figure(0, energy_figure)
    de.plot_energy_video(timestep, energy_figure, figure_path)

    os.system(f'ffmpeg -i video/{figure_name2}.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" video/{figure_name2}.mp4')
    os.system(f'rm video/{figure_name2}.gif')
