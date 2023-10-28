# -*- coding: utf-8 -*-


import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm

def plot_energy_figure(timestep, energy_figure, figure_path=None):

    imgs = []
    fig, ax = plt.subplots()

    cv0 = np.squeeze(energy_figure[timestep])
    cv0 = np.rot90(cv0,1)
    im = ax.imshow(cv0,vmin=0,vmax=6,cmap='Accent')
    cb = fig.colorbar(im)
    tx = ax.set_title(f'time step = {timestep}')
    if figure_path != None:
        plt.savefig('energy_{timestep}step')
    plt.show()

def plot_energy_video(timestep, energy_figure, figure_path):

    imgs = []
    fig, ax = plt.subplots()

    cv0 = np.squeeze(energy_figure[0])
    cv0 = np.rot90(cv0,1)
    im = ax.imshow(cv0,vmin=np.min(cv0),vmax=np.max(cv0),cmap='Accent')
    cb = fig.colorbar(im)
    tx = ax.set_title(f'time step = {timestep[0]}')
    # plt.show()

    def animate(i):
        arr=np.squeeze(energy_figure[i])
        arr=np.rot90(arr,1)
        im.set_data(arr)
        tx.set_text(f'time step = {timestep[i]}')

    ani = animation.FuncAnimation(fig, animate, frames=len(timestep))
    ani.save(figure_path+".gif",writer='pillow',fps=math.floor(len(timestep)/5))

def dump2energy(dump_path, num_steps=None):
    # Create site Hamiltonian energy figure from dump file with energy list

    with open(dump_path+".dump") as file:
        box_size = np.zeros(3)
        for i, line in enumerate(file):
            if i==3: num_sites = int(line)
            if i==5: box_size[0] = np.array(line.split(), dtype=float)[-1]
            if i==6: box_size[1] = np.array(line.split(), dtype=float)[-1]
            if i==7: box_size[2] = np.array(line.split(), dtype=float)[-1]
            if i==8: name_vars = line.split()[2:]
            if i>8: break
    box_size = np.ceil(box_size).astype(int) #reformat box_size
    entry_length = num_sites+9 #there are 9 header lines in each entry

    # total lines for dump
    if num_steps!=None: total_lines = num_steps*entry_length
    else: total_lines=None

    time_steps=[]
    energy_figure=[]
    with open(dump_path+".dump") as file:
        for i, line in tqdm(enumerate(file), "EXTRACTING SPPARKS DUMP (%s.dump)"%dump_path, total=total_lines):
            [entry_num, line_num] = np.divmod(i,entry_length) #what entry number and entry line number does this line number indicate
            if line_num==0: entry = np.zeros(box_size) #set the energy figure matrix
            if line_num==1: time_steps.append(int(float(line.split()[-1]))) #log the time step
            atom_num = line_num-9 #track which atom line we're on
            if atom_num>=0 and atom_num<num_sites:
                line_split = np.array(line.split(), dtype=float)
                site_x = int(line_split[name_vars.index('x')])
                site_y = int(line_split[name_vars.index('y')])
                site_z = int(line_split[name_vars.index('z')])
                entry[site_x,site_y,site_z] = line_split[name_vars.index('energy')] #record valid atom lines
            if line_num==entry_length-1:
                energy_figure.append(entry)
    energy_figure = np.array(energy_figure)
    time_steps = np.array(time_steps)

    return time_steps, energy_figure

if __name__ == '__main__':

    path  = "/Users/lin.yang/projects/SPPARKS-AGG/examples/Test_SimplifyIncE/2d_circle/"
    file_name = "spparks"
    figure_path="video/sk0_000_a00_enrgy_2low_bothB"

    timestep, energy_figure = dump2energy(path+file_name, 81)

    # plot_energy_figure(0, energy_figure)
    plot_energy_video(timestep, energy_figure, figure_path)
