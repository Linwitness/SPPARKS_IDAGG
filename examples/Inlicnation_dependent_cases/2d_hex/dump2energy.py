# -*- coding: utf-8 -*-


from turtle import color
import numpy as np
import os
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm
import sys
current_path = '/Users/lin.yang/Documents/SPYDER/SmoothingAlgos/'
sys.path.append(current_path)
import PACKAGE_MP_Linear as smooth
import myInput

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

def plot_energy_video(timestep, energy_figure, figure_path, delta = 0):

    imgs = []
    fig, ax = plt.subplots()

    cv0 = np.squeeze(energy_figure[0])
    cv0 = np.rot90(cv0,1)
    if delta == 0: colormap_map = np.max(cv0)
    else: colormap_map = (1 + delta) * 8
    im = ax.imshow(cv0,vmin=np.min(cv0),vmax=colormap_map,cmap='Accent',interpolation='none')
    cb = fig.colorbar(im)
    tx = ax.set_title(f'time step = {timestep[0]}')
    # plt.show()

    def animate(i):
        arr=np.squeeze(energy_figure[i])
        arr=np.rot90(arr,1)
        im.set_data(arr)
        tx.set_text(f'time step = {timestep[i]}')

    ani = animation.FuncAnimation(fig, animate, frames=len(timestep))
    FFMpegWriter = animation.writers['ffmpeg']
    writer = animation.FFMpegWriter(fps=math.floor(len(timestep)/5), bitrate=10000)
    ani.save(figure_path+".mp4",writer=writer)

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

def dump2npy(dump_path, num_steps=None):
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
    npy_matrix=[]
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
                entry[site_x,site_y,site_z] = line_split[name_vars.index('type')] #record valid atom lines
            if line_num==entry_length-1:
                npy_matrix.append(entry)
    npy_matrix = np.array(npy_matrix)
    time_steps = np.array(time_steps)

    return time_steps, npy_matrix


def dump2img(dump_path, num_steps=None):
    # Create grain structure figure from dump file with site ID
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
    grain_structure_figure=[]
    with open(dump_path+".dump") as file:
        for i, line in tqdm(enumerate(file), "EXTRACTING SPPARKS DUMP (%s.dump)"%dump_path[-20:], total=total_lines):
            [entry_num, line_num] = np.divmod(i,entry_length) #what entry number and entry line number does this line number indicate
            if line_num==0: entry = np.zeros(box_size) #set the energy figure matrix
            if line_num==1: time_steps.append(int(float(line.split()[-1]))) #log the time step
            atom_num = line_num-9 #track which atom line we're on
            if atom_num>=0 and atom_num<num_sites:
                line_split = np.array(line.split(), dtype=float)
                site_x = int(line_split[name_vars.index('x')])
                site_y = int(line_split[name_vars.index('y')])
                site_z = int(line_split[name_vars.index('z')])
                entry[site_x,site_y,site_z] = line_split[name_vars.index('type')] #record valid atom lines
            if line_num==entry_length-1:
                grain_structure_figure.append(entry)
    grain_structure_figure = np.array(grain_structure_figure)
    time_steps = np.array(time_steps)

    return time_steps, grain_structure_figure

def get_normal_vector(grain_structure_figure_one, grain_num):
    nx = grain_structure_figure_one.shape[0]
    ny = grain_structure_figure_one.shape[1]
    ng = np.max(grain_structure_figure_one)
    cores = 1
    loop_times = 5
    P0 = grain_structure_figure_one
    R = np.zeros((nx,ny,2))
    smooth_class = smooth.linear_class(nx,ny,ng,cores,loop_times,P0,R)

    smooth_class.linear_main("inclination")
    P = smooth_class.get_P()
    # sites = smooth_class.get_gb_list(1)
    # print(len(sites))
    # for id in range(2,grain_num+1): sites += smooth_class.get_gb_list(id)
    # print(len(sites))
    sites = smooth_class.get_all_gb_list()
    sites_together = []
    for id in range(len(sites)): sites_together += sites[id]
    print(f"Total num of GB sites: {len(sites_together)}")

    return P, sites_together

def plot_normal_vector_distribution(P, sites, step, figure_path):
    xLim = [0, 90]
    binValue = 5.01
    binNum = round((abs(xLim[0])+abs(xLim[1]))/binValue)
    xCor = np.linspace((xLim[0]+binValue/2),(xLim[1]-binValue/2),binNum)

    freqArray = np.zeros(binNum)
    degree = []
    # for graini in range(len(P_list)):
    #     P = P_list[graini]
    #     sites = sites_list[graini]
    for sitei in sites:
        [i,j] = sitei
        dx,dy = myInput.get_grad(P,i,j)

        # if dy > 0 and dx >= 0:
        #     degree.append(math.atan(dx/dy))
        # elif dy > 0 and dx < 0:
        #     degree.append(math.atan(dx/dy) + 2*math.pi)
        # elif dy == 0 and dx > 0:
        #     degree.append(math.pi/2)
        # elif dy == 0 and dx < 0:
        #     degree.append(math.pi/2*3)
        # elif dy < 0:
        #     degree.append(math.pi - math.atan(-dx/dy))
        if dx == 0:
            degree.append(math.pi/2)
        elif dy >= 0:
            degree.append(abs(math.atan(-dy/dx)))
        elif dy < 0:
            degree.append(abs(math.atan(dy/dx)))
    for i in range(len(degree)):
        if (degree[i] < 0) or (degree[i] > math.pi/2):
            print("Why??")
        freqArray[int((degree[i]/math.pi*180-xLim[0])/binValue)] += 1
    freqArray = freqArray/sum(freqArray)

    plt.clf()
    fig = plt.subplots()
    # plt.bar(xCor,freqArray,width=binValue*0.7)
    plt.plot(xCor, freqArray,'-o',linewidth=2,label='Distribution')
    plt.xlabel("degree")
    plt.ylabel("frequence")
    plt.ylim([0, 0.1])
    # fitting
    fit_coeff = np.polyfit(xCor, freqArray, 1)
    plt.plot(xCor, xCor*fit_coeff[0]+fit_coeff[1],'--',color='k',linewidth=2,label='fitting')
    plt.legend()
    plt.title(f"The slope of fitting line is {round(fit_coeff[0]*1e4,2)}e4")

    if step < 10:
        step_str = '000' + str(step)
    elif step < 100:
        step_str = '00' + str(step)
    elif step < 1000:
        step_str = '0' + str(step)

    plt.savefig(f'results/{figure_path}_step.{step_str}.png',dpi=400,bbox_inches='tight')


if __name__ == '__main__':
    path  = "/Users/lin.yang/projects/SPPARKS-AGG/examples/Test_SimplifyIncE/2d_hex/"
    file_name = [
                 "p_ori_ave_aveE_512_multiCore4_delta0.6_m2_J1_refer_1_0_0_seed56689_kt066"]
    figure_path = file_name
    num_grain = 512
    output_name = ["ave_0_6_c8_distribution"]

    for i in range(len(file_name)):
        if not os.path.exists("results/"+figure_path[i]+".npy"):
            timestep, grain_structure_figure = dump2img(path+file_name[i], 81)
            np.save("results/"+figure_path[i],grain_structure_figure)
        else:
            grain_structure_figure = np.load(path + "results/"+figure_path[i]+".npy")
            timestep = np.zeros(len(grain_structure_figure))

        if not os.path.exists("results/"+figure_path[i]): os.mkdir("results/"+figure_path[i])
        for step in range(len(timestep)):
            newplace = np.rot90(grain_structure_figure[step,:,:,:], 1, (0,1))
            # try:
            P, sites = get_normal_vector(newplace, num_grain)
            # except ZeroDivisionError:
            #     continue
            if len(sites) == 0: continue

            plot_normal_vector_distribution(P, sites, step, figure_path[i]+"/"+figure_path[i])

        os.system(f'ffmpeg -framerate 30 -i results/{figure_path[i]}/{figure_path[i]}_step.%04d.png \
                    -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p \
                    results/{figure_path[i]}/{output_name[i]}.mp4')

    # path  = "/Users/lin.yang/projects/SPPARKS-AGG/examples/Test_SimplifyIncE/2d_hex/"
    # file_name = "spparks"
    # figure_path="video/h_random_enrgy_2low"

    # timestep, npy_matrix = dump2npy(path+file_name, 101)
    # print(npy_matrix.shape)
    # np.save('hex_inc0.2_xy', npy_matrix)
