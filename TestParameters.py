#Repeatedly runs the simulation with different parameters
#And saves last frames (if generated) to images for comparison

import os
import subprocess
import json
import numpy as np
import matplotlib.pyplot as plt

def CreateIfNotExists(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def Plot2D(data, filepath, title, min, max):
    plt.imshow(data, cmap='plasma', interpolation='bicubic', origin='lower', vmin=min, vmax=max)
    plt.colorbar(label=title)
    plt.savefig(filepath)
    plt.close()

density_dir = "test/density"
velx_dir = "test/velx"
vely_dir = "test/vely"

last_step_output = "output/49.txt"

simconfig = "./simconfigs/cylinder.json"

CreateIfNotExists(density_dir)
CreateIfNotExists(velx_dir)
CreateIfNotExists(vely_dir)

#Define parameter ranges to test:

length_vals = np.arange(0.1, 1.0, 0.1)
time_vals = np.arange(0.1, 1.0, 0.1)
gravity_vals = np.arange(0.1, 1.0, 0.1)

for length in length_vals:
    for time in time_vals:
        for gravity in gravity_vals:
            #Override values in json simulation config
            with open(simconfig, 'r+') as f:
                data = json.load(f)
                data["LengthUnit"] = length 
                data["TimeStep"]   = time 
                data["Gravity"]    = gravity 

                f.seek(0)
                json.dump(data, f, indent=4)
                f.truncate()

            #Call simulation
            subprocess.call("./build/bin/Release/LatticeBoltzman.exe" + " " + "./simconfigs/cylinder.json")

            #Retrieve data from last iteration (if present)
            if not os.path.isfile(last_step_output):
                continue
            
            data = np.loadtxt(last_step_output)

            #Usual drawing plots
            densities    = data[::, 0::3].T
            velocities_x = data[::, 1::3].T
            velocities_y = data[::, 2::3].T

            filename = str(length) + "_" + str(time) + "_" + str(gravity)

            density_path = density_dir + "/" + filename + '.png'
            velx_path    = velx_dir    + "/" + filename + '.png'
            vely_path    = vely_dir    + "/" + filename + '.png'

            Plot2D(densities, density_path, "density", 0.0, 1.0)
            Plot2D(velocities_x, velx_path, "velocity x", -1.0, 1.0)
            Plot2D(velocities_y, vely_path, "velocity y", -1.0, 1.0)

            os.remove(last_step_output)

