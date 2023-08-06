import os
import numpy as np
import matplotlib.pyplot as plt

PLOT_MAX_DENSITY = 1.0
PLOT_MAX_VELOCITY = 1.0

def CreateIfNotExists(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

data_dir = 'output'

density_dir = 'plots/density'
velx_dir = 'plots/velx'
vely_dir = 'plots/vely'

CreateIfNotExists(density_dir)
CreateIfNotExists(velx_dir)
CreateIfNotExists(vely_dir)

def Plot2D(data, filepath, title, min, max):
    plt.imshow(data, cmap='plasma', interpolation='bicubic', origin='lower', vmin=min, vmax=max)
    plt.colorbar(label=title)
    plt.savefig(filepath)
    plt.close()

for filename in os.listdir(data_dir):
    filepath_in = os.path.join(data_dir, filename)

    if not os.path.isfile(filepath_in): 
        continue

    data = np.loadtxt(filepath_in)

    densities    = data[::, 0::3].T
    velocities_x = data[::, 1::3].T
    velocities_y = data[::, 2::3].T

    density_path = density_dir + '/' + os.path.splitext(filename)[0] + '.png'
    velx_path    = velx_dir    + '/' + os.path.splitext(filename)[0] + '.png'
    vely_path    = vely_dir    + '/' + os.path.splitext(filename)[0] + '.png'

    Plot2D(densities, density_path, "density", 0.0, PLOT_MAX_DENSITY)
    Plot2D(velocities_x, velx_path, "velocity x", -PLOT_MAX_VELOCITY, PLOT_MAX_VELOCITY)
    Plot2D(velocities_y, vely_path, "velocity y", -PLOT_MAX_VELOCITY, PLOT_MAX_VELOCITY)