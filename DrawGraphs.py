import os
import numpy as np
import matplotlib.pyplot as plt

data_dir = 'output'
plots_dir = 'plots'

if not os.path.exists(plots_dir):
   os.makedirs(plots_dir)

for filename in os.listdir(data_dir):
    filepath_in = os.path.join(data_dir, filename)

    if not os.path.isfile(filepath_in): 
        continue

    data = np.loadtxt(filepath_in)

    densities    = data[::, 0::3]
    velocities_x = data[::, 1::3]
    velocities_y = data[::, 2::3]

    filepath_out = plots_dir + '/' + os.path.splitext(filename)[0] + '.png'

    plt.imshow(densities, cmap='plasma', interpolation='bicubic')
    plt.colorbar()
    plt.savefig(filepath_out)
    plt.close()