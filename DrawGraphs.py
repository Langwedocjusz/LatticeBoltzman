import os
import numpy as np
import matplotlib.pyplot as plt

data_dir = 'output'
plots_dir = 'plots'

for filename in os.listdir(data_dir):
    filepath = os.path.join(data_dir, filename)

    if not os.path.isfile(filepath): 
        continue

    data = np.loadtxt(filepath)

    densities    = data[::, 0::3]
    velocities_x = data[::, 1::3]
    velocities_y = data[::, 2::3]

    #plt.imshow(data, cmap='hot', interpolation='nearest')
    #plt.show()