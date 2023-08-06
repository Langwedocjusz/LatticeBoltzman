import os
import numpy as np
import matplotlib.pyplot as plt

MAX_VELOCITY = 40.0
LATTICE_SIZEX = 128
SLICE = 64

sample_file = "output/49.txt"

positions_x = np.arange(0, LATTICE_SIZEX)

data = np.loadtxt(sample_file)

velocities_y = data[::, 2::3].T

vel_y = velocities_y[SLICE, ::]

plt.ylim([-MAX_VELOCITY, MAX_VELOCITY])
plt.plot(positions_x, vel_y)
plt.show()