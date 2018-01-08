
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from itertools import product, combinations
import sys



file_object = open(sys.argv[1],"r")
inp = file_object.readlines()
# length of simulation box
length = float(inp[0])
a = inp[1:]

positions = []
complete_positions = []
count = 0
for i in range(len(a)):
    if len(a[i].split()) < 5:
        positions_in_time_step = np.array(positions)
        complete_positions.append(np.array(positions_in_time_step))
        count += 1
        positions = []
    else:
        positions.append([float(x) for x in a[i].split()])



fig = plt.figure()
ax = fig.gca()
ax.set_facecolor('black')
for i in range(0, len(complete_positions)):
    ax.set_xlim(-length / 2, length / 2)
    ax.set_ylim(-length / 2, length / 2)
    plt.pause(0.1)
    plt.cla()
    # print(complete_positions[i][:,2:5])
    p1 = ax.scatter(complete_positions[i][:,0], complete_positions[i][:,1], c=complete_positions[i][:,2:5], marker='o', s=0.5)
    # for j in range(0, len(complete_positions[i])):
    #     [x, y, r, g, b] = complete_positions[i][j]
    #     ax.plot(x, y, c='w', marker='o', markersize=1)

plt.show()
