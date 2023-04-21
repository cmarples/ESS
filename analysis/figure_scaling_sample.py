"""
Created on 10/12/21
Last Modified on 07/02/23

Show a full set of sampled points on the ellipsoid surface, generated by naive
scaling and by gradient rejection.

Figure used to illustrate how naive scaling results in a 
non-uniform distribution.

NOT USED IN PAPER!

@author: Callum Marples
"""

import csv
import math
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
from figure_size import set_size

save_on = True

fig_width = 345.0
fig, axes = plt.subplots(1, 3, figsize=set_size(fig_width, height=1.0, subplots=(1, 3)))
ax = axes.flat

N = 1000

def read_file(file_name, N):
    with open(file_name) as f:
        XYZ = np.zeros([N, 3])
        row_no = 0
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            XYZ[row_no] = row
            row_no += 1
        return XYZ

def plot_points(shape, N, xyz, axes):
    if shape == 'sphere':
        a = 1.0
        b = 1.0
        c = 1.0
    elif shape == 'ellipsoid':
        a = 3.0
        b = 3.0
        c = 1.5
    else:
        raise Exception('SHAPE NOT GIVEN')
    #fig = plt.figure(figsize=(fs[0], fs[0]))
    #ax = fig.add_subplot(111, projection='3d')
    # Plot points
    for i in range(N):
        axes.scatter(xyz[i, 0], xyz[i, 1], xyz[i, 2], marker=".", c="k", s=0.05)
    # Set axis limits - make all of them the same to get same scaling
    m = max([a, b, c]) + 0.1*c
    #axes.set_xlim3d(-m, m)
    #axes.set_ylim3d(-m, m)
    #axes.set_zlim3d(-m, m)
    # Set viewing angle
    axes.view_init(elev=20, azim=0)
    # Remove axes
    axes.set_axis_off()
    #plt.show()
    return True
    
ax[0].remove()
ax[0] = fig.add_subplot(1, 3, 1, projection='3d')    
sphere_file = "../data/sample_sphere_" + str(N) + ".csv"
xyz_sphere = read_file(sphere_file, N)
plot_points('sphere', N, xyz_sphere, ax[0])

ax[1].remove()
ax[1] = fig.add_subplot(1, 3, 2, projection='3d')
naive_file = "../data/sample_naive_" + str(N) + ".csv"
xyz_naive = read_file(naive_file, N)
plot_points('ellipsoid', N, xyz_naive, ax[1])

ax[2].remove()
ax[2] = fig.add_subplot(1, 3, 3, projection='3d')
ellipsoid_file = "../data/sample_correct_" + str(N) + ".csv"
xyz_correct = read_file(ellipsoid_file, N)
plot_points('ellipsoid', N, xyz_correct, ax[2])

if save_on:
    fig_name = "Img/im_sample_points.pdf"
    fig.savefig(fig_name, format='pdf', bbox_inches='tight')

