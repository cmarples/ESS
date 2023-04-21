# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 14:06:26 2023

@author: Callum Marples

Plot binned random point distributions on the unit sphere.

This is Figure 2 of the Geodesic Paper.
"""

import csv
import math
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from figure_size import set_size

save_on = True

# Plot size
plt.style.use('tex')
fig_width = 345.0
#fig_width = 500.0
fs = set_size(fig_width)

fig, axes = plt.subplots(1, 4, figsize=set_size(fig_width, height=2.0, subplots=(1, 4)))
ax = axes.flat

# Grid information
n_th = 181
n_ph = 360
n_patches = (n_th - 2)*n_ph + 2

a = 1.0 # Sphere radius

### Analysis routines
# Convert a 1D array to a 2D array in theta and phi.
def fold_array(C):
    D = np.zeros([n_th, n_ph+1])
    count = 0
    for i in range(n_th):
        if i == 0 or i == n_th-1:
            D[i, :] = C[count]
            count += 1
        else:
            for j in range(n_ph+1):
                if j == n_ph:
                    D[i, j] = D[i, 0]
                else:
                    D[i, j] = C[count]
                    count += 1
    return D

### Read data
# Read patch areas.
patch_areas = np.zeros(n_patches)
area_file = "../data/patch_areas.csv"
with open(area_file) as f:
    i = -2
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if i >= 0:
            patch_areas[i] = row[0] # Sphere row
        i += 1

# Read data for Marsaglia's method, with Lagged Fibonacci.
sphere_contacts = np.zeros([n_patches])
dist_file = "../data/distributions_laggedfib.csv"

# Read distribution file for Lagged Fibonacci
with open(dist_file) as f:
    i = -4
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if i >= 0:
            sphere_contacts[i] = row[1] # Marsaglia's method with Lagged Fibonacci.
        elif i == -3:
            N_str = row[0]
            N = N_str.split(" ")
            N = int(N[2])        
        i += 1

# Determine contact probability per unit area.
min_contacts = np.min(sphere_contacts)
max_contacts = np.max(sphere_contacts)
map_sphere = sphere_contacts / patch_areas
map_sphere /= np.sum(map_sphere)
min_map = np.min(map_sphere)
max_map = np.max(map_sphere)

# Convert to 2D arrays
map_sphere = fold_array(map_sphere)
map_sphere_contacts = fold_array(sphere_contacts)

### Plotting functions
def plot_shape(M, a, b, c, map_name, v, ax_in, norm_in=[0,0]):
    
    u, v = np.mgrid[0:np.pi:181j, 0:2*np.pi:360j]
    #v -= 10*np.pi/180
    
    norm=mpl.colors.Normalize(vmin = norm_in[0],
                              vmax = norm_in[1], clip = False)
    
    x = a * np.sin(u) * np.cos(v)
    y = b * np.sin(u) * np.sin(v)
    z = c * np.cos(u)
    
    ax_in.plot_surface(x, y, z, rstride=1, cstride=1, cmap=plt.get_cmap(map_name),
                       linewidth=0, linestyle='None', antialiased=True,
                       facecolors=plt.get_cmap(map_name)(norm(M)),
                       vmin = v[0], vmax = v[1])
    
    ax_in.axis('equal')

    ax_in.axis('off')
    ax_in.view_init(elev=35, azim=20)
    return norm

### Draw figures
vc = [0, 2600]
vp = [8.4e-6, 2.7e-5]

ax[0].remove()
ax[0] = fig.add_subplot(1, 2, 1, projection='3d')
norm1 = plot_shape(map_sphere_contacts, a, a, a, 'inferno', vc, ax[0], vc)
ax[0].set_title('(a)', x=0.54, y=-0.15, fontsize=10)

ax[2].remove()
ax[2] = fig.add_subplot(1, 2, 2, projection='3d')
norm2 = plot_shape(map_sphere, a, a, a, 'coolwarm', vc, ax[2], vp)
ax[2].set_title('(b)', x=0.54, y=-0.15, fontsize=10)

### Include colourbars
ax[1].remove()
ax[1] = fig.add_axes([0.46, 0.21, 0.01, 0.6])
newcmp = mpl.cm.inferno
norm = mpl.colors.Normalize(vmin=0, vmax=2600)
cb1 = mpl.colorbar.ColorbarBase(ax[1], cmap=newcmp, norm=norm1, orientation='vertical',
                                ticks=[0, 500, 1000, 1500, 2000, 2500])

ax[3].remove()
ax[3] = fig.add_axes([0.88, 0.21, 0.01, 0.6])
newcmp = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=0.0000085, vmax=0.000025)
cb1 = mpl.colorbar.ColorbarBase(ax[3], cmap=newcmp, norm=norm2, orientation='vertical',
                                ticks=[0.00001, 0.000015, 0.00002, 0.000025])


plt.show()

if save_on:
    fig_name = "Img/im_sphere_distributions.pdf"
    fig.savefig(fig_name, format='pdf', bbox_inches='tight')
    