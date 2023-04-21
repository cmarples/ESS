# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 18:17:41 2023

@author: Callum Marples

Plot binned random point distributions on the ellipsoid.

This is Figure 3 of the Geodesic Paper.
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
fs = set_size(fig_width)

fig, axes = plt.subplots(2, 4, figsize=set_size(fig_width, height=2.0, subplots=(2, 4)))
ax = axes.flat

# Grid information
n_th = 181
n_ph = 360
n_patches = (n_th - 2)*n_ph + 2

a = 3.0 
b = 2.0 # Ellipsoid axes
c = 1.0

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
# Read ellipsoid patch areas.
patch_areas = np.zeros(n_patches)
area_file = "../data/patch_areas.csv"
with open(area_file) as f:
    i = -2
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if i >= 0:
            patch_areas[i] = row[3] # Triaxial row
        i += 1

# Read data for the naive method on a triaxial ellipsoid.
con_naive = np.zeros([n_patches])
con_grad  = np.zeros([n_patches])
dist_file = "../data/distributions_laggedfib.csv"

# Read distribution file for Lagged Fibonacci
with open(dist_file) as f:
    i = -4
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if i >= 0:
            con_naive[i] = row[22] # Naive
            con_grad[i]  = row[23] # Gradient rejection
        elif i == -3:
            N_str = row[0]
            N = N_str.split(" ")
            N = int(N[2])        
        i += 1
        
# Determine contact probability per unit area.
min_naive = np.min(con_naive)
max_naive = np.max(con_naive)
map_naive = con_naive / patch_areas
map_naive /= np.sum(map_naive)
min_map_naive = np.min(map_naive)
max_map_naive = np.max(map_naive)

min_grad = np.min(con_grad)
max_grad = np.max(con_grad)
map_grad = con_grad / patch_areas
map_grad /= np.sum(map_grad)
min_map_grad = np.min(map_grad)
max_map_grad = np.max(map_grad)

# Convert to 2D arrays
map_naive = fold_array(map_naive)
map_naive_contacts = fold_array(con_naive)

map_grad = fold_array(map_grad)
map_grad_contacts = fold_array(con_grad)

### Plotting functions
def plot_shape(M, a, b, c, map_name, vm, ax_in, norm_in=[0,0]):
    
    u, v = np.mgrid[0:np.pi:181j, 0:2*np.pi:360j]
    
    norm=mpl.colors.Normalize(vmin = norm_in[0],
                              vmax = norm_in[1], clip = False)
    
    x = a * np.sin(u) * np.cos(v)
    y = b * np.sin(u) * np.sin(v)
    z = c * np.cos(u)
    
    ax_in.plot_surface(x, y, z, rstride=1, cstride=1, cmap=plt.get_cmap(map_name),
                       linewidth=0, linestyle='None', antialiased=True,
                       facecolors=plt.get_cmap(map_name)(norm(M)),
                       vmin = vm[0], vmax = vm[1])
    
    ax_in.axis('equal')
    lim = 2.0
    ax_in.set_xlim(-lim, lim)
    ax_in.set_ylim(-lim, lim)
    ax_in.set_zlim(-lim, lim)
    ax_in.axis('off')
    ax_in.view_init(elev=15, azim=20)
    return norm        
        
### Draw figures
vc = [0, 3000]
vp = [0.0, 3.3e-5]

# Naive
ax[0].remove()
ax[0] = fig.add_subplot(2, 2, 1, projection='3d')
norm1 = plot_shape(map_naive_contacts, a, b, c, 'inferno', vc, ax[0], vc)
ax[0].set_title('(a)', x=0.54, y=-0.05, fontsize=10)

ax[2].remove()
ax[2] = fig.add_subplot(2, 2, 2, projection='3d')
norm2 = plot_shape(map_naive, a, b, c, 'coolwarm', vc, ax[2], vp)
ax[2].set_title('(b)', x=0.54, y=-0.05, fontsize=10)

# Gradient rejection
ax[4].remove()
ax[4] = fig.add_subplot(2, 2, 3, projection='3d')
norm3 = plot_shape(map_grad_contacts, a, b, c, 'inferno', vc, ax[4], vc)
ax[4].set_title('(c)', x=0.54, y=-0.05, fontsize=10)

ax[6].remove()
ax[6] = fig.add_subplot(2, 2, 4, projection='3d')
norm4 = plot_shape(map_grad, a, b, c, 'coolwarm', vc, ax[6], vp)
ax[6].set_title('(d)', x=0.54, y=-0.05, fontsize=10)

### Include colourbars
ax[1].remove()
ax[1] = fig.add_axes([0.46, 0.6, 0.01, 0.25])
newcmp = mpl.cm.inferno
norm = mpl.colors.Normalize(vmin=vc[0], vmax=vc[1])
cb1 = mpl.colorbar.ColorbarBase(ax[1], cmap=newcmp, norm=norm1, orientation='vertical',
                                ticks=[0, 500, 1000, 1500, 2000, 2500, 3000])

ax[3].remove()
ax[3] = fig.add_axes([0.88, 0.6, 0.01, 0.25])
newcmp = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=vp[0], vmax=vp[1])
cb1 = mpl.colorbar.ColorbarBase(ax[3], cmap=newcmp, norm=norm2, orientation='vertical',
                                ticks=[0.0, 0.5e-5, 1.0e-5, 1.5e-5, 2.0e-5, 2.5e-5, 3.0e-5, 3.5e-5])

ax[5].remove()
ax[5] = fig.add_axes([0.46, 0.15, 0.01, 0.25])
newcmp = mpl.cm.inferno
norm = mpl.colors.Normalize(vmin=vc[0], vmax=vc[1])
cb1 = mpl.colorbar.ColorbarBase(ax[5], cmap=newcmp, norm=norm3, orientation='vertical',
                                ticks=[0, 500, 1000, 1500, 2000, 2500, 3000])

ax[7].remove()
ax[7] = fig.add_axes([0.88, 0.15, 0.01, 0.25])
newcmp = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=vp[0], vmax=vp[1])
cb1 = mpl.colorbar.ColorbarBase(ax[7], cmap=newcmp, norm=norm4, orientation='vertical',
                                ticks=[0.0, 0.5e-5, 1.0e-5, 1.5e-5, 2.0e-5, 2.5e-5, 3.0e-5, 3.5e-5])


plt.show()

if save_on:
    fig_name = "Img/im_ellipsoid_distributions.pdf"
    fig.savefig(fig_name, format='pdf', bbox_inches='tight')        
        
        