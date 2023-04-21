# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 15:58:02 2021

@author: Callum Marples
"""

import csv
import math
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
from my_plot import set_size

save_on = True
shape = 'sphere'
rng = 'laggedfib'

def plot_map(M, map_name):
    fig, ax = plt.subplots(1, 1, figsize=(fs[0], fs[0]))
    im = plt.imshow(M, cmap=plt.get_cmap(map_name)) 
    xt = np.arange(0, 420, 60)
    xx = np.arange(0, 420, 60)
    plt.xticks(xx, xt)
    yt = np.arange(0, 210, 30)
    yy = np.arange(0, 210, 30)
    plt.yticks(yy, yt)
    plt.xlabel('$\phi$')
    plt.ylabel('$\\theta$')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax)

    return fig

def plot_bar(M, map_name, v):
    fig, ax = plt.subplots(1, 1, figsize=(fs[0], fs[0]))
    
    im = plt.imshow(M, cmap=plt.get_cmap(map_name), vmin = v[0], vmax = v[1])
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax)
    ax.remove()
    
    return fig

def plot_bar_discrete(M):
    fig, ax = plt.subplots(1, 1, figsize=(fs[0], fs[0]))
    cmap = colors.ListedColormap([[0,1,1], [0.6,0.6,0.6], [1,0,1]])
    bounds=[-1.5, -0.5, 0.5, 1.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    im = plt.imshow(M, interpolation='nearest', origin='lower', cmap=cmap, norm=norm)
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(im, cax=cax, cmap=cmap, norm=norm, boundaries=bounds, ticks=[-1, 0, 1])
    ax.remove()
    
    return fig

def plot_shape(M, a, b, c, map_name, v):
    fig = plt.figure(figsize=(fs[0], fs[0]))
    ax = fig.add_subplot(111, projection='3d')
    u, v = np.mgrid[0:np.pi:181j, 0:2*np.pi:360j]
    #v -= 10*np.pi/180
    
    strength = M
    norm=colors.Normalize(vmin = np.min(strength),
                          vmax = np.max(strength), clip = False)
    
    x = a * np.sin(u) * np.cos(v)
    y = b * np.sin(u) * np.sin(v)
    z = c * np.cos(u)
    
    ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=plt.get_cmap(map_name),
                           linewidth=0, antialiased=True,
                           facecolors=plt.get_cmap(map_name)(norm(strength)),
                           vmin = v[0], vmax = v[1])
    def axisEqual3D(ax):
        extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
        sz = extents[:,1] - extents[:,0]
        centers = np.mean(extents, axis=1)
        maxsize = max(abs(sz))
        r = maxsize/2
        for ctr, dim in zip(centers, 'xyz'):
            getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
    axisEqual3D(ax)
    plt.axis('off')
    ax.view_init(elev=35, azim=20)
    plt.show()
    return fig

def plot_shape_discrete(M, a, b, c):
    fig = plt.figure(figsize=(fs[0], fs[0]))
    ax = fig.add_subplot(111, projection='3d')
    u, v = np.mgrid[0:np.pi:181j, 0:2*np.pi:360j]
    
    cmap = colors.ListedColormap([[0,1,1], [0.6,0.6,0.6], [1,0,1]])
    bounds=[-1.5, -0.5, 0.5, 1.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    # tell imshow about color map so that only set colors are used
    #img = plt.imshow(prob, interpolation='nearest', origin='lower',
    #                    cmap=cmap, norm=norm)
    
    strength = M
    norm=colors.Normalize(vmin = np.min(strength),
                          vmax = np.max(strength), clip = False)
    
    x = a * np.sin(u) * np.cos(v)
    y = b * np.sin(u) * np.sin(v)
    z = c * np.cos(u)
    
    ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cmap, norm=norm, linewidth=0,
                    facecolors=plt.get_cmap(cmap)(norm(strength)))
    
    def axisEqual3D(ax):
        extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
        sz = extents[:,1] - extents[:,0]
        centers = np.mean(extents, axis=1)
        maxsize = max(abs(sz))
        r = maxsize/2
        for ctr, dim in zip(centers, 'xyz'):
            getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
            
    axisEqual3D(ax)
    plt.axis('off')
    ax.view_init(elev=35, azim=20)
    plt.show()
    return fig  




  
n_th = 181
n_ph = 360
n_patches = (n_th - 2)*n_ph + 2

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

def deviation(C):
    m = np.mean(C)
    s = np.std(C)
    D = np.zeros(n_patches)
    for i in range(n_patches):
        if math.fabs(C[i] - m) > s:
            if C[i] > m:
                D[i] = 1
            else:
                D[i] = -1
        else:
            D[i] = 0     
    D = fold_array(D)
    return D
    
axes = np.zeros([4, 3])
shape_file = "../data/ellipsoid_shapes.csv"
with open(shape_file) as f:
    r_no = 1
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if r_no == 2:    # Sphere axes
            axes[0][0] = row[0]
            axes[0][1] = row[1]
            axes[0][2] = row[2]
        elif r_no == 4: # Oblate axes
            axes[1][0] = row[0]
            axes[1][1] = row[1]
            axes[1][2] = row[2]
        elif r_no == 6: # Prolate axes
            axes[2][0] = row[0]
            axes[2][1] = row[1]
            axes[2][2] = row[2]
        elif r_no == 8: # Ellipsoid axes
            axes[3][0] = row[0]
            axes[3][1] = row[1]
            axes[3][2] = row[2]
        r_no += 1    
            
patch_areas = np.zeros([4, n_patches])
area_file = "../data/patch_areas.csv"
with open(area_file) as f:
    i = -2
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if i >= 0:
            patch_areas[0][i] = row[0]
            patch_areas[1][i] = row[1]
            patch_areas[2][i] = row[2]
            patch_areas[3][i] = row[3]
        i += 1

sphere_contacts = np.zeros([6, n_patches])
oblate_contacts = np.zeros([8, n_patches])
prolate_contacts = np.zeros([8, n_patches])
triaxial_contacts = np.zeros([8, n_patches])
distribution_file = "../data/distributions_" + rng + ".csv"
with open(distribution_file) as f:
    i = -4
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if i >= 0:
            for j in range(6):
                sphere_contacts[j][i] = row[j]
            for j in range(8):
                oblate_contacts[j][i] = row[j+6]
            for j in range(8):
                prolate_contacts[j][i] = row[j+14]
            for j in range(8):
                triaxial_contacts[j][i] = row[j+22]
        elif i == -3:
            N_str = row[0]
            N = N_str.split(" ")
            N = int(N[2])        
        i += 1

plt.style.use('tex')
tex_font = {
    "font.size": 8
    }
plt.rcParams.update(tex_font)
width = 345.0
fs = set_size(width)

if shape == 'sphere':
    map_sphere = sphere_contacts[1] / patch_areas[0]
    map_sphere /= np.sum(map_sphere)
    prob = deviation(map_sphere)
    
    map_sphere = fold_array(map_sphere)
    map_sphere_contacts = fold_array(sphere_contacts[1])
    s1 = "sphere_contacts_Marsaglia_" + rng
    s2 = "sphere_map_Marsaglia_" + rng
    s3 = "sphere_difference_Marsaglia_" + rng
    s4 = "sphere_contacts_bar_Marsaglia_" + rng
    s5 = "sphere_map_bar_Marsaglia_" + rng
    s6 = "sphere_difference_bar_Marsaglia_" + rng
    
    vc = [0, 2600]
    vp = [8.4e-6, 2.7e-5]
    
    fig1 = plot_shape(map_sphere_contacts, axes[0][0], axes[0][1], axes[0][2], 'inferno', vc)
    fig2 = plot_shape(map_sphere, axes[0][0], axes[0][1], axes[0][2], 'cool', vp)
    fig3 = plot_shape_discrete(prob, axes[0][0], axes[0][1], axes[0][2])
    
    if save_on:
        fig4 = plot_bar(map_sphere_contacts, 'inferno', vc)
        fig5 = plot_bar(map_sphere, 'cool', vp)
        fig6 = plot_bar_discrete(prob)
        fig_name = "im_" + s1 + ".png"
        fig1.savefig(fig_name, format='png', bbox_inches='tight', dpi=500)
        fig_name = "im_" + s2 + ".png"
        fig2.savefig(fig_name, format='png', bbox_inches='tight', dpi=500)
        fig_name = "im_" + s3 + ".png"
        fig3.savefig(fig_name, format='png', bbox_inches='tight', dpi=500)
        fig_name = "im_" + s4 + ".png"
        fig4.savefig(fig_name, format='png', bbox_inches='tight', dpi=500)
        fig_name = "im_" + s5 + ".png"
        fig5.savefig(fig_name, format='png', bbox_inches='tight', dpi=500)
        fig_name = "im_" + s6 + ".png"
        fig6.savefig(fig_name, format='png', bbox_inches='tight', dpi=500)
       
if shape == 'ellipsoid':
    map_oblate = oblate_contacts[1] / patch_areas[1]
    map_prolate = prolate_contacts[1] / patch_areas[2]
    map_triaxial = triaxial_contacts[1] / patch_areas[3]
    
    map_oblate /= np.sum(map_oblate)
    map_prolate /= np.sum(map_prolate)
    map_triaxial /= np.sum(map_triaxial)
    
    dev_oblate = deviation(map_oblate)
    dev_prolate = deviation(map_prolate)
    dev_triaxial = deviation(map_triaxial)
    
    map_oblate = fold_array(map_oblate)
    map_prolate = fold_array(map_prolate)
    map_triaxial = fold_array(map_triaxial)
    map_oblate_contacts = fold_array(oblate_contacts[1])
    map_prolate_contacts = fold_array(prolate_contacts[1])
    map_triaxial_contacts = fold_array(triaxial_contacts[1])
    
    s1 = "oblate_contacts_" + rng
    s2 = "oblate_map_" + rng
    s3 = "oblate_difference_" + rng
    
    s4 = "prolate_contacts_" + rng
    s5 = "prolate_map_" + rng
    s6 = "prolate_difference_" + rng
    
    s7 = "triaxial_contacts_" + rng
    s8 = "triaxial_map_" + rng
    s9 = "triaxial_difference_" + rng
    
    s10 = "oblate_contacts_bar"
    s11 = "prolate_contacts_bar"
    s12 = "triaxial_contacts_bar"
    s13 = "oblate_map_bar"
    s14 = "prolate_map_bar"
    s15 = "triaxial_map_bar"
    
    s = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15]
    
    vOc = [np.min(map_oblate_contacts), np.max(map_oblate_contacts)]
    vOm = [np.min(map_oblate), np.max(map_oblate)]
    vPc = [np.min(map_prolate_contacts), np.max(map_prolate_contacts)]
    vPm = [np.min(map_prolate), np.max(map_prolate)]
    vTc = [np.min(map_triaxial_contacts), np.max(map_triaxial_contacts)]
    vTm = [np.min(map_triaxial), np.max(map_triaxial)]
    
    fig1 = plot_shape(map_oblate_contacts, axes[1][0], axes[1][1], axes[1][2], 'inferno', vOc)
    fig2 = plot_shape(map_oblate, axes[1][0], axes[1][1], axes[1][2], 'cool', vOm)
    fig3 = plot_shape_discrete(dev_oblate, axes[1][0], axes[1][1], axes[1][2])
    
    fig4 = plot_shape(map_prolate_contacts, axes[2][0], axes[2][1], axes[2][2], 'inferno', vPc)
    fig5 = plot_shape(map_prolate, axes[2][0], axes[2][1], axes[2][2], 'cool', vPm)
    fig6 = plot_shape_discrete(dev_prolate, axes[2][0], axes[2][1], axes[2][2])
    
    fig7 = plot_shape(map_triaxial_contacts, axes[3][0], axes[3][1], axes[3][2], 'inferno', vTc)
    fig8 = plot_shape(map_triaxial, axes[3][0], axes[3][1], axes[3][2], 'cool', vTm)
    fig9 = plot_shape_discrete(dev_triaxial, axes[3][0], axes[3][1], axes[3][2])
    
    if save_on:
        fig10 = plot_bar(map_oblate_contacts, 'inferno', vOc)
        fig11 = plot_bar(map_prolate_contacts, 'inferno', vPc)
        fig12 = plot_bar(map_triaxial_contacts, 'inferno', vTc)
        fig13 = plot_bar(map_oblate, 'cool', vOm)
        fig14 = plot_bar(map_prolate, 'cool', vPm)
        fig15 = plot_bar(map_triaxial, 'cool', vTm)
        
        f = [fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12, fig13, fig14, fig15]
        
        for i in range(15):
            fig_name = "im_" + s[i] + ".png"
            f[i].savefig(fig_name, format='png', bbox_inches='tight', dpi=500)
         