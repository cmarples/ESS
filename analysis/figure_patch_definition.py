"""
Created on 09/09/21
Last Modified on 07/02/23

Plot a patch on the surface of a triaxial ellipsoid.

Figure used in Introduction to explain definition of
a 'patch'.

This is Figure 1 of the Sampling Paper.

@author: Callum Marples
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
import math
from figure_size import set_size

save_on = True

# Ellipsoid Axes
a = 3.0
b = 2.0
c = 1.0

# Plot size
plt.style.use('tex')
width = 345.0
fs = set_size(width)

# Plot ellipsoid
fig = plt.figure(figsize=(fs[0], fs[0]))
ax = fig.add_subplot(111, projection='3d')
th = np.linspace(0, math.pi, 100)
ph = np.linspace(0, 2*math.pi, 100)
TH, PH = np.meshgrid(th, ph)
X = a * np.sin(TH) * np.cos(PH)
Y = b * np.sin(TH) * np.sin(PH)
Z = c * np.cos(TH)
ax.plot_surface(X, Y, Z, color=(0.3, 0.3, 0.3), linewidth=0, alpha=0.1, shade=False)

# Plot patch
th3 = 30*math.pi/180
th4 = 60*math.pi/180
th2 = np.linspace(th3, th4, 36)
ph3 = 240*math.pi/180
ph4 = 300*math.pi/180
ph2 = np.linspace(ph3, ph4, 36)
TH2, PH2 = np.meshgrid(th2, ph2)
X2 = (a+0.0) * np.sin(TH2) * np.cos(PH2)
Y2 = (b+0.0) * np.sin(TH2) * np.sin(PH2)
Z2 = (c+0.0) * np.cos(TH2)
ax.plot_surface(X2, Y2, Z2, color='r', linewidth=0, alpha=1, shade=False)

# Plot grid lines
X3 = a * math.sin(th3) * np.cos(ph)
Y3 = b * math.sin(th3) * np.sin(ph)
Z3 = np.array(100*[1]) * c * math.cos(th3)
ax.plot(X3, Y3, Z3, 'k')
X3 = a * math.sin(th4) * np.cos(ph)
Y3 = b * math.sin(th4) * np.sin(ph)
Z3 = np.array(100*[1]) * c * math.cos(th4)
ax.plot(X3, Y3, Z3, 'k')
thx = np.linspace(0, 2*math.pi, 100) 
X3 = a * np.sin(th) * math.cos(ph3)
Y3 = b * np.sin(th) * math.sin(ph3)
Z3 = c * np.cos(th)
ax.plot(X3, Y3, Z3, 'k', linestyle='dashed')
X3 = a * np.sin(th) * math.cos(ph4)
Y3 = b * np.sin(th) * math.sin(ph4)
Z3 = c * np.cos(th)
ax.plot(X3, Y3, Z3, 'k', linestyle='dashed')

# Plot theta axis
Xp = 2 * [0]
Yp = 2 * [0]
Zp = [-3*c, 3*c]
ax.plot(Xp, Yp, Zp, color = 'k', linestyle='dotted')

# Add text
ax.text(a*math.sin(th3)*math.cos(ph3-0.3), b*math.sin(th3)*math.sin(ph3-0.3), c*math.cos(th3)+0.1, r"$\theta_0$")
ax.text(a*math.sin(th4)*math.cos(ph3-0.3), b*math.sin(th4)*math.sin(ph3-0.3), c*math.cos(th4)+0.1, r"$\theta_1$")
ax.text(a*math.sin(th4+0.6)*math.cos(ph3-0.15), b*math.sin(th4+0.6)*math.sin(ph3-0.15), c*math.cos(th4+0.6), r"$\phi_0$")
ax.text(a*math.sin(th4+0.6)*math.cos(ph4+0.05), b*math.sin(th4+0.6)*math.sin(ph4+0.05), c*math.cos(th4+0.6)+0.1, r"$\phi_1$")


# Set axis limits - make all of them the same to get same scaling
m = 1.95
ax.set_xlim3d(-m, m)
ax.set_ylim3d(-m, m)
ax.set_zlim3d(-m, m)

# Set viewing angle
ax.view_init(elev=34, azim=-90)
    
# Remove axes
plt.axis('off')
plt.show()    

# Note: will need to crop the image after it has been made for it to fit
# nicely in a LaTeX generated document.
if save_on:
    fig_name = "Img/im_patch_definition_uncropped.pdf"
    fig.savefig(fig_name, format='pdf', bbox_inches='tight', pad_inches=0)