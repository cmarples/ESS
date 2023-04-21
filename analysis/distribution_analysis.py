"""
Created on 29/11/21
Last Modified on 06/02/23

@author: Callum Marples

Analyse distribution data from ESS and compute:
    - relative standard deviation
    - chi squared
"""

import csv
import numpy as np
from scipy import stats

#n_th = 5
#n_ph = 4

n_th = 181
n_ph = 360
n_patches = (n_th - 2)*n_ph + 2

def fold_array(C):
    D = np.zeros([n_th, n_ph])
    count = 0
    for i in range(n_th):
        if i == 0 or i == n_th-1:
            D[i, :] = C[count]
            count += 1
        else:
            for j in range(n_ph):
                D[i, j] = C[count]
                count += 1
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

distribution_file = ["../data/distributions_laggedfib.csv", "../data/distributions_yarn.csv"]
rng_name = ["laggedfib", "yarn"]
for rng_index in range(2):

    sphere_contacts = np.zeros([6, n_patches])
    oblate_contacts = np.zeros([8, n_patches])
    prolate_contacts = np.zeros([8, n_patches])
    triaxial_contacts = np.zeros([8, n_patches])
    
    with open(distribution_file[rng_index]) as f:
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
    
    # chi-squared tests for uniformity
    sphere_chi2 = np.zeros(6)
    oblate_chi2 = np.zeros(8)
    prolate_chi2 = np.zeros(8)
    triaxial_chi2 = np.zeros(8)
    p_values_sphere = np.zeros(6)
    p_values_oblate = np.zeros(8)
    p_values_prolate = np.zeros(8)
    p_values_triaxial = np.zeros(8)
    
    # mean and standard deviation
    sphere_mean = np.zeros(6)
    sphere_sigma = np.zeros(6)
    prolate_mean = np.zeros(8)
    prolate_sigma = np.zeros(8)
    oblate_mean = np.zeros(8)
    oblate_sigma = np.zeros(8)
    triaxial_mean = np.zeros(8)
    triaxial_sigma = np.zeros(8)
    
    alpha = 0.05 # 95% confidence level.
    nu = n_patches
    crit_val = stats.chi2.ppf(1.0-alpha, nu)
    for i in range(6):
        # Spheres
        expected_spheres = N * patch_areas[0] / sum(patch_areas[0])   # Expected contacts
        sphere_chi2[i] = sum( (sphere_contacts[i] - expected_spheres)**2 / expected_spheres )
        p_values_sphere[i] = 1.0 - stats.chi2.cdf(sphere_chi2[i], nu)
        
        prob = sphere_contacts[i] / patch_areas[0]
        sphere_mean[i] = np.mean(prob)
        sphere_sigma[i] = np.std(prob)
        
    for i in range(8):    
        # Oblate 
        expected_oblate = N * patch_areas[1] / sum(patch_areas[1])  # Expected contacts
        oblate_chi2[i] = sum( (oblate_contacts[i] - expected_oblate)**2 / expected_oblate )
        p_values_oblate[i] = 1.0 - stats.chi2.cdf(oblate_chi2[i], nu)
        # Spheroids
        expected_prolate = N * patch_areas[2] / sum(patch_areas[2])  # Expected contacts
        prolate_chi2[i] = sum( (prolate_contacts[i] - expected_prolate)**2 / expected_prolate )
        p_values_prolate[i] = 1.0 - stats.chi2.cdf(prolate_chi2[i], nu)
        # Ellipsoids
        expected_triaxial = N * patch_areas[3] / sum(patch_areas[3]) # Expected contacts
        triaxial_chi2[i] = sum( (triaxial_contacts[i] - expected_triaxial)**2 / expected_triaxial )
        p_values_triaxial[i] = 1.0 - stats.chi2.cdf(triaxial_chi2[i], nu)
        
        prob = oblate_contacts[i] / patch_areas[1]
        oblate_mean[i] = np.mean(prob)
        oblate_sigma[i] = np.std(prob)
        
        prob = prolate_contacts[i] / patch_areas[2]
        prolate_mean[i] = np.mean(prob)
        prolate_sigma[i] = np.std(prob)
        
        prob = triaxial_contacts[i] / patch_areas[3]
        triaxial_mean[i] = np.mean(prob)
        triaxial_sigma[i] = np.std(prob)
    
    sphere_rsd = sphere_sigma / sphere_mean
    oblate_rsd = oblate_sigma / oblate_mean
    prolate_rsd = prolate_sigma / prolate_mean
    triaxial_rsd = triaxial_sigma / triaxial_mean
    
    ### Write to .csv files
    chi_file = "../data/chi2_tests_" + rng_name[rng_index] + ".csv"  
    with open(chi_file, 'w') as writer:
        writer.write('chi Squared Statistics\n')
        writer.write('Spheres\n')
        writer.write('Cubic Rejection,Marsaglia,Trig,Gaussian,Cook,Area Rejection\n')
        for i in sphere_chi2:
            writer.write(str(i))
            writer.write(',')
        writer.write('\n')
        writer.write('Ellipsoids\n')
        writer.write('Naive,Gradient Reject,Grad Rej (Trig),Grad Rej (Pols),Area Rejection,Area Rej (Pols),Area Rej (Mer),Generic\n')
        writer.write('Oblate\n')
        for i in oblate_chi2:
            writer.write(str(i))
            writer.write(',')
        writer.write('\n')
        writer.write('Prolate\n')
        for i in prolate_chi2:
            writer.write(str(i))
            writer.write(',')
        writer.write('\n')
        writer.write('Triaxial\n')
        for i in triaxial_chi2:
            writer.write(str(i))
            writer.write(',')
        writer.write('\n')
        writer.write('Critical Value\n')
        writer.write(str(crit_val))
    
    rsd_file = "../data/rsd_" + rng_name[rng_index] + ".csv"  
    with open(rsd_file, 'w') as writer:
        writer.write('Relative standard deviation\n')
        writer.write('Spheres\n')
        writer.write('Cubic Rejection,Marsaglia,Trig,Gaussian,Cook,Area Rejection\n')
        for i in sphere_rsd:
            writer.write(str(i))
            writer.write(',')
        writer.write('\n')
        writer.write('Ellipsoids\n')
        writer.write('Naive,Gradient Reject,Grad Rej (Trig),Grad Rej (Pols),Area Rejection,Area Rej (Pols),Area Rej (Mer),Generic\n')
        writer.write('Oblate\n')
        for i in oblate_rsd:
            writer.write(str(i))
            writer.write(',')
        writer.write('\n')
        writer.write('Prolate\n')
        for i in prolate_rsd:
            writer.write(str(i))
            writer.write(',')
        writer.write('\n')
        writer.write('Triaxial\n')
        for i in triaxial_rsd:
            writer.write(str(i))
            writer.write(',')
