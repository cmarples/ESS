# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 14:38:53 2023

Reads the data from the various files and writes out .txt files containing
LaTeX code for the tables.

This generates Tables 1, 2 and 3 of the Sampling Paper.

@author: Callum Marples
"""

import numpy as np
import csv

### Read data
# Shapes
axes = np.zeros([4, 3])
file_name = "../data/ellipsoid_shapes.csv"
with open(file_name) as f:
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

# Run-times and Acceptance Rates
rng_row = 3
sph = [7, 12]
obt = [16, 23]
prt = [27, 34]
tri = [38, 45]
t_sph = np.zeros([2, 6])
t_obt = np.zeros([2, 8])
t_prt = np.zeros([2, 8])
t_tri = np.zeros([2, 8])
r_sph = np.zeros([2, 6])
r_obt = np.zeros([2, 8])
r_prt = np.zeros([2, 8])
r_tri = np.zeros([2, 8])
file_name = ["../data/timing_laggedfib.csv", "../data/timing_yarn.csv"]
for i in range(2):
    with open(file_name[i]) as f:
        r_no = 1
        j = 0
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            # RNG
            if row == rng_row:
                t_lf = row[1]
            # Sphere
            elif r_no >= sph[0] and r_no <= sph[1]:
                t_sph[i, j] = row[1]
                r_sph[i, j] = row[2]
                if r_no == sph[1]:
                    j = 0
                else:
                    j += 1
            # Oblate
            elif r_no >= obt[0] and r_no <= obt[1]:
                t_obt[i, j] = row[1]
                r_obt[i, j] = row[2]
                if r_no == obt[1]:
                    j = 0
                else:
                    j += 1
            # Prolate
            elif r_no >= prt[0] and r_no <= prt[1]:
                t_prt[i, j] = row[1]
                r_prt[i, j] = row[2]
                if r_no == prt[1]:
                    j = 0
                else:
                    j += 1
            # Triaxial
            elif r_no >= tri[0] and r_no <= tri[1]:
                t_tri[i, j] = row[1]
                r_tri[i, j] = row[2]
                if r_no == tri[1]:
                    j = 0
                else:
                    j += 1
            r_no += 1

# Relative standard deviation
sph = 4
obt = 8
prt = 10
tri = 12
rsd_sph = np.zeros([2, 6])
rsd_obt = np.zeros([2, 8])
rsd_prt = np.zeros([2, 8])
rsd_tri = np.zeros([2, 8])
file_name = ["../data/rsd_laggedfib.csv", "../data/rsd_yarn.csv"]
for i in range(2):
    with open(file_name[i]) as f:
        r_no = 1
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            if r_no == sph:
                for j in range(6):
                    rsd_sph[i, j] = row[j]
            elif r_no == obt:
                for j in range(8):
                    rsd_obt[i, j] = row[j]
            elif r_no == prt:
                for j in range(8):
                    rsd_prt[i, j] = row[j]
            elif r_no == tri:
                for j in range(8):
                    rsd_tri[i, j] = row[j] 
            r_no += 1

# chi-squared statistics
chi2_sph = np.zeros([2, 6])
chi2_obt = np.zeros([2, 8])
chi2_prt = np.zeros([2, 8])
chi2_tri = np.zeros([2, 8])
file_name = ["../data/chi2_tests_laggedfib.csv", "../data/chi2_tests_yarn.csv"]
for i in range(2):
    with open(file_name[i]) as f:
        r_no = 1
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            if r_no == sph:
                for j in range(6):
                    chi2_sph[i, j] = row[j]
            elif r_no == obt:
                for j in range(8):
                    chi2_obt[i, j] = row[j]
            elif r_no == prt:
                for j in range(8):
                    chi2_prt[i, j] = row[j]
            elif r_no == tri:
                for j in range(8):
                    chi2_tri[i, j] = row[j] 
            elif r_no == 14:
                chi2_crit = float(row[0])
            r_no += 1

### Prepare rows of the tables
sph_algo = ['Cubic Rej', 
            'Marsaglia',
            'Trig     ',
            'Gaussian ',
            'Cook     ',
            'Area Rej ']
ell_algo = ['Naive Scale',
            'Grad Rej   ',
            'Grad (Trig)',
            'Grad (Pol) ',
            'Area Rej   ',
            'Area (Pol) ',
            'Area (Merc)',
            'Ray Method ']
sph_lf = np.zeros([6, 6])
sph_yn = np.zeros([6, 6])
for i in range(6):
    # Lagged Fibonacci
    sph_lf[i][0] = t_sph[0][i]          # Run-time
    sph_lf[i][1] = t_sph[0][0] / t_sph[0][i] # Speed relative to cubic rejection
    sph_lf[i][2] = r_sph[0][i]          # Acceptance rate
    sph_lf[i][3] = rsd_sph[0][i] * 100  # Relative standard deviation
    sph_lf[i][4] = chi2_sph[0][i]       # chi^2 statistic (as a percenatge)
    # YARN
    sph_yn[i][0] = t_sph[1][i]          # Run-time
    sph_yn[i][1] = t_sph[1][0] / t_sph[1][i] # Speed relative to cubic rejection
    sph_yn[i][2] = r_sph[1][i]          # Acceptance rate
    sph_yn[i][3] = rsd_sph[1][i] * 100  # Relative standard deviation
    sph_yn[i][4] = chi2_sph[1][i]       # chi^2 statistic (as a percenatge)

obt_lf = np.zeros([8, 6])
obt_yn = np.zeros([8, 6])
for i in range(8):
    # Lagged Fibonacci
    obt_lf[i][0] = t_obt[0][i]          # Run-time
    obt_lf[i][1] = t_obt[0][0] / t_obt[0][i] # Speed relative to cubic rejection
    obt_lf[i][2] = r_obt[0][i]          # Acceptance rate
    obt_lf[i][3] = rsd_obt[0][i] * 100  # Relative standard deviation
    obt_lf[i][4] = chi2_obt[0][i]       # chi^2 statistic (as a percenatge)
    # YARN
    obt_yn[i][0] = t_obt[1][i]          # Run-time
    obt_yn[i][1] = t_obt[1][0] / t_obt[1][i] # Speed relative to cubic rejection
    obt_yn[i][2] = r_obt[1][i]          # Acceptance rate
    obt_yn[i][3] = rsd_obt[1][i] * 100  # Relative standard deviation
    obt_yn[i][4] = chi2_obt[1][i]       # chi^2 statistic (as a percenatge)    

tri_lf = np.zeros([8, 6])
tri_yn = np.zeros([8, 6])
for i in range(8):
    # Lagged Fibonacci
    tri_lf[i][0] = t_tri[0][i]          # Run-time
    tri_lf[i][1] = t_tri[0][0] / t_tri[0][i] # Speed relative to cubic rejection
    tri_lf[i][2] = r_tri[0][i]          # Acceptance rate
    tri_lf[i][3] = rsd_tri[0][i] * 100  # Relative standard deviation
    tri_lf[i][4] = chi2_tri[0][i]       # chi^2 statistic (as a percenatge)
    # YARN
    tri_yn[i][0] = t_tri[1][i]          # Run-time
    tri_yn[i][1] = t_tri[1][0] / t_tri[1][i] # Speed relative to cubic rejection
    tri_yn[i][2] = r_tri[1][i]          # Acceptance rate
    tri_yn[i][3] = rsd_tri[1][i] * 100  # Relative standard deviation
    tri_yn[i][4] = chi2_tri[1][i]       # chi^2 statistic (as a percenatge) 
    
    
### Write to .txt files in the form of LaTeX tables

# Spheres
with open("Tables/sphere_table.txt", "w") as file:
    file.write(r'\begin{table}[!ht]' + '\n')
    file.write('\t' + r'\centering' + '\n')
    file.write('\t' + r'\begin{tabular}{|l|l l l l l l|}' + '\n')
    file.write('\t\t' + r'\hline' + '\n')
    file.write('\t\t' + r'\multicolumn{7}{|c|}{\textbf{Sphere}: $(a,b,c) = (1,1,1)$} \\ \hline' + '\n')
    # Lagged Fibonacci
    file.write('\t\t' + r'\textbf{RNG} & \multicolumn{6}{|c|}{\textbf{Lagged Fibonacci} (lagfib4xor)} \\' + '\n')
    file.write('\t\t' + r'\textbf{Algorithm} & $t_\mathrm{run}$ (s) & speed & $r$ & RSD (\%) & $\chi^2$ & $\chi^2 < \chi_\mathrm{crit}^2$ \\ \hline' + '\n')
    for i in range(6):
        file.write('\t\t' + sph_algo[i] + ' & ')
        file.write(str(sph_lf[i][0]) + ' & ')
        file.write(str(round(sph_lf[i][1], 3)) + ' & ')
        file.write(str(round(sph_lf[i][2], 5)) + ' & ')
        file.write(str(round(sph_lf[i][3], 3)) + ' & ')
        file.write(str(round(sph_lf[i][4], 1)) + ' & ')
        if sph_lf[i][4] < chi2_crit:  # Has the chi^2 test passed?
            passed = "Yes"
        else:
            passed = "No"
        if i == 5:
            file.write(passed + r' \\ \hline' + '\n')
        else:
            file.write(passed + r' \\' + '\n')
    # YARN
    file.write('\t\t' + r'\textbf{RNG} & \multicolumn{6}{|c|}{\textbf{YARN5s} (yarn5s)} \\' + '\n')
    file.write('\t\t' + r'\textbf{Algorithm} & $t_\mathrm{run}$ (s) & speed & $r$ & RSD (\%) & $\chi^2$ & $\chi^2 < \chi_\mathrm{crit}^2$  \\ \hline' + '\n')
    for i in range(6):
        file.write('\t\t' + sph_algo[i] + ' & ')
        file.write(str(sph_yn[i][0]) + ' & ')
        file.write(str(round(sph_yn[i][1], 3)) + ' & ')
        file.write(str(round(sph_yn[i][2], 5)) + ' & ')
        file.write(str(round(sph_yn[i][3], 3)) + ' & ')
        file.write(str(round(sph_yn[i][4], 1)) + ' & ')
        if sph_yn[i][4] < chi2_crit:  # Has the chi^2 test passed?
            passed = "Yes"
        else:
            passed = "No"
        if i == 5:
            file.write(passed + r' \\ \hline' + '\n')
        else:
            file.write(passed + r' \\' + '\n')
    file.write('\t' + r'\end{tabular}' + '\n')
    file.write('\t' + r'\caption{ENTER CAPTION HERE! $\chi_\mathrm{crit}^2 = 65033.6$}' + '\n')
    file.write('\t' + r'\label{LABEL}' + '\n')
    file.write(r'\end{table}' + '\n')
    
# Oblate
with open("Tables/oblate_table.txt", "w") as file:
    file.write(r'\begin{table}[!ht]' + '\n')
    file.write('\t' + r'\centering' + '\n')
    file.write('\t' + r'\begin{tabular}{|l|l l l l l l|}' + '\n')
    file.write('\t\t' + r'\hline' + '\n')
    file.write('\t\t' + r'\multicolumn{7}{|c|}{\textbf{Oblate Spheroid}: $(a,b,c) = (3,3,1.5)$} \\ \hline' + '\n')
    # Lagged Fibonacci
    file.write('\t\t' + r'\textbf{RNG} & \multicolumn{6}{|c|}{\textbf{Lagged Fibonacci} (lagfib4xor)} \\' + '\n')
    file.write('\t\t' + r'\textbf{Algorithm} & $t_\mathrm{run}$ (s) & speed & $r$ & RSD (\%) & $\chi^2$ & $\chi^2 < \chi_\mathrm{crit}^2$ \\ \hline' + '\n')
    for i in range(8):
        file.write('\t\t' + ell_algo[i] + ' & ')
        file.write(str(obt_lf[i][0]) + ' & ')
        file.write(str(round(obt_lf[i][1], 3)) + ' & ')
        file.write(str(round(obt_lf[i][2], 5)) + ' & ')
        file.write(str(round(obt_lf[i][3], 3)) + ' & ')
        file.write(str(round(obt_lf[i][4], 1)) + ' & ')
        if obt_lf[i][4] < chi2_crit:  # Has the chi^2 test passed?
            passed = "Yes"
        else:
            passed = "No"
        if i == 7:
            file.write(passed + r' \\ \hline' + '\n')
        else:
            file.write(passed + r' \\' + '\n')
    # YARN
    file.write('\t\t' + r'\textbf{RNG} & \multicolumn{6}{|c|}{\textbf{YARN5s} (yarn5s)} \\' + '\n')
    file.write('\t\t' + r'\textbf{Algorithm} & $t_\mathrm{run}$ (s) & speed & $r$ & RSD (\%) & $\chi^2$ & $\chi^2 < \chi_\mathrm{crit}^2$  \\ \hline' + '\n')
    for i in range(8):
        file.write('\t\t' + ell_algo[i] + ' & ')
        file.write(str(obt_yn[i][0]) + ' & ')
        file.write(str(round(obt_yn[i][1], 3)) + ' & ')
        file.write(str(round(obt_yn[i][2], 5)) + ' & ')
        file.write(str(round(obt_yn[i][3], 3)) + ' & ')
        file.write(str(round(obt_yn[i][4], 1)) + ' & ')
        if obt_yn[i][4] < chi2_crit:  # Has the chi^2 test passed?
            passed = "Yes"
        else:
            passed = "No"
        if i == 7:
            file.write(passed + r' \\ \hline' + '\n')
        else:
            file.write(passed + r' \\' + '\n')
    file.write('\t' + r'\end{tabular}' + '\n')
    file.write('\t' + r'\caption{ENTER CAPTION HERE! $\chi_\mathrm{crit}^2 = 65033.6$}' + '\n')
    file.write('\t' + r'\label{LABEL}' + '\n')
    file.write(r'\end{table}' + '\n')

# Triaxial
with open("Tables/triaxial_table.txt", "w") as file:
    file.write(r'\begin{table}[!ht]' + '\n')
    file.write('\t' + r'\centering' + '\n')
    file.write('\t' + r'\begin{tabular}{|l|l l l l l l|}' + '\n')
    file.write('\t\t' + r'\hline' + '\n')
    file.write('\t\t' + r'\multicolumn{7}{|c|}{\textbf{Triaxial Ellipsoid}: $(a,b,c) = (3,2,1)$} \\ \hline' + '\n')
    # Lagged Fibonacci
    file.write('\t\t' + r'\textbf{RNG} & \multicolumn{6}{|c|}{\textbf{Lagged Fibonacci} (lagfib4xor)} \\' + '\n')
    file.write('\t\t' + r'\textbf{Algorithm} & $t_\mathrm{run}$ (s) & speed & $r$ & RSD (\%) & $\chi^2$ & $\chi^2 < \chi_\mathrm{crit}^2$ \\ \hline' + '\n')
    for i in range(8):
        file.write('\t\t' + ell_algo[i] + ' & ')
        file.write(str(tri_lf[i][0]) + ' & ')
        file.write(str(round(tri_lf[i][1], 3)) + ' & ')
        file.write(str(round(tri_lf[i][2], 5)) + ' & ')
        file.write(str(round(tri_lf[i][3], 3)) + ' & ')
        file.write(str(round(tri_lf[i][4], 1)) + ' & ')
        if tri_lf[i][4] < chi2_crit:  # Has the chi^2 test passed?
            passed = "Yes"
        else:
            passed = "No"
        if i == 7:
            file.write(passed + r' \\ \hline' + '\n')
        else:
            file.write(passed + r' \\' + '\n')
    # YARN
    file.write('\t\t' + r'\textbf{RNG} & \multicolumn{6}{|c|}{\textbf{YARN5s} (yarn5s)} \\' + '\n')
    file.write('\t\t' + r'\textbf{Algorithm} & $t_\mathrm{run}$ (s) & speed & $r$ & RSD (\%) & $\chi^2$ & $\chi^2 < \chi_\mathrm{crit}^2$  \\ \hline' + '\n')
    for i in range(8):
        file.write('\t\t' + ell_algo[i] + ' & ')
        file.write(str(tri_yn[i][0]) + ' & ')
        file.write(str(round(tri_yn[i][1], 3)) + ' & ')
        file.write(str(round(tri_yn[i][2], 5)) + ' & ')
        file.write(str(round(tri_yn[i][3], 3)) + ' & ')
        file.write(str(round(tri_yn[i][4], 1)) + ' & ')
        if tri_yn[i][4] < chi2_crit:  # Has the chi^2 test passed?
            passed = "Yes"
        else:
            passed = "No"
        if i == 7:
            file.write(passed + r' \\ \hline' + '\n')
        else:
            file.write(passed + r' \\' + '\n')
    file.write('\t' + r'\end{tabular}' + '\n')
    file.write('\t' + r'\caption{ENTER CAPTION HERE! $\chi_\mathrm{crit}^2 = 65033.6$}' + '\n')
    file.write('\t' + r'\label{LABEL}' + '\n')
    file.write(r'\end{table}' + '\n')










