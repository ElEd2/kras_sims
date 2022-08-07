# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:56:49 2021

@author: edward
"""
from sim_diff_mut_w_fission import run_sim
import numpy as np
#########################################
# sim params
loci  = 2.7e9 # From Clevers paper
alpha = 1.25037e-08

t_vals = np.arange(30) #[5, 10, 25, 50, 100, 150, 200, 250, 300]
cells_tp = 1000

clonal_threshold = 0.99
one_cell_organoid = False 
#########################################


# WT no fission
N = 7; lamb = 1.3; rho = 0
file_sims = "sims/human_no_fission_WT.csv"
run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid)

# WT fission
N = 7; lamb = 1.3; rho = 0.007 
file_sims = "sims/human_fission_WT.csv"
run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid)

# Use lamb gain  1.3 * 0.263/0.207
# KRAS no fission
N = 5; lamb = 1.3 * 0.263/0.207; rho = 0
file_sims = "sims/human_no_fission_kras_prop.csv"
run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid)

# KRAS fission
N = 5; lamb = 1.3 * 0.263/0.207; rho = 0.007*11
file_sims = "sims/human_fission_kras_prop.csv"
run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid)

# KRAS max lamb fission, mouse value ...
N = 5; lamb = 0.263*365; rho = 0.007*11
file_sims = "sims/human_fission_kras_max.csv"
run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid)


