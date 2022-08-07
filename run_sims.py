#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:56:49 2021

@author: edward
"""
from sim_diff_mut_w_fission import run_sim
#########################################
# sim params old
# loci  = 26557613 # Bases covered in an organoid from sequencing file
# alpha = 4e-8

# sim params upadetd
loci  = 30e6 # Bases covered in an organoid from sequencing file
alpha = 4e-8

t_vals = [5, 10, 25, 50, 100, 150, 200, 250, 300]
cells_tp = 1000

clonal_threshold = 0.99
one_cell_organoid = False 
#########################################


# WT no fission
N = 7; lamb = 0.207; rho = 0
file_sims = "sims/no_fission_WT.csv"
run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid)

# WT fission
N = 7; lamb = 0.207; rho = 0.032/30 
file_sims = "sims/fission_WT.csv"
run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid)

# KRAS no fission
N = 5; lamb = 0.263; rho = 0
file_sims = "sims/no_fission_kras.csv"
run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid)

# KRAS fission
N = 5; lamb = 0.263; rho = 0.103/30 
file_sims = "sims/fission_kras.csv"
run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid)



