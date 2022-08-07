#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 15:33:33 2021

@author: edward
"""

import numpy as np
import pandas as pd
from class_crypt_mut import Crypt_with_mutatons
import copy

## Simulate crypts with fission
## Patch is a list of crypt objects. crypt class at every iteration updates mutations and time step 

# #########################################
# # sim params
# loci  = 26557613
# alpha = 4e-8

# # WT no fission
# N = 7; lamb = 0.027; rho = 0
# file_sims = "sims/no_fission_WT.csv"

# # WT fission
# # N = 7; lamb = 0.027; rho = 0.027/30 
# # file_sims = "sims/fission_WT.csv"

# # KRAS fission
# # rho = 0.1
# # WT rho = 0.032
# # lamb_WT , N_WT
# # lamb_KRAS, N_KRAS


# # t_vals = [25, 100, 200, 300]
# # cells_tp = 200

# t_vals = [5, 10, 25, 50, 100, 150, 200, 250, 300]
# cells_tp = 1000

# clonal_threshold = 0.99
# one_cell_organoid = False 
# #########################################

def run_sim(loci, alpha, N, lamb, rho, file_sims, t_vals, cells_tp, clonal_threshold, one_cell_organoid):
            
    tp_measure = np.concatenate([ti*np.ones(cells_tp) for ti in t_vals])
    vec_out    = np.zeros((len(tp_measure), 4))
    
    for indx, t2 in enumerate(tp_measure):
        crypt_list = [Crypt_with_mutatons(lamb, rho, N, loci, alpha, t2, one_cell_organoid)]
        all_crypts_finished = 0
        # Store some sim values 
        while all_crypts_finished != 1:        
            fission_list = []
            for indx_crypt, crypt_i in enumerate(crypt_list):
                # Run mutation sim for crypt i until t2 or fission event. if t2 reached, update will do nothing
                crypt_i.update()
                # Check for fission 
                if crypt_i.is_fission():
                    fission_list.append(indx_crypt)
                    crypt_i.reset_fission()
            
            # If any fissioned make copy and add to list
            for indx_crypt in fission_list:
                new_crypt = copy.deepcopy(crypt_list[indx_crypt])
                crypt_list.append(new_crypt)
                
            # Proportion of crypts finished
            all_crypts_finished = np.mean([crypt_i.is_finished() for crypt_i in crypt_list])
        # summarise and store ------------
        vec_out[indx,1] = t2
        # Record patch size 
        vec_out[indx,2] = len(crypt_list)
        
        # Mutations: mean, sum, different number
        all_mut = [crypt_i.get_mutations() for crypt_i in crypt_list]
        all_mut = np.concatenate(all_mut)
        
        mut_ids, mut_count = np.unique(all_mut, return_counts=True)
        if not one_cell_organoid:
            num_mut_i = np.sum(mut_count > clonal_threshold*N)                            
        else:
            num_mut_i = len(mut_count)                                    
        vec_out[indx,0] = num_mut_i
        vec_out[indx,3] = len(mut_count) ## Any mutation
        
    
    df = pd.DataFrame(data=vec_out, columns=["num_mut_clonal", "t_measure", "patch_size", "mut_any_VAF"])
    df["N"] = N
    df["lamb"] = lamb
    df["loci"] = loci
    df["alpha"] = alpha
    df["det_thresh"] = clonal_threshold
    df["one_cell_organoid"] = one_cell_organoid
    
    df.to_csv(file_sims,   index = False)




