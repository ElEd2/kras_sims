#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 16:10:24 2021

@author: edward
"""
import numpy as np

class Crypt_with_mutatons:
    def __init__(self, lamb, rho, N, Nloci, alpha_loci, t2, one_cell_organoid=False):
        self.sim_time = 0
        self.fission  = False
        self.finished = False
        self.cell_mut_list = [[] for i in range(N)]
        self.lamb, self.rho, self.N, self.Nloci, self.alpha_loci, self.t2 = lamb, rho, N, Nloci, alpha_loci, t2
        self.one_cell_organoid = one_cell_organoid
        self.total_rate = N*lamb + rho
        
        # Make event matrix, first value fission, the rest replacemnet: from to
        range_N = np.arange(0, N)
        self.event_list = np.vstack([[-1, 0], 
                                     np.array([range_N, range_N +1]).T, 
                                     np.array([range_N, range_N -1]).T])
        self.event_list[  self.event_list[:,1] < 0 , 1]   = N - 1
        self.event_list[self.event_list[:,1]   > N-1 , 1] = 0
        
        # Prob vec
        self.prob_vec = np.concatenate([[rho], 0.5*lamb*np.ones(2*N)])
        self.prob_vec = self.prob_vec/sum(self.prob_vec)
        # Event vec 
        self.indx_vec = np.arange(0,2*N+1)
        
    # Update until end or until fission event
    def update(self):
        while self.sim_time <= self.t2 and not self.fission:
            r  =  np.random.uniform(0, 1)
            dt = -np.log(r)/self.total_rate
            self.sim_time += dt
            
            # Fission or replacement, plus which replacment
            event_num = np.random.choice(self.indx_vec, p = self.prob_vec) 
            
            # Check to see if fission
            if event_num ==0:
                # Fission, activate fission flag and exit sim
                self.fission = True
                break
                           
            # Choose cell
            from_N = self.event_list[event_num, 0]
            # Choose left or right replace
            to_N   = self.event_list[event_num, 1]
            # to_N = from_N
            # Choose if mutate
            dmut     = np.random.binomial(self.Nloci, self.alpha_loci)
            indx_mut = np.random.randint(1, self.Nloci, size = dmut)

            if self.sim_time < self.t2 and event_num!=0:
                # Replacement and mutation
                self.cell_mut_list[to_N] = self.cell_mut_list[from_N].copy() 
                if dmut>0: self.cell_mut_list[to_N].append(indx_mut)
            
        # Might have ended due to fission before t2
        if self.sim_time > self.t2:
            self.finished = True
                
    def is_fission(self):
        return self.fission
        
    def reset_fission(self):
        self.fission = False
        
    def is_finished(self):
        return self.finished
        
    def get_mutations(self):
        all_mut = [np.concatenate(i) if len(i) >0 else [] for i in self.cell_mut_list]
        if self.one_cell_organoid:
            # Choose one
            all_mut = all_mut[np.random.randint(self.N)]
        else:
            all_mut = np.concatenate(all_mut)
        return all_mut




# To test sims
#########################################
# sim params
# N     = 7
# loci  = 26557613
# alpha = 4e-8
# lamb  = 0.21
# rho   = 0
# t2    = 200
# est_num = loci*alpha*lamb*t2
# one_cell_organoid = False

# crypt_sim = Crypt_with_mutatons(lamb, rho, N, loci, alpha, t2, one_cell_organoid)    

# crypt_sim.update()
# all_mut = crypt_sim.get_mutations()

# mut_ids, mut_count = np.unique(all_mut, return_counts=True)
# print(est_num)
# print(mut_count.shape)



