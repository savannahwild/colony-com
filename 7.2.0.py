# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:57:35 2021

@author: savan
"""

#target everywhere to show motility not growth affected

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf
import matplotlib.pyplot as plt

def main():
    ## experimental parameters
    D = 3E-3        # nutrient diffusion coeff (#mm2/min)
    rho_n = 0.3     # consumption rate of nutrients by X
    rc = 3.5E-2       # growth rate of X on N
    Dc = 1E-5       # cell diffusion coefficient
    w = 1   
    Da = 0.03   #0.05 1/sec
    rho_A = 0.1       # production rate of AHL synthesis rate 3,300nM/hr: degradation rate 0.108/hour
    Dt = 0     #random difusion rate of target molecule #K 25 - 50 nM AHL

    environment_size = (59, 59)
    colours = ['b', 'r', 'g', 'y', 'k']
    concs = [0,1]
    fig, axs = plt.subplots(1,3)
    plt.suptitle('Change in concentration of species over time with inducer present')

    def S_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            s = Dc*hf.ficks(species['S'], w)*hf.leaky_hill(s=species['T'], K=0.6, lam = 2, max=3.96, min=1.68) + species['S']*hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
            return s
    def T_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            t = Dt * hf.ficks(species['T'], w)
            #t = 0
            return t
    def N_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            n = D * hf.ficks(species['N'], w) - species['S'] * rho_n * hf.leaky_hill(s=species['N'], K=0.15, lam=1, max=rc, min=0)
            return n
        
    for col, conc in enumerate(concs):
        plate = Plate(environment_size)
    
        ## add nutrient to the plate
        U_N = np.ones(environment_size)
        N = Species("N", U_N)
        
        N.set_behaviour(N_behaviour)
        plate.add_species(N)
        
        ## add sender strain to the plate
        U_S = np.zeros(environment_size)
        for i in np.linspace(29, 29, 1):
            for j in np.linspace(29, 29, 1):
                U_S[int(i), int(j)] = 0.001
        S = Species("S", U_S)
        
        S.set_behaviour(S_behaviour)
        plate.add_species(S)
    
       ## add target to plate
        U_T = np.zeros(environment_size)
        for i in np.linspace(0, 58, 59):
            for j in np.linspace(0, 58, 59):
                U_T[int(i), int(j)] = conc
        T = Species("T", U_T)
        T.set_behaviour(T_behaviour)
        plate.add_species(T)
        ## run the experiment
        params = (D, rho_n, Dc, rc, w, rho_A, Da, Dt)
        sim = plate.run(t_final = 150*60, #hours times 60
                        dt = 1.,
                        params = params)
        print("run simulation")

        #plate.plot_simulation(sim, 10)
 
        S = plate.get_all_species()
        colour = colours[col]
        
        #plate.plot_conc_distribution(sim, S, 10, fig, axs, colour)
        plate.compare_species(sim, S, 10, fig, axs, colour)
    fig.legend(labels=concs, title='Initial concentration of inducer')
    fig.show()
main()