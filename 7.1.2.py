# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:41:10 2021

@author: savan
"""
#Basic model, Different inital conc S

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf
import matplotlib.pyplot as plt

def main():
## experimental parameters
    D = 3E-3
    rho_n = 3E-1     
    rc = 3.5E-2
    Dc = 1E-5       
    w = 1           

    environment_size = (59, 59)
    

## add nutrient to the plate
    
## add sender strain to the plate
    
    
    colours = ['b', 'r', 'g', 'y', 'k']
    fig, axs = plt.subplots(1, 2)
    plt.suptitle('Change in concentration of species over time, with varying initial concentration of sender')
    plt.style.use('ggplot')
    concs = [0, 3/10, 1]
    
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        n = D * hf.ficks(species['N'], w) - species['S'] * rho_n * hf.leaky_hill(s=species['N'], K=0.15, lam=1, max=rc, min=0)
        return n
    
        
    def S_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        sp = Dc*hf.ficks(species['S'], w)*1.68 + species['S']*hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
        return sp
    fig2, axs2 = plt.subplots(1, 3)
    for col, conc in enumerate(concs):
        
        #plt.suptitle('Change in average concentration of each species over time in each plate quadrant')
        plate = Plate(environment_size)
        
        
        U_N = np.ones(environment_size)
        N = Species("N", U_N)
        N.set_behaviour(N_behaviour)
        plate.add_species(N)

        U_S = np.zeros(environment_size)
        for i in np.linspace(0, 58, 59):
            for j in np.linspace(0, 58, 59):
                U_S[int(i), int(j)] = conc
        S = Species("S", U_S)
        
        S.set_behaviour(S_behaviour)
        plate.add_species(S)

    ## run the experiment
        params = (D, rho_n, Dc, rc, w)
        sim = plate.run(t_final = 2*60,
                    dt = 1.,
                    params = params)
    
        #plate.plot_simulation(sim, 3)
      
        S = plate.get_all_species()
        colour = colours[col]
        #plate.plot_conc_distribution(sim, S, 10, fig2, axs2[col], colour)
        plate.compare_species(sim, S, 10, fig, axs, colour)
    fig.legend(labels=concs, title='Initial concentration of S', loc='center right')
    fig.show()
main()