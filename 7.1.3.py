# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 22:03:17 2021

@author: savan
"""
#Basic model, Different yield vals

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf
import matplotlib.pyplot as plt

def main():
## experimental parameters
    D = 3E-3
    #rho_n = 3E-1     
    rc = 3.5E-2
    Dc = 1E-5       
    w = 1           

    environment_size = (59, 59)
    

## add nutrient to the plate
    
## add sender strain to the plate
    colours = ['b', 'r', 'g', 'y', 'k']
    #fig = 
    fig, axs = plt.subplots(1, 2)
    plt.suptitle('Change in concentration of species over time, with varying yield of sender:nutrient')
    plt.style.use('ggplot')
    yields = [1/10,2/10,3/10]
    
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        n = D * hf.ficks(species['N'], w) - species['S'] * rho_n * hf.leaky_hill(s=species['N'], K=0.15, lam=1, max=rc, min=0)
        return n
    
        
    def S_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        sp = Dc * hf.ficks(species['S'], w)*1.68 + species['S'] * hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
        return sp
       
    for col, yi in enumerate(yields):
        plate = Plate(environment_size)
        rho_n = yi
        U_N = np.ones(environment_size)
        N = Species("N", U_N)
        N.set_behaviour(N_behaviour)
        plate.add_species(N)

        U_S = np.zeros(environment_size)
        for i in np.linspace(29, 29, 1):
            U_S[int(i), int(i)] = 1
        S = Species("S", U_S)
        
        S.set_behaviour(S_behaviour)
        plate.add_species(S)

    ## run the experiment
        params = (D, rho_n, Dc, rc, w)
        sim = plate.run(t_final = 150*60,
                    dt = 1.,
                    params = params)
    
        #plate.plot_simulation(sim, 10)
        
        Sp = plate.get_all_species()
        #plate.plot_conc_distribution(sim, Sp, 10)
        colour = colours[col]
        plate.compare_species(sim, Sp, 10, fig, axs, colour)
    
    fig.legend(labels=yields, title='Yield: species/nutrient', loc='center right')
 
    fig.show()
main()