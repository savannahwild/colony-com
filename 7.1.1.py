# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 12:32:07 2021

@author: savan
"""
#Basic model, N and S

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
    plate = Plate(environment_size)
    
    colours = ['b', 'r', 'g', 'y', 'k']
    fig, axs = plt.subplots(1, 2)
    plt.suptitle('Change in concentration of species over time')
    fig2, axs2 = plt.subplots(1, 2)
    plt.suptitle('Change in concentration of species over time in each plate quadrant')
## add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        n = D * hf.ficks(species['N'], w) - species['S'] * rho_n * hf.leaky_hill(s=species['N'], K=0.15, lam=1, max=rc, min=0)
        return n
    N.set_behaviour(N_behaviour)
    plate.add_species(N)

## add sender strain to the plate
    U_S = np.zeros(environment_size)
    for i in np.linspace(29, 29, 1):
        U_S[int(i), int(i)] = 1
    S = Species("S", U_S)
    def S_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        s = Dc * hf.ficks(species['S'], w) + species['S'] * hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
        return s
    S.set_behaviour(S_behaviour)
    plate.add_species(S)

## run the experiment
    params = (D, rho_n, Dc, rc, w)
    sim = plate.run(t_final = 150*60,
                    dt = 1.,
                    params = params)
    
    
    print("run simulation")
    plate.plot_simulation(sim, 10)

    print("plot conc dist")
    S = plate.get_all_species()
    colour = colours[0]
    plate.compare_species(sim, S, 10, fig, axs, colour)
    plate.plot_conc_distribution(sim, S, 10, fig2, axs2, colour)

main()