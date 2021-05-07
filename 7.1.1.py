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

    plate.plot_simulation(sim, 10)
    
    colour = colours[0]
    S = plate.get_all_species()
    plate.compare_species(sim, S, 10, colour)
    plate.plot_conc_target(sim, S, 10)

main()