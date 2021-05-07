# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 23:32:40 2021

@author: savan
"""
#introduction of target everywhere, different concs

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf

def main():
    ## experimental parameters
    D = 0.03
    rho_n = 0.5     
    rc = 3.33e-2 #6e-4
    Dc = 1e-5       
    w = 0.75 
    Da = 2.94e-6
    rho_A = 0.01    
    Dt = 0
    
    environment_size = (59, 59)
    concs = np.linspace(0,10,5)

    def S_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            s = Dc*hf.ficks(species['S'], w)*hf.leaky_hill(s=species['T'], K=0.6, lam=2, max=2.75, min=1) + species['S'] * hf.leaky_hill(s=species['N'], K=80, lam = 1, max=rc, min=0)
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
            n = D * hf.ficks(species['N'], w) - species['S'] * rho_n * hf.leaky_hill(s=species['N'], K=80, lam=1, max=rc, min=0)
            return n
        
    for col, conc in enumerate(concs):
        plate = Plate(environment_size)
        
        ## add nutrient to the plate
        U_N = np.zeros(environment_size)
    
        for i in np.linspace(0, 58, 59):
            for j in np.linspace(0,58,59):
                U_N[int(i),int(j)]=100
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
        sim = plate.run(t_final = 200*60,
                        dt = 1.,
                        params = params)

        #plate.plot_simulation(sim, 3)
        S = plate.get_all_species()
        plate.compare_species(sim, S, 10, col)

main()