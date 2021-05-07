# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:41:10 2021

@author: savan
"""
#Basic model, sender majority

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
    concs=[1e-1,1e-2,1e-3,1e-4,1e-5]
    
    for col, conc in enumerate(concs):
        
        #plt.suptitle('Change in average concentration of each species over time in each plate quadrant')
        plate = Plate(environment_size)
        ##add nutrient to the plate
        U_N = np.zeros(environment_size)
        for i in np.linspace(0, 58, 59):
            for j in np.linspace(0,58,59):
                U_N[int(i),int(j)]=100
        N = Species("N", U_N)
        def N_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            n = D*hf.ficks(species['N'], w) - (species['R']+species['S'])*rho_n*hf.leaky_hill(s=species['N'], K=80, lam=1, max=rc, min=0) #rho_n * species['N'] * (species['R'])
            return n
        N.set_behaviour(N_behaviour)
        plate.add_species(N)
    
        ##add sender strain to the plate
        U_S = np.zeros(environment_size)
        for i in np.linspace(29, 29, 1):
            U_S[int(i), int(i)] = conc*0.001
        S = Species("S", U_S)
        def S_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            s = Dc*hf.ficks(species['S'], w) + species['S']*hf.leaky_hill(s=species['N'], K=80, lam = 1, max=rc, min=0)
            return s
        S.set_behaviour(S_behaviour)
        plate.add_species(S)
        
        ##add uninformed strain to the plate
        U_R = np.zeros(environment_size)
        for i in np.linspace(29, 29, 1):
            for j in np.linspace(29,29, 1):
                U_R[int(i), int(j)] = 0.001
        R = Species("R", U_R)
        def R_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            r = Dc*hf.ficks(species['R'], w) + species['R']*hf.leaky_hill(s=species['N'], K=80, lam = 1, max=rc, min=0)
            return r
        R.set_behaviour(R_behaviour)
        plate.add_species(R)

        ## run the experiment
        params = (D, rho_n, Dc, rc, w, rho_A, Da, Dt)
      
        S = plate.get_all_species()
        sim = plate.run(t_final = 200*60,
                    dt = 1.,
                    params = params)

        #plate.plot_simulation(sim, 10)
        S = plate.get_all_species()
        idx=col
        plate.compare_species(sim, S, 10, idx)
        #plate.plot_conc_target(sim, S, 10, idx)
main()