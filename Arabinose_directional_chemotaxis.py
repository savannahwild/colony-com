# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 17:38:20 2021

@author: savan
"""

#ahl production from sender + induced receiver motility via chemotaxis 
#gradient of inducer arabinose for sender directional chemotaxis

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
    plate = Plate(environment_size)                
    concs = np.linspace(0,1,2)
    for col, conc in enumerate(concs):
        ##add nutrient to the plate
        U_N = np.zeros(environment_size)
        
        for i in np.linspace(0, 58, 59):
            for j in np.linspace(0,58,59):
                U_N[int(i),int(j)]=100
        N = Species("N", U_N)
        def N_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            n = D*hf.ficks(species['N'], w) - (species['R']+species['S'])*rho_n*hf.leaky_hill(s=species['N'], K=80, lam=1, max=rc, min=0)
            return n
        N.set_behaviour(N_behaviour)
        plate.add_species(N)
    
    
        ##add sender strain to the plate
        U_S = np.zeros(environment_size)
        for i in np.linspace(29, 29, 1):
            U_S[int(i), int(i)] = 0.001
        S = Species("S", U_S)
        def S_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            s = Dc * hf.ficks(species['S'], w)*hf.leaky_hill(s=species['T'], K=0.6, lam=2, max=2.75, min=1) + species['S']*hf.leaky_hill(s=species['N'], K=80, lam = 1, max=rc, min=0)
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
            r = Dc*hf.leaky_hill(s=species['A'], K=40e-9, lam=2, max=2.75, min=1)*hf.ficks(species['R'], w) + species['R']*hf.leaky_hill(s=species['N'], K=80, lam = 1, max=rc, min=0)
            return r
        R.set_behaviour(R_behaviour)
        plate.add_species(R)
        
        ##add AHL to plate
        U_A = np.zeros(environment_size)
        A = Species("A", U_A)
        def A_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            a = conc*(Da * hf.ficks(species['A'], w) + rho_A*species['S'])
            return a
        A.set_behaviour(A_behaviour)
        plate.add_species(A)
    
        
        ##add target to plate
        U_T = np.zeros(environment_size)
        for j in np.linspace(0,58, environment_size[0]):
                for i in np.linspace(30, 58, 29):
                    U_T[int(i), int(j)] = (int(i)-29)/2.9
        T = Species("T", U_T)
        def T_behaviour(species, params):
            ## unpack params
            D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
            t = Dt*hf.ficks(species['T'], w)
            return t
        T.set_behaviour(T_behaviour)
        plate.add_species(T)
        ##target need induce more motility than ahl? different k 2/3rd than ahl
       
        params = (D, rho_n, Dc, rc, w, rho_A, Da, Dt)
        sim = plate.run(t_final = 200*60,
                        dt = 1.,
                        params = params)
        
    
        #plate.plot_simulation(sim, 3)
        idx=col
        S = plate.get_all_species()
        plate.plot_conc_target(sim, S, 10,idx)
        plate.compare_species(sim, S, 10,idx)
main()