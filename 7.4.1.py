
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:42:37 2021

@author: savan
"""
#ahl production, receiver on plate but no detection

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf
import matplotlib.pyplot as plt

def main():
    
## experimental parameters
    D = 3E-3        # nutrient diffusion coeff (#mm2/min) maybe?
    rho_n = 0.3     #consumption rate of nutrients by X calc?
    rc = 1E-2     # growth rate of X divisions is max 0.0352
    Dc = 1E-5       # cell diffusion coefficient? calc that 0.03
    w = 1           #w = time?
    Da = 0.0294     #mm2/min
    rho_A = 0.01    #production rate of AHL
    Dt = 6e-2    #diffusion rate of target? change

    environment_size = (59, 59)
    plate = Plate(environment_size)
    #fig, axs = plt.subplots(1,5)
    fig2, axs2 = plt.subplots(1,4)
    
    ##add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        n = D*hf.ficks(species['N'], w) - (species['R']+species['S'])*rho_n*hf.leaky_hill(s=species['N'], K=0.15, lam=1, max=rc, min=0) #rho_n * species['N'] * (species['R'])
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
        s = Dc*hf.ficks(species['S'], w)*1.68 + species['S']*hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
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
        r = Dc*hf.ficks(species['R'], w)*1.68 + species['R']*hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
        return r
    R.set_behaviour(R_behaviour)
    plate.add_species(R)
    
    
    ##add AHL to plate
    U_A = np.zeros(environment_size)
    A = Species("A", U_A)
    def A_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        #a = Da * hf.ficks(species['A'], w) + hf.leaky_hill(s=species['T'], K=1, lam=2, max=rho_A, min=0)
        a = Da*hf.ficks(species['A'], w) + rho_A*species['S']
        return a
    A.set_behaviour(A_behaviour)
    plate.add_species(A)

    params = (D, rho_n, Dc, rc, w, rho_A, Da, Dt)
    sim = plate.run(t_final = 150*60,
                    dt = 1.,
                    params = params)

    #plate.plot_simulation(sim, 10)
    colour = 'b'
    S = plate.get_all_species()
    #plate.plot_conc_distribution(sim, S, 10,fig,axs,colour)
    plate.compare_species(sim, S, 10,fig2,axs2,colour)

main()