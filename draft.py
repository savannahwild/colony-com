# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 09:16:26 2021

@author: savan
"""

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf


def main():
    ## experimental parameters
    D = 3E-3        # nutrient diffusion coeff (#mm2/min)
    rho_n = 0.3     # consumption rate of nutrients by X
    rc = 6E-3       # growth rate of X on N
    Dc = 1E-5       # cell diffusion coefficient
    w = 1
    Da = 0.03       #diffusion ahl
    rho_A = 0.1       # production rate of AHL
    Dt = 0.05   #diffusion from target

    environment_size = (50, 50)
    plate = Plate(environment_size)

    ## add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        n = D * hf.ficks(species['N'], w) - rho_n * species['N'] * (species['R'] + species['S'])
        return n
    N.set_behaviour(N_behaviour)
    plate.add_species(N)

    ## add uninformed strain to the plate
    U_R = np.zeros(environment_size)
    for i in np.linspace(5, 45, 9):
        for j in np.linspace(5, 45, 9):
            if (i == 25) & (j == 25):
                continue
            U_R[int(i), int(j)] = 0.001
    R = Species("R", U_R)
    def R_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        r = Dc * hf.ficks(species['R'], w) + rc * species['N'] * species['R']
        return r
    R.set_behaviour(R_behaviour)
    plate.add_species(R)

    ## add informed strain to the plate
    U_S = np.zeros(environment_size)
    U_S[25, 25] = 0.001
    S = Species("S", U_S)
    def S_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        r = Dc * hf.ficks(species['S'], w) + rc * species['N'] * species['S']
        return r
    S.set_behaviour(S_behaviour)
    plate.add_species(S)

    ##add target to plate
    U_T = np.zeros(environment_size)
    grad_values = np.logspace(-4, 5, environment_size[0])
    for idx, value in enumerate(grad_values):
        U_T[idx,:] = value

    T = Species("T", U_T)
    def T_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        t = Dt * hf.ficks(species['T'], w)
        #t = 0
        return t
    T.set_behaviour(T_behaviour)
    plate.add_species(T)
    
    ## add AHL to plate
    U_A = np.zeros(environment_size)
    A = Species("A", U_A)
    def A_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        a = Da * hf.ficks(species['A'], w) + rho_A * species['S']
        return a
    A.set_behaviour(A_behaviour)
    plate.add_species(A)

    ## add GFP FOR RECEIVER TO plate
    U_G = np.zeros(environment_size)
    G = Species("G", U_G)
    def G_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        g = Dc * hf.ficks(species['G'], w) + hf.leaky_hill(s = species['A'], K = 1E-3, lam = 2, min = 1E-3, max=1) * species['R']
        return g
    G.set_behaviour(G_behaviour)
    plate.add_species(G)

    ## run the experiment
    params = (D, rho_n, Dc, rc, w, rho_A, Da, Dt)
    sim = plate.run(t_final = 28*60,
                    dt = .1,
                    params = params)

    ## plotting
    tp = np.arange(0, 28, 4)
    plate.plot_simulation(sim, tp)

main()