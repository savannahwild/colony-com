4E4# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 19:13:54 2021

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
    rc = 3.5E-2       # growth rate of X on N
    Dc = 1E-5       # cell diffusion coefficient
    w = 1
    Da = 0.03
    rho_A = 0.1       # production rate of AHL
    Dt = 0.01     #random difusion rate of target molecule

    environment_size = (60, 60)
    plate = Plate(environment_size)

    ## add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        n = D * hf.ficks(species['N'], w) - species['S'] * rho_n * hf.leaky_hill(s=species['N'], K=0.15, lam=1, max=rc, min=0)
        return n
    N.set_behaviour(N_behaviour)
    plate.add_species(N)

    ## add sender strain to the plate
    U_S = np.zeros(environment_size)
    for i in np.linspace(30, 30, 1):
        U_S[int(i), int(i)] = 1
    S = Species("S", U_S)
    def S_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        s = Dc * hf.ficks(species['S'], w)* hf.leaky_hill(s=species['T'], K=0.15, lam = 2, max=2, min=0) + species['S'] * hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
        return s
    S.set_behaviour(S_behaviour)
    plate.add_species(S)

    ## add target to plate
    U_T = np.ones(environment_size)
    T = Species("T", U_T)
    def T_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        #t = Dt * hf.ficks(species['T'], w)
        t = 0
        return t
    T.set_behaviour(T_behaviour)
    plate.add_species(T)

    ## run the experiment
    params = (D, rho_n, Dc, rc, w, rho_A, Da, Dt)
    sim = plate.run(t_final = 200*60,
                    dt = 2.,
                    params = params)
    print("run simulation")

    ## plotting
    print("run simulation")
    plate.plot_simulation(sim, 10)

    print("plot conc dist")
    S = plate.get_all_species()
    plate.plot_conc_distribution(sim, S, 10)
    
    plate.compare_species(sim, 10)

main()