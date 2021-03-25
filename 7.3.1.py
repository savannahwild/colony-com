# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 01:16:51 2021

@author: savan
"""

#inducer as point/line , no diffusion

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf
import matplotlib.pyplot as plt

def main():
    ## experimental parameters
    D = 3E-3        # nutrient diffusion coeff (#mm2/min)
    rho_n = 0.3     # consumption rate of nutrients by X
    rc = 3.5E-2       # growth rate of X on N
    Dc = 1E-5       # cell diffusion coefficient
    w = 1
    Da = 0.03
    rho_A = 0.1       # production rate of AHL
    Dt = 0# 6e-2   random difusion rate of target molecule
    #298.2 0.773, 313.2 1.11, m2/s = 6e-2
    environment_size = (59, 59)
    plate = Plate(environment_size)
    fig, axs = plt.subplots(1,3)
    plt.suptitle('Change in concentration of species over time with inducer at a point')
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
    for i in np.linspace(29, 29, 1):
        U_S[int(i), int(i)] = 0.001
    S = Species("S", U_S)
    def S_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        sp = Dc * hf.ficks(species['S'], w)*hf.leaky_hill(s=species['T'], K=0.6, lam=2, max=3.96, min=1.68) + species['S'] * hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
        return sp
    S.set_behaviour(S_behaviour)
    plate.add_species(S)

    ## add target to plate
    U_T = np.zeros(environment_size)
    for i in np.linspace(58,58,3):
        for j in np.linspace(29,29,1):
            U_T[int(i), int(j)] = 10
    T = Species("T", U_T)
    def T_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        t = Dt * hf.ficks(species['T'], w)
        return t
    T.set_behaviour(T_behaviour)
    plate.add_species(T)

    ## run the experiment
    params = (D, rho_n, Dc, rc, w, rho_A, Da, Dt)
    sim = plate.run(t_final = 150*60,
                    dt = 1.,
                    params = params)

    #plate.plot_simulation(sim, 10)
    colour = 'b'
    spec = plate.get_all_species()
    plate.plot_conc_distribution(sim, spec, 10, fig, axs, colour)

main()