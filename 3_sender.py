# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 15:59:31 2021

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
    rc = 1E-2       # growth rate of X on N
    Dc = 1E-5       # cell diffusion coefficient
    w = 1
    Da = 0.03
    rho_A = 0.1       # production rate of AHL
    Dt = 0.01     #difusion rate of target

    environment_size = (60, 60)
    plate = Plate(environment_size)

    ## add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        n = D * hf.ficks(species['N'], w) - rho_n * species['N'] * (species['R'])
        return n
    N.set_behaviour(N_behaviour)
    plate.add_species(N)

    ## add receiver strain to the plate
    U_R = np.zeros(environment_size)
    for i in np.linspace(0, 59, environment_size[0]):
        if (i == 1) or (i == 30) or (i == 59):
            U_R[15, int(i)] = 0.001
    R = Species("R", U_R)
    def R_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        
        r = Dc * hf.leaky_hill(s=species['T'], K=1, lam=2, max=1e2, min=1) * hf.ficks(species['R'], w) + rc * species['N'] * species['R']
        #r = Dc * hf.ficks(species['R'], w) + rc * species['N'] * species['R']
        #########
        return r
    R.set_behaviour(R_behaviour)
    plate.add_species(R)

    ## add target to plate
    U_T = np.zeros(environment_size)
    grad_values = np.logspace(-4, 5, environment_size[0])
    for idx, value in enumerate(grad_values):
        U_T[idx,:] = value

    T = Species("T", U_T)
    def T_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        #t = Dt * hf.ficks(species['T'], w)
        t = 0
        return t
    T.set_behaviour(T_behaviour)
    plate.add_species(T)

    ## add AHL to plate
    U_A = np.zeros(environment_size)
    #grad_values = np.logspace(-4, 5, environment_size[0])
    #for idx, value in enumerate(grad_values):
        #U_A[idx,:] = value
    A = Species("A", U_A)
    def A_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da, Dt = params
        #D, rho_n, Dc, rc, w, rho_A, Da = params
        #a = Da * hf.ficks(species['A'], w)
        a = Da * hf.ficks(species['A'], w) + rho_A * species['R']
        #a = 0
        return a
    A.set_behaviour(A_behaviour)
    plate.add_species(A)

    #plate.plot_plate()
    #print("plotted first figure")

    ## run the experiment
    params = (D, rho_n, Dc, rc, w, rho_A, Da, Dt)
    sim = plate.run(t_final = 150*60,
                    dt = 2.,
                    params = params)
    print("run simulation")

    ## plotting
    plate.plot_simulation(sim, 10)

    print("final plot")

main()