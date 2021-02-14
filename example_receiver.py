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
    Da = 0.03
    rho_A = 0.1       # production rate of AHL

    environment_size = (50, 50)
    plate = Plate(environment_size)

    ## add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da = params
        n = D * hf.ficks(species['N'], w) - rho_n * species['N'] * species['X']
        return n
    N.set_behaviour(N_behaviour)
    plate.add_species(N)

    ## add one strain to the plate
    U_X = np.zeros(environment_size)
    for i in np.linspace(5, 45, 9):
        for j in np.linspace(5, 45, 9):
            U_X[int(i), int(j)] = 0.001
    strain = Species("X", U_X)
    def X_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da = params
        x = Dc * hf.ficks(species['X'], w) + rc * species['N'] * species['X']
        return x
    strain.set_behaviour(X_behaviour)
    plate.add_species(strain)

    ## add AHL to plate
    U_A = np.zeros(environment_size)
    U_A[25,25] = 1
    A = Species("A", U_A)
    def A_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da = params
        a = Da * hf.ficks(species['A'], w)
        return a
    A.set_behaviour(A_behaviour)
    plate.add_species(A)

    ## add GFP to plate
    U_G = np.zeros(environment_size)
    G = Species("G", U_G)
    def G_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, rho_A, Da = params
        g = Dc * hf.ficks(species['G'], w) + hf.leaky_hill(s = species['A'], K = 1E-3, lam = 2, min = 1E-3, max=1) * species['X']
        return g
    G.set_behaviour(G_behaviour)
    plate.add_species(G)

    ## run the experiment
    params = (D, rho_n, Dc, rc, w, rho_A, Da)
    sim = plate.run(t_final = 28*60,
                    dt = .1,
                    params = params)

    ## plotting
    tp = np.arange(0, 28, 4)
    plate.plot_simulation(sim, tp)

main()