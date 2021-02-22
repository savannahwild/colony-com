from plate import Plate
from species import Species
import numpy as np
import helper_functions


def main():
    ## experimental parameters
    D = 3E-3        # nutrient diffusion coeff (#mm2/min)
    rho_n = 0.3     # consumption rate of nutrients by X
    rc = 6E-3       # growth rate of X on N
    Dc = 1E-5       # cell diffusion coefficient
    w = 1
    Da = 0.03
    rho_A = 0.1       # production rate of AHL
    omega = 1E-1

    environment_size = (85, 125)
    plate = Plate(environment_size)

    ## add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, omega = params
        n = D * helper_functions.ficks(species['N'], w) - rho_n * species['N'] * species['X']
        return n
    N.set_behaviour(N_behaviour)
    plate.add_species(N)

    ## add lawn strain to the plate
    U_X = np.ones(environment_size) * 0.001
    strain = Species("X", U_X)
    def X_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, omega = params
        x = Dc * helper_functions.ficks(species['X'], w) + rc * species['N'] * species['X'] - omega * species['B'] * species['X']
        return x
    strain.set_behaviour(X_behaviour)
    plate.add_species(strain)

    ## add bacteriocin to the plate
    U_B = np.zeros(environment_size)
    U_B[[20,40,60], 15] = 0.1
    U_B[[20,40,60], 30] = 0.2
    U_B[[20,40,60], 45] = 0.5
    U_B[[20,40,60], 60] = 1
    U_B[[20,40,60], 75] = 2
    U_B[[20,40,60], 90] = 5
    U_B[[20,40,60], 105] = 10
    B = Species('B', U_B)
    def B_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w, omega = params
        b = D * helper_functions.ficks(species['B'], w) - omega * species['B'] * species['X']
        return b
    B.set_behaviour(B_behaviour)
    plate.add_species(B)

    ## run the experiment
    params = (D, rho_n, Dc, rc, w, omega)
    sim = plate.run(t_final = 20*60,
                    dt = .1,
                    params = params)

    ## plotting
    tp = np.arange(0, 18, 2)
    plate.plot_simulation(sim, tp)

main()