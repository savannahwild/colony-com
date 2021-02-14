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

    environment_size = (20, 20)
    plate = Plate(environment_size)

    ## add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        n = D * helper_functions.ficks(species['N'], w) - rho_n * species['N'] * species['X']
        return n
    N.set_behaviour(N_behaviour)
    plate.add_species(N)

    ## add one strain to the plate
    U_X = np.zeros(environment_size)
    X_pos = [[10., 10.]]
    X_radius = 0.75
    X_coordinates = helper_functions.get_node_coordinates(X_pos,
                                                          X_radius,
                                                          environment_size[0],
                                                          environment_size[1],
                                                          w)
    rows = X_coordinates[:, 0]
    cols = X_coordinates[:, 1]
    U_X[rows, cols] = 0.001
    strain = Species("X", U_X)
    def X_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        x = Dc * helper_functions.ficks(species['X'], w) + rc * species['N'] * species['X']
        return x
    strain.set_behaviour(X_behaviour)
    plate.add_species(strain)

    ## run the experiment
    params = (D, rho_n, Dc, rc, w)
    sim = plate.run(t_final = 1000,
                    dt = .1,
                    params = params)

    ## plotting
    tp = np.arange(0, 18, 2)
    plate.plot_simulation(sim, tp)

main()