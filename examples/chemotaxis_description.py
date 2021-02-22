
from plate import Plate
from species import Species
import numpy as np 
import helper_functions as hf


def main():
    D = 3E-3        # nutrient diffusion coeff (#mm2/min)
    rho_n = 0.3     # consumption rate of nutrients by X
    rc = 1E-2       # growth rate of X on N
    Dc = 1E-5       # cell diffusion coefficient
    w = 1
    Da = 0.03
    rho_A = 0.1       # production rate of AHL
    #parameters
    
    environment_size = (20, 20)
    plate = Plate(environment_size)
    #plate.size = (20,20)
    
    #N is nutrient
    U_N = np.ones(environment_size)
    #location of N = 20x20 pixel grid of value one each

    N = Species("N", U_N)
    #N.name = "N"
    #N.U = U_N is grid of one in every pixel
    def N_behaviour(species, params): 
        #species is "N"
        D, rho_n, Dc, rc, w, rho_A, Da = params #defined for N in each location
        #those = params = U_N so = grid of ones
        n = D * hf.ficks(species['N'], w) - rho_n * species['N'] * (species['R'])
        #n = diffusion - consumption
        return n
        #return n which is movement
    N.set_behaviour(N_behaviour)
    #set behaviour as n movement
    plate.add_species(N)
    #plate.species = [N]

    #R is the bacteria
    U_R = np.zeros(environment_size)
    #initial location R = empty grid
    for i in np.linspace(10, 5, 8): #(spaced points on grid)
        U_R[int(i), int(i)] = 0.001
        #conc of R in those coord = 0.001
        #U_R[(0,0),(1,1),etc] = 0.001
    R = Species("R", U_R)
    #R.name = "R"
    #0 in everywhere 0.001 in position diagonal
    def R_behaviour(species, params):
        D, rho_n, Dc, rc, w, rho_A, Da = params
        r = Dc * hf.leaky_hill(s=species['A'], K=1, lam=2, max=1e2, min=1) * hf.ficks(species['R'], w) + rc * species['N'] * species['R']
        #movement = (diffusion coef x hill x ficks) + growth(rate x conc of N x number of R)
        return r
        #return movement
    R.set_behaviour(R_behaviour)
    plate.add_species(R)

    #A is AHL
    U_A = np.zeros(environment_size)
    #initial location A = empty grid
    grad_values = np.logspace(-4, 5, environment_size[0])
    #linear diffusion spaced on log scale (start,end,x-axis)
    for idx, value in enumerate(grad_values):
        U_A[idx,:] = value
        #conc A for pixel in grid = value(log) in grad_values same position
    A = Species("A", U_A)
    def A_behaviour(species, params):
        D, rho_n, Dc, rc, w, rho_A, Da = params
        # a = Da * hf.ficks(species['A'], w)
        a = 0
        return a
        #movement 0, AHL static
    A.set_behaviour(A_behaviour)
    plate.add_species(A)
##
    plate.plot_plate()
    print("plotting plate")
    #run the experiment
    params = (D, rho_n, Dc, rc, w, rho_A, Da)
    sim = plate.run(t_final = 50*60,
                    dt = .1,
                    params = params)
    #final time 50 x 60mins, dt?, parameters
    tp = np.arange(0, 50, 4)
    #time point = 
    plate.plot_simulation(sim, tp)
    print("plotting complete")
main()