<<<<<<< Updated upstream
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:41:10 2021

@author: savan
"""

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf

def main():
## experimental parameters
    D = 3E-3
    rho_n = 3E-1     
    rc = 3.5E-2
    Dc = 1E-5       
    w = 1           

    environment_size = (60, 60)
    plate = Plate(environment_size)

## add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        n = D * hf.ficks(species['N'], w) - species['S'] * rho_n * hf.leaky_hill(s=species['N'], K=0.15, lam=1, max=rc, min=0)
        return n
    N.set_behaviour(N_behaviour)
    plate.add_species(N)

## add sender strain to the plate
    U_S = np.zeros(environment_size)
    for i in np.linspace(30, 30, 1):
        U_S[int(i), int(i)] = 0.001
    S = Species("S", U_S)
    def S_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        s = Dc * hf.ficks(species['S'], w) + species['S'] * hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
        return s
    S.set_behaviour(S_behaviour)
    plate.add_species(S)

    U_H = np.zeros(environment_size)
    for i in np.linspace(30,59, 30):
        for j in np.linspace(0, 59, 60):
              if U_S[int(j), int(i)] > 0:
                  U_H[int(j),(int(i)-30)] = 1

   # for i in np.linspace(0,59, environment_size[0]):
    #    for i in np.linspace(30,59, 30):
     #         if U_S.any() > 0:
      #             ConcS = ConcS + U_S        
                   
    #U_D[i, math.exp(ConcS)] = 1
    H = Species("H", U_H)
    def H_behaviour(species, params):
        D, rho_n, Dc, rc, w = params
        h = 0
        return h
    H.set_behaviour(H_behaviour)
    plate.add_species(H)
    
## run the experiment
    params = (D, rho_n, Dc, rc, w)
    sim = plate.run(t_final = 100*60,
                    dt = 1.,
                    params = params)
    print("run simulation")
    
## plotting
    plate.plot_simulation(sim, 4)
    print("final plot")
    
    #concS = 0
    #for t in sim:
     #   for i in np.linspace(0, 59, environment_size):
      #      for j in np.linspace(0, 59, environment_size):
       #         if U_S[i, j] > 0:
        #            concS = concS + U_S[i, j]
    #plt.plot(t, math.log(concS))
        
=======
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:41:10 2021

@author: savan
"""

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf

def main():
## experimental parameters
    D = 3E-3
    rho_n = 3E-1     
    rc = 3.5E-2
    Dc = 1E-5       
    w = 1           

    environment_size = (60, 60)
    plate = Plate(environment_size)

## add nutrient to the plate
    U_N = np.ones(environment_size)
    N = Species("N", U_N)
    def N_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        n = D * hf.ficks(species['N'], w) - species['S'] * rho_n * hf.leaky_hill(s=species['N'], K=0.15, lam=1, max=rc, min=0)
        return n
    N.set_behaviour(N_behaviour)
    plate.add_species(N)

## add sender strain to the plate
    U_S = np.zeros(environment_size)
    for i in np.linspace(30, 30, 1):
        U_S[int(i), int(i)] = 0.001
    S = Species("S", U_S)
    def S_behaviour(species, params):
        ## unpack params
        D, rho_n, Dc, rc, w = params
        s = Dc * hf.ficks(species['S'], w) + species['S'] * hf.leaky_hill(s=species['N'], K=0.15, lam = 1, max=rc, min=0)
        return s
    S.set_behaviour(S_behaviour)
    plate.add_species(S)

    U_H = np.zeros(environment_size)
    for i in np.linspace(30,59, 30):
        for j in np.linspace(0, 59, 60):
              if U_S[int(j), int(i)] > 0:
                  U_H[int(j),(int(i)-30)] = 1

   # for i in np.linspace(0,59, environment_size[0]):
    #    for i in np.linspace(30,59, 30):
     #         if U_S.any() > 0:
      #             ConcS = ConcS + U_S        
                   
    #U_D[i, math.exp(ConcS)] = 1
    H = Species("H", U_H)
    def H_behaviour(species, params):
        D, rho_n, Dc, rc, w = params
        h = 0
        return h
    H.set_behaviour(H_behaviour)
    plate.add_species(H)
    
## run the experiment
    params = (D, rho_n, Dc, rc, w)
    sim = plate.run(t_final = 100*60,
                    dt = 1.,
                    params = params)
    print("run simulation")
    
## plotting
    plate.plot_simulation(sim, 4)
    print("final plot")
    
    #concS = 0
    #for t in sim:
     #   for i in np.linspace(0, 59, environment_size):
      #      for j in np.linspace(0, 59, environment_size):
       #         if U_S[i, j] > 0:
        #            concS = concS + U_S[i, j]
    #plt.plot(t, math.log(concS))
        
>>>>>>> Stashed changes
main()