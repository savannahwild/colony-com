import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.integrate import solve_ivp

#plate = Plate(size)

class Plate: 
    def __init__(self, size):
        self.size = size
        #plate.size = size
        self.species = []
        #plate.species = empty list
        #[X,Y,Z]
    def get_size(self):
        return self.size
        #plate.get_size returns plate.size

    def get_num_species(self):
        return len(self.species)
        #plate.get_num_species returns number of items in plate.species
        
    def get_all_species(self):
        return self.species
        #plate.get_all_species returns plate.species
    
    def get_species_by_name(self, name):
        for s in self.species:
            if s.get_name() == name:
                return s
        else:
            return None
        #plate.get_species_by_name("Y") returns X,Y,Z in plate.species[] where s.get_name == "Y"
        #returns Y
    
    def get_all_species_U(self):
        U = np.zeros((self.get_num_species(), self.size[0], self.size[1]))
        #U is a empty grid, plate area by number of species
        #3 x 20 x 20
        for idx, s in enumerate(self.species): 
            #s in enumerate is (X.name,X.U),(Y.name,Y.U) etc
            U[idx] = s.get_U()
            #U[idx] = idx.__U__()
            #U[X,Y,Z] = (X.U,Y.U,Z.U)
        return U
        #Y.get_all_species_U return grid of 0 or X,Y,Z conc

   #plate.add_species(Y) 
    def add_species(self, new_species):
        self.species.append(new_species)
        #plate.species.append(Y)

    #plate.set_species(X)
    def set_species(self, species): #reset species list?
        self.species = species
        #plate.species = [X]
##
    #plate.model(t?,y?,params)
    def model(self, t, y, params):
        U = y.reshape(self.get_all_species_U().shape)
        #U = 
        dU = np.zeros(U.shape)
        species_dict = {}
        behaviour_dict = {}
        for idx, s in enumerate(self.species):
            species_dict[s.get_name()] = U[idx]
            behaviour_dict[s.get_name()] = s.behaviour
        for idx, s in enumerate(self.species):
            dU[idx] = behaviour_dict[s.get_name()](species_dict, params)
        return dU.flatten()

    def run(self, t_final, dt, params):
        t = np.arange(0, t_final, dt)
        U_init = self.get_all_species_U().flatten()
        sim_ivp = solve_ivp(self.model, [0, t_final], U_init,
                            t_eval=t, args=(params,))
        sim_ivp = sim_ivp.y.reshape(self.get_num_species(),
                                    self.size[0], self.size[1],
                                    int(t_final / dt))
        return sim_ivp

    def plot_simulation(self, sim, timepoints):
        for tp in timepoints:
            fig, axs = plt.subplots(int(np.ceil(len(self.species) / 3)), 3, sharex='all', sharey='all')
            for idx, (ax, s) in enumerate(zip(axs.flatten(), self.species)):
            # for idx, s in enumerate(self.species):
                im = ax.imshow(sim[idx, :, :, tp*600], interpolation="none",
                                     cmap=cm.viridis, vmin=0,
                               vmax=np.max(sim[idx, :, :, :]))
                ax.set_title(s.get_name() + ' hour: ' + str(tp))
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(im, cax=cax, shrink=0.8)
            fig.savefig('fig_hour_' + str(tp) +'.pdf')
            fig.show()

    def plot_plate(self):
        print("plotting plate")
        fig, axs = plt.subplots(int(np.ceil(len(self.species) / 3)), 3, sharex='all', sharey='all')
        for idx, s in enumerate(self.species):
            im = axs[idx].imshow(s.get_U(), interpolation="none", cmap=cm.viridis, vmin=0)
            axs[idx].set_title(s.get_name())
            divider = make_axes_locatable(axs[idx])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax, shrink=0.8)
        fig.savefig('fig.pdf')
        fig.show()
