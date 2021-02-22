import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.integrate import solve_ivp

class Plate:
    def __init__(self, size):
        self.size = size
        self.species = []

    def get_size(self):
        return self.size

    def get_num_species(self):
        return len(self.species)

    def get_all_species(self):
        return self.species


    def get_species_by_name(self, name):
        for s in self.species:
            if s.get_name() == name:
                return s
        else:
            return None
    
    def get_all_species_U(self):
        U = np.zeros((self.get_num_species(), self.size[0], self.size[1]))
        for idx, s in enumerate(self.species):
            U[idx] = s.get_U()
        return U
            
    def add_species(self, new_species):
        self.species.append(new_species)

    def set_species(self, species):
        self.species = species

    def model(self, t, y, params):
        U = y.reshape(self.get_all_species_U().shape)
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
    

    
    def plot_simulation(self, sim, bacteria, timepoints):
        tps = np.linspace(0, sim.shape[3] - 1, timepoints)
        x = []
        y_top = []
        y_bottom = []
        conc_top = 0
        conc_bottom = 0
        #fig, axs = plt.subplots(int(np.ceil(len(self.species) / 3)), 3, sharex='all', sharey='all')
        fig, axs = plt.subplots(int(np.ceil(len(self.species) / 3)), 3)
        for tp in tps:
            tp = int(tp)
            x.append(tp)
            for idx, (ax, s) in enumerate(zip(axs.flatten(), self.species)):
            # for idx, s in enumerate(self.species):
                im = ax.imshow(sim[idx, :, :, tp], interpolation="none",
                                     cmap=cm.viridis, vmin=0,
                               vmax=np.max(sim[idx, :, :, :]))
                ax.set_title(s.get_name() + ' timepoint: ' + str(tp))
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(im, cax=cax, shrink=0.8)
                if s.get_name()== "S":
                    for i in range(0, 29):
                        for j in range(0, 59):
                            conc_top += sim[idx, i, j, tp]
                if s.get_name()== "S":
                    for i in range(30, 59):
                        for j in range(0, 59):
                            conc_bottom += sim[idx, i, j, tp]   
            y_top.append(conc_top/3600)
            y_bottom.append(conc_bottom/3600)
            #if tp == tps[3]:
            axs[2].plot(x, y_top)
            axs[2].plot(x, y_bottom)
            axs[2].set_xlabel("fig_timepoint_")
            axs[2].set_ylabel("concentration of S")
            fig.savefig('fig_timepoint_' + str(tp) +'.pdf')
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
