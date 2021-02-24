import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.integrate import solve_ivp
import math

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
    

    
    def plot_simulation(self, sim, timepoints):
        tps = np.linspace(0, sim.shape[3] - 1, timepoints)
        for tp in tps:
            fig, axs = plt.subplots(int(np.ceil(len(self.species) / 3)), 3, sharex='all', sharey='all')
            tp = int(tp)
            for idx, (ax, s) in enumerate(zip(axs.flatten(), self.species)):
                im = ax.imshow(sim[idx, :, :, tp], interpolation="none",
                                     cmap=cm.viridis, vmin=0,
                               vmax=np.max(sim[idx, :, :, :]))
                ax.set_title(s.get_name() + ' timepoint: ' + str(tp))
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(im, cax=cax, shrink=0.8)
            fig.savefig('fig_timepoint_' + str(tp) +'.pdf')
            fig.show()
    
    def plot_conc_distribution(self, sim, species, timepoints, sections):
        tps = np.linspace(0, sim.shape[3] - 1, timepoints)
        fig, axs = plt.subplots(1, len(species))
        newsecs = int(math.sqrt(sections))
        for idx, s in enumerate(species):
            ji = []
            jj = []
            count = 0
            for sectionx in range(0, newsecs-1):
                i1 = int(sectionx*(self.size[0]/newsecs))
                i2 = int(sectionx+(self.size[0]/newsecs)-1)
                ji.append(i1) #[0, 20, 40]
                jj.append(i2) #[19, 39, 59]
            for sectiony in range(0, int(newsecs-1)):    #[0,0, 0, 20, 20, 20, 40, 40, 40, 40]#[19, 39, 59]
                for i in range(ji[sectiony], jj[sectiony]):
                    for jjj in range(0, newsecs-1):
                        x = []
                        y = []
                        for tp in tps:
                            tp = int(tp)
                            x.append(tp)
                            summ = 0
                            for j in range(ji[jjj], jj[jjj]):
                                summ += sim[idx, i, j, tp]
                            y.append(summ)
                        count += 1
                        axs[idx].plot(x, y)
                        print(count)
            axs[idx].set_xlabel("fig_timepoint_")
            axs[idx].set_ylabel("concentration of " + str(species[idx].get_name()))
        fig.show()        
                    
            
        
    def compare_species(self, sim, timepoints):
        tps = np.linspace(0, sim.shape[3] - 1, timepoints)
        fig, axs = plt.subplots(1, 2)
        for idx, s in enumerate(self.species):
            y = []
            x = []
            for pos, tp in enumerate(tps):
                y.append(0)
                tp = int(tp)
                x.append(tp)
                for i in range(0, self.size[0]):
                    for j in range(0, self.size[1]):
                        y[pos] += sim[idx, i, j, tp]
            axs[1].plot(x, y)
            axs[1].set_xlabel("fig_timepoint_")
            axs[1].set_ylabel("concentration")
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
