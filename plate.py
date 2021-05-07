# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 22:05:57 2021

@author: savan
"""
#Class Plate
#Simulation, growth rate plot, concentration distribution plot
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
    
    def plot_simulation(self, sim, timepoints):
        tps = np.linspace(0, sim.shape[3] - 1, timepoints)
        for tp in tps:
            #if int(tp) == 0 or int(tp) == 4999 or int(tp) == 8999:
                fig, axs = plt.subplots(int(np.ceil(len(self.species)/3)), 3) #sharex='all', sharey='all')
                tp = int(tp)
                for idx, (ax, s) in enumerate(zip(axs.flatten(), self.species)):
                    im = ax.imshow(sim[idx, :, :, tp], interpolation="none",
                                         cmap=cm.viridis, vmin=0,
                                   vmax=np.max(sim[idx, :, :, :]))
                    ax.set_title('Species:'+s.get_name())
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    fig.colorbar(im, cax=cax, shrink=0.8)
                    ax.set_xlabel('Concentration (mm^2/min)')
                    ax.set_ylabel('Concentration (mm^2/min)')
                fig.suptitle('Plate simulation at ' + str(tp) + ' minutes')
                #fig.savefig('fig_timepoint_' + str(tp) +'.pdf')
                plt.style.use('ggplot')
                plt.tight_layout()
                fig.show()
    
    def plot_conc_target(self, sim, species, timepoints, loop):
        tps = np.linspace(0, sim.shape[3] - 1, timepoints)
        xlabels=['1:100','1:10','1:1','10:1','100:1']
        if loop==0:
            fig, axs = plt.subplots(1, 2)
            #xlabels=[]
        if loop>0:
            fig=plt.gcf()
            axs = plt.gcf().get_axes()
        y=[]
        y1=[0,0]
        y2=[0,0]
        y3=[0,0]
        y4=[0,0]
        tp = 5999
        print(tp)
        labels = ['quadrant 1','quadrant 2','quadrant 3','quadrant 4']
        
        for idx, s in enumerate(self.species):
            if idx == 1 or idx==2:
                for i in range(0, 14): 
                    for j in range(0, 59):
                        y1[idx-1]+= sim[idx, i, j, tp]/870.25                   
                for i in range(15, 29): 
                    for j in range(0, 59):
                        y2[idx-1] += sim[idx, i, j, tp]/870.25
                for i in range(30, 44): 
                    for j in range(0, 59):
                        y3[idx-1]+= sim[idx, i, j, tp]/870.25                   
                for i in range(45, 59): 
                    for j in range(0, 59):
                        y4[idx-1] += sim[idx, i, j, tp]/870.25
                y.append(y1)
                y.append(y2)
                y.append(y3)
                y.append(y4)
                x=np.arange(0,5,1)
                #print(y[16])
                #print(y[17])
                #print(y[18])
                #print(y[19])
                axs[idx-1].bar(loop-0.3, y[0],width=0.2,color='b')  
                axs[idx-1].bar(loop-0.1,y[1],width=0.2,color='r')
                axs[idx-1].bar(loop+0.1,y[2],width=0.2,color='g')
                axs[idx-1].bar(loop+0.3,y[3],width=0.2,color='y')
                axs[idx-1].set_xticks(x)
                
                
                axs[idx-1].set_title(s.get_name())
                axs[idx-1].set_xticklabels(['1:100','1:10','1:1','10:1','100:1'])#axs[idx-1].set_ylim(0,3)
                axs[idx-1].set_xlabel('ratio S:R')
                axs[idx-1].set_ylabel('Concentration (mm^2/min)')
       
        #axs[1].set_xticklabels(xlabels)
        #xlabels.append(loop)
       # if len(xlabels) == 5:    
          #  axs[0].set_xticklabels(xlabels)
           # axs[1].set_xticklabels(xlabels)
        fig.legend(labels, title='Plate section', loc='upper left')
        plt.style.use('ggplot')
        plt.tight_layout()
        fig.show()
        
    def compare_species(self, sim, species, timepoints, loop):
        tps = np.linspace(0, sim.shape[3] - 1, timepoints)
        colours = ['b', 'r', 'g', 'y', 'k']
        labels = []
        if loop==0:
            if self.get_num_species() < 4:
                fig, axs = plt.subplots(1, self.get_num_species())
            else:
                fig, axs = plt.subplots(1, 3)
        if loop>0:
            axs = plt.gcf().get_axes()
        for idx, s in enumerate(self.species):
            if idx==0 or idx == 1 or idx==2:
                    y = []
                    x = []
                    for pos, tp in enumerate(tps):
                        y.append(0)
                        tp = int(tp)
                        x.append(tp)
                        for i in range(0, self.size[0]):
                            for j in range(0, self.size[1]):
                                y[pos] += (sim[idx, i, j, tp]/3481)
                    axs[idx].plot(x, y, colours[loop],label=str(sim[idx, 29, 29, 0]))
                    axs[idx].set_xlabel("time (min)")
                    #[idx].set_ylim(0,3)
                    axs[idx].set_ylabel("concentration of " + str(s.get_name()) + " (mm^2/min)")
                    axs[idx].set_title(str(s.get_name()))
                    #axs[idx].set_xlim(0,9000)
                    #axs[idx].legend(title='initial concentration')
        plt.legend([0,0.0025,0.005,0.0075,0.01], title='rho_A')
        plt.tight_layout()
        plt.style.use('ggplot')
        #fig.show()

