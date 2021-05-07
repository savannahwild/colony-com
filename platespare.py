# -*- coding: utf-8 -*-
"""
Created on Sun May  2 15:48:21 2021

@author: savan
"""

    def plot_conc_distribution(self, sim, species, timepoints, fig, axs,colour):
        tps = np.linspace(0, sim.shape[3] - 1, timepoints)
        #section=['upper half','lower half']
         
        y=[]
        y1=[0,0]
        y2=[0,0]
        tp = int(tps[5])
        for idx, s in enumerate(self.species):
            if idx==1: 
                for i in range(0, 29): 
                    for j in range(0, 59):
                        y1[0]+= sim[idx, i, j, tp]/3481                   
                for i in range(30, 59): 
                    for j in range(0, 59):
                        y2[0] += sim[idx, i, j, tp]/3481              
                y.append(y1)
                y.append(y2)
                labels=['Sender']
                x=np.arange(1,2,1)
            if idx==2 and s.get_name()=='R':
                for i in range(0, 29): 
                    for j in range(0, 59):
                        y1[1]+=sim[idx, i, j, tp]/3481                   
                for i in range(30, 59): 
                    for j in range(0, 59):
                        y2[1] += sim[idx, i, j, tp]/3481
                labels=['Sender', 'Receiver']
                x=np.arange(1,2,0.5) 
                y.append(y1)
                y.append(y2)     
        
        axs.bar(x-0.06, y[0],width=0.1,color='b',label='upper half')  
        axs.bar(x+0.06,y[1],width=0.1,color='g',label='lower half')
        axs.set_xticks(x)
        #set_horizontalalignment('center')
        axs.set_xticklabels(labels)
        fig.suptitle('Differences in the average concentration of bacteria between the upper and lower half of the plate')
        #axs.legend()
        axs.set_ylim(0,1.5)
        axs.set_xlabel('E.coli strain')
        axs.set_ylabel('Concentration (mm^2/min)')
        plt.style.use('ggplot')
        plt.tight_layout()
        fig.show()        

      
    #not sued
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
        #fig.show()