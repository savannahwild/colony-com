# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 23:46:11 2021

@author: savan
"""
#show max and min motility

from plate import Plate
from species import Species
import numpy as np
import helper_functions as hf
import matplotlib.pyplot as plt

def main():
    
    ## experimental parameters
    
    conc_arab = np.linspace(0,9.12,100)
    conc_ahl = np.linspace(0,5.7e-7,100)
    
    fig, axs = plt.subplots(1,2)
    plt.style.use('ggplot')
    
    axs[0].set_title('Motility vs arabinose concentration')
    axs[0].plot(conc_arab, hf.leaky_hill(conc_arab, K=0.6, lam=2, max=2.75, min=1))

    #n = 1 ± 0.6, K = 0.6 ± 0.4M #normalised activity from 0 to 1
    axs[0].set_xlabel('arabinose concentration (M)')
    axs[0].set_ylabel('Magnitute of motility induction') #?
    
    axs[1].set_title('Motility vs AHL concentration')
    axs[1].plot(conc_ahl, hf.leaky_hill(conc_ahl, K=40e-9, lam=2, max=2.75, min=1))
    
    axs[1].set_xlabel('AHL concentration (M)')
    axs[1].set_ylabel('Magnitude of motility induction') #?

    ratio = (9.12/5.7)*10**7
    print(str(ratio)+'T : A')
    
    fig.suptitle('Dependence of araC and luxR promoter controlled motility on arabinose and AHL concentration, respectively')
    fig.show()
main()