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
    axs[0].set_title('Dependence of araC promoter controlled motility on arabinose concentration')
    axs[0].plot(conc_arab, hf.leaky_hill(conc_arab, K=0.6, lam=2, max=3.96, min=1.68))
    #n = 1 ± 0.6, K = 0.6 ± 0.4M #normalised activity from 0 to 1
    axs[0].set_xlabel('arabinose concentration (M)')
    axs[0].set_ylabel('Swim velocity (mm/min)') #?
    axs[1].set_title('Dependence of luxR promoter controlled motility on AHL concentration')
    axs[1].plot(conc_ahl, hf.leaky_hill(conc_ahl, K=37.5e-9, lam=2, max=3.96, min=1.68)) #K 25 - 50 nM AHL = 
    axs[1].set_xlabel('AHL concentration (M)')
    axs[1].set_ylabel('Swim velocity (mm/min)') #?
    #print(hf.leaky_hill(conc_ahl, K=37.5e-9, lam=2, max=3.96, min=1.68))
    #print(hf.leaky_hill(conc_arab, K=0.6, lam=2, max=3.96, min=1.68))
    ratio = (9.12/5.7)*10**7
    print(str(ratio)+'T : A')
    fig.show()
main()