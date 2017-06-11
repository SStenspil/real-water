import pybel as pb
import molecule as m
import openbabel as ob
import numpy as np
import os, sys
import random_konf
import molecule as m


b = open("energy_hastings_all.csv", 'w') #match with name of simulation
T = 300.0
R = 0.00198588 #kcal/mol
c = 0
n_gen = 2000
simulations=100

#This function has to be here. Can't get it to run from molecule.
def save_data(n_gen, energy_list,simulations):
    #saves the average of the simulation energies
    for j in range(n_gen):
        b.write(str(j) + ',' + str(energy_list[j]/simulations))
        b.write('\n')


#RUN SIMULATION
for i in range(simulations):
    print "Simulation", i+1, "of", simulations
    m.random_konf()
    
    #Choose one function to run:
    energy_list = m.hastings_all(n_gen,T,R,c,)
    #energy_list = m.hastings_hydrogen(n_gen,T,R,c,)
    

save_data(n_gen, energy_list, simulations)
