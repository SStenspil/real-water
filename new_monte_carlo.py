import pybel as pb
import molecule as m
import openbabel as ob
import numpy as np
import os, sys


f = open("energy_newmonte.csv", 'w')

#Variables
forcefield = ob.OBForceField.FindForceField("MMFF94")
c = 0
n_gen = 2000

for j in range(n_gen):
    #read file 
    mol_old = pb.readfile("xyz","water_HE.xyz").next()
    w = open('water_HE.xyz')
    lines_after_2 = w.readlines()[2:]
    
    no_atom = np.random.randint(18)
    energy_old = m.get_energy(mol_old)
    

    atom = []
    x = []
    y = []
    z = []

    for line in lines_after_2:
        line = line.split()
        atom.append(line[0])
        x.append(float(line[1]))
        y.append(float(line[2]))
        z.append(float(line[3]))
    
    f.write(str(j) + ',' + str(energy_old))
    f.write('\n')
    
    
    theta = np.random.uniform(0.0,np.pi*2.0, 6)
    
    x, y, z= m.perturb(x,y,z, atom, no_atom) #perturbation in first H2O
     
    
    mol_new = pb.readfile("xyz","water_temp.xyz").next()
    #calculate MMFF94 force field energies
    energy_old = m.get_energy(mol_old)
    energy_new = m.get_energy(mol_new)
    
        #overwrite first file if the energy is lower
    if energy_new < energy_old:
        os.system("mv water_temp.xyz water_HE.xyz")
        print j, "file overwritten:", energy_old, energy_new
        

    
    
    
    



