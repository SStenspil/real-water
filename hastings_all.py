import pybel as pb
import molecule as m
import openbabel as ob
import numpy as np
import os, sys



f = open("energy_hastings.csv", 'w')
T = 300.0
R = 0.00198588 #kcal/mol
c = 0
n_gen = 2000

for j in range(n_gen):
    #read file 
    mol_old = pb.readfile("xyz","water_HE.xyz").next()
    w = open('water_HE.xyz')
    lines_after_2 = w.readlines()[2:]
    # defines forcefield MMFF94
    forcefield = ob.OBForceField.FindForceField("MMFF94")
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
    dE = energy_new - energy_old
    #print energy_new, energy_old, dE
    randomtal = np.random.random()
    
    print np.exp(- np.absolute(dE)/1.0*R*T), randomtal
    if dE <= 0:
        os.system("mv water_temp.xyz water_HE.xyz")
        #print "rotation:", energy_new
        print "Accepted", j
    elif np.exp(- np.absolute(dE)/1.0*R*T) >= randomtal: #np.random.random():  downscale energy difference
    #elif np.random.choice(2, 1, p=[1-float(np.exp(-dE/R*T)), float(np.exp(-dE/R*T)])) == 1:
        #x, y, z = m.perturb(x, y, z, atom, no_atom)
        os.system("mv water_temp.xyz water_HE.xyz")
        print "Hastings"
        c += 1
    else:
        print "nothing changed" 

print c
    
