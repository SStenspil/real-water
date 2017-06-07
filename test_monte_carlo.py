import pybel as pb
import molecule as m
import openbabel as ob
import numpy as np
import os, sys


"""
Problems
"""
f = open("energy_evo5.csv", 'w')

c = 0
n_gen = 20
for j in range(n_gen):
    #read file 
    mol_old = pb.readfile("xyz","water_HE.xyz").next()
    w = open('w6_copy.xyz')
    lines_after_2 = w.readlines()[2:]
    # defines forcefield MMFF94
    forcefield = ob.OBForceField.FindForceField("MMFF94")
    no_atom = np.random.randint(3)
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
    if energy_old > energy_new:
        os.system("mv water_temp.xyz water_HE.xyz")
        print "rotation:", energy_new
        print j
        continue
    
    x, y, z = m.perturb(x,y,z, atom, no_atom+3) #pert 2nd
    mol_new = pb.readfile("xyz","water_temp.xyz").next()
    #calculate MMFF94 force field energies
    energy_old = m.get_energy(mol_old)
    energy_new = m.get_energy(mol_new)
    if energy_old > energy_new:
        os.system("mv water_temp.xyz water_HE.xyz")
        print "rotation:", energy_new
        print j
        continue
            
    x, y, z= m.perturb(x,y,z, atom, no_atom+6)  # pert 3rd
    mol_new = pb.readfile("xyz","water_temp.xyz").next()
    #calculate MMFF94 force field energies
    energy_old = m.get_energy(mol_old)
    energy_new = m.get_energy(mol_new)
    if energy_old > energy_new:
        os.system("mv water_temp.xyz water_HE.xyz")
        print "rotation:", energy_new
        print j
        continue

    x, y, z = m.perturb(x,y,z, atom, no_atom+9) #pert 4
    mol_new = pb.readfile("xyz","water_temp.xyz").next()
    #calculate MMFF94 force field energies
    energy_old = m.get_energy(mol_old)
    energy_new = m.get_energy(mol_new)
    if energy_old > energy_new:
        os.system("mv water_temp.xyz water_HE.xyz")
        print "rotation:", energy_new
        print j
        continue


    x, y, z= m.perturb(x,y,z, atom, no_atom+12) #pert 5
    mol_new = pb.readfile("xyz","water_temp.xyz").next()
    #calculate MMFF94 force field energies
    energy_old = m.get_energy(mol_old)
    energy_new = m.get_energy(mol_new)
    if energy_old > energy_new:
        os.system("mv water_temp.xyz water_HE.xyz")
        print "rotation:", energy_new
        print j
        continue

    x, y, z = m.perturb(x,y,z, atom, no_atom+15) #pert 6
    mol_new = pb.readfile("xyz","water_temp.xyz").next()
    #calculate MMFF94 force field energies
    energy_old = m.get_energy(mol_old)
    energy_new = m.get_energy(mol_new)
    if energy_old > energy_new:
        os.system("mv water_temp.xyz water_HE.xyz")
        print "rotation:", energy_new
        print j
        continue
    
    
    
    #if j % n_gen/10 == 0:
        #c += 10 
        #print c



