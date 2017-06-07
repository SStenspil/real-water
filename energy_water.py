import pybel as pb
import molecule as m
import openbabel as ob
import numpy as np
import matplotlib.pyplot as plt


#read file 
mol = pb.readfile("xyz","w6.xyz").next()
w = open('w6.xyz')
lines_after_2 = w.readlines()[2:] #removes to first lines
# defines forcefield MMFF94
forcefield = ob.OBForceField.FindForceField("MMFF94")

#calculate MMFF94 force field energy
energy = m.get_energy(mol)

no_atoms = 4

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


"""
Run it here
"""

print m.get_energy(mol)

mol = m.set_dihedral(mol, forcefield,x,y,z)

print m.get_energy(mol)



"""
#Its important to call this function before the others as it generates the global variable mol
mol.generate_chain(no_atoms)
# Number of dihedral angles corresponding
# to the molecule alkane chain


# Create a list of random dihedral angles
# between 0.0 and 360.0 degrees
dihedral_list = np.random.uniform(0.0, 360.0, no_dihedral)

mol.set_dihedral(dihedral_list)

energy = mol.get_energy()


mol.save_molecule('butane.xyz')

print "The energy for \omega = ", dihedral_list, "is", energy
"""
