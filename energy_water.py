import pybel as pb
import molecule as m
import openbabel as ob
import numpy as np

#read file 
mol = pb.readfile("xyz","w6.xyz").next()
w = open('w6.xyz')
# defines forcefield MMFF94
forcefield = ob.OBForceField.FindForceField("MMFF94")

# Number of dihedral angles corresponding to water molecules
no_dihedral = 1


#calculate MMFF94 force field energy
energy = m.get_energy(mol)


for line in w:
    line = line.split()
    


mols =
print mols

"""
# dihedral angles
angles = np.random.uniform(0.0, 120.0, no_dihedral)
m.set_dihedral(angles, mol, forcefield)

m.save_molecule(water2.xyz)
"""



# rotate H1 about O-H2 bond

def R(theta, u):
    return [[np.cos(theta) + u[0]**2 * (1-np.cos(theta)), 
             u[0] * u[1] * (1-np.cos(theta)) - u[2] * np.sin(theta), 
             u[0] * u[2] * (1 - np.cos(theta)) + u[1] * np.sin(theta)],
            [u[0] * u[1] * (1-np.cos(theta)) + u[2] * np.sin(theta),
             np.cos(theta) + u[1]**2 * (1-np.cos(theta)),
             u[1] * u[2] * (1 - np.cos(theta)) - u[0] * np.sin(theta)],
            [u[0] * u[2] * (1-np.cos(theta)) - u[1] * np.sin(theta),
             u[1] * u[2] * (1-np.cos(theta)) + u[0] * np.sin(theta),
             np.cos(theta) + u[2]**2 * (1-np.cos(theta))]]

def Rotate(pointToRotate, point1, point2, theta):


    u= []
    squaredSum = 0
    for i,f in zip(point1, point2):
        u.append(f-i)
        squaredSum += (f-i) **2

    u = [i/squaredSum for i in u]

    r = R(theta, u)
    rotated = []

    for i in range(3):
        rotated.append(round(sum([r[j][i] * pointToRotate[j] for j in range(3)])))

    return rotated


point = []
p1 = []
p2 = []

#print Rotate(point, p1, p2, np.pi)



