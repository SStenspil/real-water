import pybel as pb
import molecule as m
import openbabel as ob
import numpy as np

"""
Problems:
Not sure that the rotation is correct --> understand rotate func
Save as a molecule/object to calculate energy (problem reading the first line in new file) - cannot read the same file that is made in one go
Perhaps rotate the other hydrogen afterwards?

"""


#read file 
mol = pb.readfile("xyz","w6.xyz").next()
w = open('w6.xyz')
lines_after_2 = w.readlines()[2:]
# defines forcefield MMFF94
forcefield = ob.OBForceField.FindForceField("MMFF94")

# Number of dihedral angles corresponding to water molecules
#no_dihedral = 1


#calculate MMFF94 force field energy
energy = m.get_energy(mol)

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
# dihedral angles
angles = np.random.uniform(0.0, 120.0, no_dihedral)
m.set_dihedral(angles, mol, forcefield)

m.save_molecule(water2.xyz)
"""

# function to rotate H1 about O-H2 bond

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

def Rotate(pointToRotate, point1, point2, theta): #What does this do?????

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



#Rotation of 6/multiple hydrogen atoms lidt bedre!!
theta = np.random.uniform(0.0,np.pi*0.66, 6)
print theta
for n in range(len(theta)):
    for i in range(0,18,3): # 18 = number of atoms
        point = [x[1+i],y[1+i],z[1+i]] # hydrogen to be rotated first no 1
        p1 = [x[0+i],y[0+i],z[0+1]] # oxygen first no. 0
        p2 = [x[2+i],y[2+i],z[2+i]] # other hydrogen no 2
        new_point = Rotate(point, p1, p2, theta[n])
        x[1+i]=new_point[0]
        y[1+i]=new_point[1]
        z[1+i]=new_point[2]
#constrants...
constraints = ob.OBFFConstraints()
constraints.AddDistanceConstraint(1, 10, 0.9)       # Angstrom
constraints.AddAngleConstraint(1, 2, 3, 120.0)      # Degrees
constraints.AddTorsionConstraint(1, 2, 3, 4, 180.0) # Degrees

"""        


point = [x[1],y[1],z[1]] # hydrogen to be rotated
p1 = [x[0],y[0],z[0]] # oxygen
p2 = [x[2],y[2],z[2]] # other hydrogen


new = Rotate(point,p1,p2,np.pi*2)
x[1] = new[0]
y[1] = new[1]
z[1] = new[2]


"""
#write new coordinates to new file
xyz = open('new_water4.xyz','w')

noatoms = str(len(x))
xyz.write('18\ntest\n')
for n in range(len(x)):
    xyz.write('%s         %8s        %8s       %8s\n' % (atom[n],x[n],y[n],z[n]))


#read new molecule/object and calculate enrgy 
mol_new = pb.readfile("xyz","new_water.xyz").next()
energy_new = m.get_energy(mol_new)

print energy, energy_new




