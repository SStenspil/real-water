import pybel as pb
import molecule as m
import openbabel as ob
import numpy as np
import os, sys

"""
Problems:

Perhaps rotate the other hydrogen afterwards?

"""


#read file 
mol_old = pb.readfile("xyz","w6_copy.xyz").next()
w = open('w6_copy.xyz')
lines_after_2 = w.readlines()[2:]
# defines forcefield MMFF94
forcefield = ob.OBForceField.FindForceField("MMFF94")


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

# funtion for rotatin one hydrogen atom
def Rot(RotationPoint, p1, p2, theta):
    x = RotationPoint[0] 
    y = RotationPoint[1]
    z = RotationPoint[2]
    a = p1[0] #arbitrary point, oxygen
    b = p1[1]
    c = p1[2]
    u = p2[0] - a # direction vector
    v = p2[1] - b 
    w = p2[2] - c
    L = u**2 + v**2 + w**2 # length of direction vector
    npx = ((a*(v**2 +w**2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-np.cos(theta))+L*x*np.cos(theta)+np.sqrt(L)*(-c*v+b*w-w*y+v*z)*np.sin(theta))/L
    npy = ((b*(u**2 +w**2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-np.cos(theta))+L*y*np.cos(theta)+np.sqrt(L)*(c*u-a*w+w*x-u*z)*np.sin(theta))/L
    npz = ((c*(u**2 +v**2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-np.cos(theta))+L*z*np.cos(theta)+np.sqrt(L)*(-b*u+a*v-v*x+u*y)*np.sin(theta))/L
    return npx, npy, npz 

c = 0  
for j in range(500):
    
    theta = np.random.uniform(0.0,np.pi, 6)
    for n in range(len(theta)):
        for i in range(0,18,3): # 18 = number of atoms
            point = [x[1+i],y[1+i],z[1+i]] # hydrogen to be rotated, first no 1
            p1 = [x[0+i],y[0+i],z[0+i]] # oxygen, first no. 0
            p2 = [x[2+i],y[2+i],z[2+i]] # other hydrogen no 2
            new_point_x,new_point_y,new_point_z = Rot(point,p1,p2,theta[n])
            #print new_point_x
            x[1+i]=new_point_x
            y[1+i]=new_point_y
            z[1+i]=new_point_z
    #write new coordinates to new temporary file
    xyz = open('water_temp.xyz','w')

    noatoms = str(len(x))
    xyz.write('18\ntest\n')
    for n in range(len(x)):
        xyz.write('%s         %8s        %8s       %8s\n' % (atom[n],x[n],y[n],z[n]))
    xyz.close()


    mol_new = pb.readfile("xyz","water_temp.xyz").next()
    #calculate MMFF94 force field energies
    energy_old = m.get_energy(mol_old)
    energy_new = m.get_energy(mol_new)
    
    #overwrite first file if the energy is lower
    if energy_new < energy_old:
        os.system("mv water_temp.xyz w6_copy.xyz")
        print "file overwritten"
    if j%5 == 0 :
        c += 1
        print c
    
#constrants...
#constraints = ob.OBFFConstraints()
#constraints.AddDistanceConstraint(1, 10, 0.9)       # Angstrom
#constraints.AddAngleConstraint(1, 2, 3, 120.0)      # Degrees
#constraints.AddTorsionConstraint(1, 2, 3, 4, 180.0) # Degrees


"""
#write new coordinates to new temporary file
xyz = open('water_temp.xyz','w')

noatoms = str(len(x))
xyz.write('18\ntest\n')
for n in range(len(x)):
    xyz.write('%s         %8s        %8s       %8s\n' % (atom[n],x[n],y[n],z[n]))
xyz.close()


mol_new = pb.readfile("xyz","water_temp.xyz").next()
#calculate MMFF94 force field energies
energy_old = m.get_energy(mol_old)
energy_new = m.get_energy(mol_new)

#overwrite first file if the energy is lower
if energy_new < energy_old:
    os.system("mv water_temp.xyz w6_copy.xyz")
    print "file overwritten"

print energy_old
print energy_new
"""
