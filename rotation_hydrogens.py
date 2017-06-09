import pybel as pb
import molecule as m
import openbabel as ob
import numpy as np
import os, sys

"""
Problems:
never returns a lower energy than -49

"""

f = open("energy_rot.csv", 'w')


   
for j in range(2000):        
    #for n in range(len(theta)):
    for i in range(0,18,3): # 18 = number of atoms
        #read file 
        mol_old = pb.readfile("xyz","water_HE.xyz").next()
        w = open('water_HE.xyz')
        lines_after_2 = w.readlines()[2:]
        # defines forcefield MMFF94
        forcefield = ob.OBForceField.FindForceField("MMFF94")
        energy_old = m.get_energy(mol_old)
        
        f.write(str(j) + ',' + str(energy_old))
        f.write('\n')


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
        theta = np.random.uniform(0.0,np.pi*2.0)
        theta2 = np.random.uniform(0.0,np.pi*2.0)
        point = [x[1+i],y[1+i],z[1+i]] # 1st hydrogen to be rotated, first no 1
        p1 = [x[0+i],y[0+i],z[0+i]] # oxygen, first no. 0
        p2 = [x[2+i],y[2+i],z[2+i]] # other hydrogen no 2
        new_point_x,new_point_y,new_point_z = m.Rot(point,p1,p2,theta)
        x[1+i]=new_point_x
        y[1+i]=new_point_y
        z[1+i]=new_point_z
            
        point = [x[2+i],y[2+i],z[2+i]] #2nd hydrogen to be rotated
        p2 = [x[1+i],y[1+i],z[1+i]]
        new_point_x,new_point_y,new_point_z = m.Rot(point,p1,p2,theta2)
        x[2+i]=new_point_x
        y[2+i]=new_point_y
        z[2+i]=new_point_z
        
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
            os.system("mv water_temp.xyz water_HE.xyz")
            print j, "file overwritten:", energy_old, energy_new
    #if j%10 == 0 : 
        #print energy_old, energy_new
   


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
