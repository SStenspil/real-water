import openbabel as ob
import pybel as pb
import numpy as np
import os, sys

global mol
global forcefield


def get_energy(mol):
    """ Returns the GAFF energy of the molecule in kcal/mol.
    """

#    global mol

    FF = ob.OBForceField.FindForceField("MMFF94")
    FF.Setup(mol.OBMol)

    return FF.Energy()


def find_local_min(mol, forcefield):
    """
    """


    # tested both 5 and 25 steps. As far as I remember, 5 steps were enough
    forcefield.SteepestDescent(5)
    forcefield.GetCoordinates(mol.OBMol)
    return mol #optimeret mol


#own functions


# funtion for rotating one hydrogen atom
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

def perturb(x,y,z, atom, no_atom):
    theta = np.random.uniform(0.0,np.pi*2.0)
    if no_atom % 3 == 0:
        point = [x[no_atom],y[no_atom],z[no_atom]] #atom to be rotated
        p1 = [x[no_atom+1],y[no_atom+1],z[no_atom+1]]
        p2 = [x[no_atom+2],y[no_atom+2],z[no_atom+2]]
                
    if no_atom % 3 == 1:
        point = [x[no_atom],y[no_atom],z[no_atom]] #atom to be rotated
        p1 = [x[no_atom-1],y[no_atom-1],z[no_atom-1]]
        p2 = [x[no_atom+1],y[no_atom+1],z[no_atom+1]]
    if no_atom % 3 == 2:
        point = [x[no_atom],y[no_atom],z[no_atom]] #atom to be rotated
        p1 = [x[no_atom-2],y[no_atom-2],z[no_atom-2]]
        p2 = [x[no_atom-1],y[no_atom-1],z[no_atom-1]]
            
    new_point_x,new_point_y,new_point_z = Rot(point,p1,p2,theta)
    x[no_atom]=new_point_x
    y[no_atom]=new_point_y
    z[no_atom]=new_point_z
    #write new coordinates to new temporary file
    xyz = open('water_temp.xyz','w')

    noatoms = str(len(x))
    xyz.write('18\ntest\n')
    for n in range(len(x)):
        xyz.write('%s         %8s        %8s       %8s\n' % (atom[n],x[n],y[n],z[n]))
    xyz.close()        
            
    return x, y, z


def random_konf():
#for i in range(1):
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

    for j in range(1): #Does it even have to be more than once?
        
        theta = np.random.uniform(0.0,np.pi*2.0, 6)
        theta2 = np.random.uniform(0.0,np.pi*2.0, 6)
        for n in range(len(theta)):
            for i in range(0,18,3): # 18 = number of atoms
                point = [x[1+i],y[1+i],z[1+i]] # 1st hydrogen to be rotated, first no 1
                p1 = [x[0+i],y[0+i],z[0+i]] # oxygen, first no. 0
                p2 = [x[2+i],y[2+i],z[2+i]] # other hydrogen no 2
                new_point_x,new_point_y,new_point_z = Rot(point,p1,p2,theta[n])
                #print new_point_x
                x[1+i]=new_point_x
                y[1+i]=new_point_y
                z[1+i]=new_point_z
                point = [x[2+i],y[2+i],z[2+i]] #2nd hydrogen to be rotated
                p2 = [x[1+i],y[1+i],z[1+i]]
                new_point_x,new_point_y,new_point_z = Rot(point,p1,p2,theta2[n])
                x[2+i]=new_point_x
                y[2+i]=new_point_y
                z[2+i]=new_point_z
            
        #write new coordinates to new temporary file
        xyz = open('water_HE.xyz','w')

        noatoms = str(len(x))
        xyz.write('18\ntest\n')
        for n in range(len(x)):
            xyz.write('%s         %8s        %8s       %8s\n' % (atom[n],x[n],y[n],z[n]))
        xyz.close()

        mol_new = pb.readfile("xyz","water_HE.xyz").next()
        #calculate MMFF94 force field energies
        energy_old = get_energy(mol_old)
        energy_new = get_energy(mol_new)

              
        
        
        return energy_new


def hastings_all(n_gen,T,R,c,energy_list):
    percent = -10 #percentage counter
    
    os.system("cp water_HE.xyz water_HE_temp.xyz")
    
    for j in range(n_gen):
        #read file 
        
        if j == 0:
            w = open('water_HE.xyz')
            mol_old = pb.readfile("xyz","water_HE.xyz").next()
        
        elif j > 0:
            w = open('water_HE_temp.xyz')
            mol_old = pb.readfile("xyz","water_HE_temp.xyz").next()
            
        
        
        lines_after_2 = w.readlines()[2:]
        # defines forcefield MMFF94
        forcefield = ob.OBForceField.FindForceField("MMFF94")
        no_atom = np.random.randint(18)
        energy_old = get_energy(mol_old)
        

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
        
        energy_list[j] += energy_old # for data plotting
        
        theta = np.random.uniform(0.0,np.pi*2.0, 6)
        
        x, y, z= perturb(x,y,z, atom, no_atom) #perturbation in first H2O
        mol_new = pb.readfile("xyz","water_temp.xyz").next()
        
        #calculate MMFF94 force field energies
        energy_old = get_energy(mol_old)
        energy_new = get_energy(mol_new)
        
        #overwrite first file if the energy is lower
        dE = energy_new - energy_old
        #print energy_new, energy_old, dE
        randomtal = np.random.random()
        
        #print np.exp(- np.absolute(dE)/1.0*R*T), randomtal
        
        if dE <= 0:
            os.system("mv water_temp.xyz water_HE_temp.xyz")

            #print "rotation:", energy_new
            #print "Accepted", j
        elif np.exp(- np.absolute(dE)/1.0*R*T) >= randomtal: #np.random.random():  downscale energy difference
        #elif np.random.choice(2, 1, p=[1-float(np.exp(-dE/R*T)), float(np.exp(-dE/R*T)])) == 1:
            #x, y, z = m.perturb(x, y, z, atom, no_atom)
            os.system("mv water_temp.xyz water_HE_temp.xyz")
            #print "Hastings"
            c += 1
        #else:
            #print "nothing changed"
        
        #percentage counter
        if j % (n_gen/10) == 0:
            percent += 10
            print percent
            
    return energy_list
        
        



def hastings_hydrogen(n_gen,T,R,c,energy_list):
    percent = -10 #percentage counter
    
    os.system("cp water_HE.xyz water_HE_temp.xyz")
    
    for j in range(n_gen):
        #read file 
        
        if j == 0:
            w = open('water_HE.xyz')
            mol_old = pb.readfile("xyz","water_HE.xyz").next()
        
        elif j > 0:
            w = open('water_HE_temp.xyz')
            mol_old = pb.readfile("xyz","water_HE_temp.xyz").next()
    
        #for i in range(0,18,3):
            
        lines_after_2 = w.readlines()[2:]
        # defines forcefield MMFF94
        forcefield = ob.OBForceField.FindForceField("MMFF94")
        no_atom = np.random.randint(18)
        while no_atom % 3 == 0:
            no_atom = np.random.randint(18)
            
        energy_old = get_energy(mol_old)
        
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
            
            #f.write flyttet her fra
            
        x,y,z = perturb(x,y,z,atom,no_atom)
       
        #New energy write approach
        energy_list[j] += energy_old
        
        mol_new = pb.readfile("xyz","water_temp.xyz").next()
        
        #calculate MMFF94 force field energies
        energy_old = get_energy(mol_old)
        energy_new = get_energy(mol_new)
        
        #overwrite first file if the energy is lower
        dE = energy_new - energy_old

        randomtal = np.random.random()

        if dE <= 0:
            os.system("mv water_temp.xyz water_HE_temp.xyz")

        elif np.exp(- np.absolute(dE)/1.0*R*T) >= randomtal: 
            os.system("mv water_temp.xyz water_HE_temp.xyz")
 
    
        if j % (n_gen/10) == 0:
            percent += 10
            print percent
                
    return energy_list

#DONE
def rotation_hydrogen(n_gen,T,R,c,energy_list):
    percent = -10 #percentage counter
    os.system("cp water_HE.xyz water_HE_temp.xyz")
    for j in range(n_gen):   
        #read file 
        
        if j == 0:
            w = open('water_HE.xyz')
            mol_old = pb.readfile("xyz","water_HE.xyz").next()
        
        elif j > 0:
            w = open('water_HE_temp.xyz')
            mol_old = pb.readfile("xyz","water_HE_temp.xyz").next()
        
        no_atom = np.random.randint(18)
        while no_atom % 3 == 0:
            no_atom = np.random.randint(18)

        #read file 
        lines_after_2 = w.readlines()[2:]
        # defines forcefield MMFF94
        forcefield = ob.OBForceField.FindForceField("MMFF94")
        energy_old = m.get_energy(mol_old)
        
        energy_list[j] += energy_old
        
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
            
        x,y,z = perturb(x,y,z,atom,no_atom)
        #write new coordinates to new temporary file
        xyz = open('water_temp.xyz','w')

        mol_new = pb.readfile("xyz","water_temp.xyz").next()
        #calculate MMFF94 force field energies
        energy_old = m.get_energy(mol_old)
        energy_new = m.get_energy(mol_new)
        
        #overwrite first file if the energy is lower
        if energy_new < energy_old:
            os.system("mv water_temp.xyz water_HE_temp.xyz")
            print j, "file overwritten:", energy_old, energy_new
        
        if j % (n_gen/10) == 0:
            percent += 10
            print percent

    return energy_list


def rotation_all(n_gen,T,R,c,energy_list):
    percent = -10
    
    for j in range(n_gen):
        #read file 
        if j == 0:
            w = open('water_HE.xyz')
            mol_old = pb.readfile("xyz","water_HE.xyz").next()
        
        elif j > 0:
            w = open('water_HE_temp.xyz')
            mol_old = pb.readfile("xyz","water_HE_temp.xyz").next()
            
        lines_after_2 = w.readlines()[2:]
        
        no_atom = np.random.randint(18)
        energy_old = get_energy(mol_old)
        
        energy_list[j] += energy_old
        
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
        
        
        
        theta = np.random.uniform(0.0,np.pi*2.0, 6)
        
        x, y, z= perturb(x,y,z, atom, no_atom) #perturbation in first H2O
         
        
        mol_new = pb.readfile("xyz","water_temp.xyz").next()
        #calculate MMFF94 force field energies
        energy_old = get_energy(mol_old)
        energy_new = get_energy(mol_new)
        
            #overwrite first file if the energy is lower
        if energy_new < energy_old:
            os.system("mv water_temp.xyz water_HE_temp.xyz")
        
        if j % (n_gen/10) == 0:
            percent += 10
            print percent
    return energy_list

