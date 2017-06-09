import openbabel as ob
import pybel as pb
import numpy as np

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

