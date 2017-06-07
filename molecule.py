import openbabel as ob
import pybel as pb
import numpy as np

#global mol
#global forcefield


def get_energy(mol):
    """ Returns the GAFF energy of the molecule in kcal/mol.
    """

#    global mol

    FF = ob.OBForceField.FindForceField("MMFF94")
    FF.Setup(mol.OBMol)

    return FF.Energy()

def choose_atoms(): #some kind of random but systematic choices of water molecules
    #for i in range(10):
        random = np.random.randint(2)
        O1 = np.random.randint(6)*3 #random O
        if random == 1:
            H1 = O1+1
            H2 = O1+2
        else:
            H1 = O1+2
            H2 = O1+1
        O2 = np.random.randint(6)*3
        while O1 == O2:
            O2 = np.random.randint(6)*3
        
        return O1, H1, H2, O2

#calculate dihedral angle https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
# angle in radians
def find_angle(n_O1, n_H1, n_H2, n_O2, x, y, z):

    O1 = np.array([x[n_O1], y[n_O1], z[n_O1]])
    H1 = np.array([x[n_H1], y[n_H1], z[n_H1]])
    H2 = np.array([x[n_H2], y[n_H2], z[n_H2]])
    O2 = np.array([x[n_O2], y[n_O2], z[n_O2]])

    #vector coordinates by vector substraction
    b1 = O1 - H1
    b2 = H2 - O1
    b3 = O2 - H2

    # normalvector
    c1 = np.cross(b1,b2)
    c2 = np.cross(b2,b3)

    n1 = c1/np.linalg.norm(c1)
    n2 = c2/np.linalg.norm(c2)

    b2_norm = b2/np.linalg.norm(b2)

    m1 = np.cross(n1,b2_norm)

    #coordinates on orthonormal frame by n1, b2_norm and m1
    x = np.dot(n1,n2)
    y = np.dot(m1,n2)

    #dihedral angle
    angle = np.arctan2(x,y)
    
    return angle

def set_dihedral(mol, forcefield,x,y,z):

#    global mol
#    global forcefield

    constraints = ob.OBFFConstraints()
    # constraints.AddDistanceConstraint(1, 10, 3.4)       # Angstroms
    # constraints.AddAngleConstraint(1, 2, 3, 120.0)      # Degrees
    # constraints.AddTorsionConstraint(1, 2, 3, 4, 180.0) # Degrees
    
    for i in range(2):
        n_O1, n_H1, n_H2, n_O2 = choose_atoms() #random oxygens
        angle = find_angle(n_O1, n_H1, n_H2, n_O2 ,x,y,z) #initial dihedral angle

        a = mol.OBMol.GetAtom(n_H1)
        b = mol.OBMol.GetAtom(n_O1)
        c = mol.OBMol.GetAtom(n_H2)
        d = mol.OBMol.GetAtom(n_O2)

        mol.OBMol.SetTorsion(a, b, c, d, angle)
        anglep = mol.OBMol.GetTorsion(a, b, c, d)
        

        # Define constraint. fjernet fordi "ai" ikke giver mening
        # constraints.AddTorsionConstraint(ai, bi, ci, di, anglep) # Degrees
        
        print "rotate"

    # Setup the force field with the constraints
    # forcefield = ob.OBForceField.FindForceField("MMFF94")
    forcefield.Setup(mol.OBMol)
    forcefield.SetConstraints(constraints)
    forcefield.EnableCutOff(True) # VDW cutoff
    forcefield.SetElectrostaticCutOff(0) # Remove electrostatics

    # Use forcefield to reoptimize bondlengths+angles
    find_local_min(mol, forcefield)
    return mol


def find_local_min(mol, forcefield):


#    global mol
#    global forcefield


    # tested both 5 and 25 steps. As far as I remember, 5 steps were enough
    forcefield.SteepestDescent(5)
    forcefield.GetCoordinates(mol.OBMol)
    return mol #optimeret mol


def save_molecule(filename):
    """ Save the molecule as <filename> in XYZ format.
    """

 #   global mol

    mol.write('w6.xyz', filename, overwrite="True")


if __name__ == '__main__':

    #generate_chain(6)

    #print get_energy()

    # save_molecule('carbon_man1.xyz')

    #set_dihedral([180.0, -180.0, 180.0]) throw away write something else

    print get_energy()

    # save_molecule('carbon_man2.xyz')

    find_local_min()

    print get_energy()

    # save_molecule('carbon_man3.xyz')
