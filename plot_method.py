import matplotlib.pyplot as plt
import numpy as np

alle = open("energy_all.csv",'r')
#dobrot = dobrot.readlines()[:250]
hydrogen = open("energy_hydrogen.csv", 'r')
#monte = monte.readlines()[:250]
#dihedral = open("energy_dihedral", 'r')
alle_hastings = open("energy_hastings_all.csv", 'r')

hydrogen_hastings = open("energy_hastings.csv", 'r')

step = []
#step2 = []
e_alle = []
e_hydrogen = []
e_alle_hastings = []
e_hydrogen_hastings = []

for line in alle:
    line = line.split(',')
    step.append(line[0])
    e_alle.append(line[1])
for line in hydrogen:
    line = line.split(',')
    e_hydrogen.append(line[1])
for line in alle_hastings:
    line = line.split(',')
    e_alle_hastings.append(line[1])
for line in hydrogen_hastings:
    #step2.append(line[0])
    line = line.split(',')
    e_hydrogen_hastings.append(line[1])

#Trying to fic x and y problem
step = np.array(step)
e_alle = np.array(e_alle)
e_hydrogen = np.array(e_hydrogen)
e_alle_hastings = np.array(e_alle_hastings)
e_hydrogen_hastings = np.array(e_hydrogen_hastings)


plt.plot(step,e_alle, 'r-',label='Random rotation of all atoms')
plt.plot(step,e_hydrogen, 'r.',label='Random rotation of all hydrogens')
plt.plot(step,e_alle_hastings, 'b-',label='Random rotation of all atoms with Hastings')
plt.plot(step,e_hydrogen_hastings, 'b.',label='Random rotation of all hydrogens with Hastings')

#plt.plot(smonte, emonte, 'g-', label='Rotation of random atom')
#plt.plot(shast,ehast, 'k-', label = 'Hastings')

#plt.ylim(ymax=40, ymin = -55)
#plt.xlim(xmax=200, xmin= 0) 
plt.legend(loc='upper right')
plt.xlabel('Steps')
plt.ylabel('Energy [kcal/mol]')

plt.savefig('plot-metoder.png')


