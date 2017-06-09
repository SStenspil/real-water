import matplotlib.pyplot as plt

dobrot = open("energy_rot.csv",'r')
#dobrot = dobrot.readlines()[:250]
monte = open("energy_newmonte.csv", 'r')
#monte = monte.readlines()[:250]
#dihedral = open("energy_dihedral", 'r')
hastings = open("energy_hastings.csv", 'r')

hastings_H = open(

sdobrot = []
edobrot = []

smonte = []
emonte = []

shast = []
ehast = []

for line in dobrot:
    line = line.split(',')
    sdobrot.append(line[0])
    edobrot.append(line[1])
for line in monte:
    line = line.split(',')
    smonte.append(line[0])
    emonte.append(line[1])
for line in hastings:
    line = line.split(',')
    shast.append(line[0])
    ehast.append(line[1])


plt.plot(sdobrot,edobrot, 'r-',label='Random rotation of all hydrogens')
plt.plot(smonte, emonte, 'g-', label='Rotation of random atom')
plt.plot(shast,ehast, 'k-', label = 'Hastings')

#plt.ylim(ymax=40, ymin = -55)
#plt.xlim(xmax=200, xmin= 0) 
plt.legend(loc='upper right')
plt.xlabel('Steps')
plt.ylabel('Energy [kcal/mol]')

plt.savefig('methods2.png')


