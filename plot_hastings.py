import matplotlib.pyplot as plt


H = open("energy_hastings.csv",'r')

s = []
e = []

for line in H:
    line = line.split(',')
    s.append(line[0])
    e.append(line[1])


plt.plot(s,e, 'r-',label='Random rotation of all hydrogens')
#plt.plot(smonte, emonte, 'g-', label='Monte Carlo')
#plt.plot(sdihedral, edihedral, 'b-', label='Dihedral')
#plt.ylim(ymax=40, ymin = -55)
#plt.xlim(xmax=200, xmin= 0) 
#plt.legend(loc='upper right')
plt.xlabel('Steps')
plt.ylabel('Energy [kcal/mol]')

plt.savefig('hastings.png')
