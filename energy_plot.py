import matplotlib.pyplot as plt

evo1 = open("energy_evo1.csv","r")
evo2 = open("energy_evo2.csv","r")
evo3 = open("energy_evo3.csv","r")
evo4 = open("energy_evo4.csv","r")
evo5 = open("energy_evo5.csv","r")

steps1 = []
energy1 = []
steps2 = []
energy2 = []
steps3 = []
energy3 = []
steps4 = []
energy4 = []
steps5 = []
energy5 = []

for line in evo1:
    line = line.split(',')
    steps1.append(line[0])
    energy1.append(line[1])
for line in evo2:
    line = line.split(',')
    steps2.append(line[0])
    energy2.append(line[1])
for line in evo3:
    line = line.split(',')
    steps3.append(line[0])
    energy3.append(line[1])
for line in evo4:
    line = line.split(',')
    steps4.append(line[0])
    energy4.append(line[1])
for line in evo5:
    line = line.split(',')
    steps5.append(line[0])
    energy5.append(line[1])



plt.plot(steps1, energy1, 'r-')
plt.plot(steps2, energy2, 'k-')
plt.plot(steps3, energy3, 'b-')
plt.plot(steps4, energy4, 'g-')
plt.plot(steps5, energy5, 'y-')
plt.xlabel("Steps")
plt.ylabel("Energy [kcal/mol]")



plt.savefig("energy_evo2.png")

evo1.close()
evo2.close()
evo3.close()
evo4.close()
evo5.close()



