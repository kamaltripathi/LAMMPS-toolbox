import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

NPoly = 204
NSol  = 102
b     = 1.122
qm    = 0.01
qs    = -0.02
atomtype     = 2
bondtype     = 1
angletype    = 1
dihedraltype = 0
impropertype = 0

masstype =[1, 2]
mass = [1.0000, 1.0000]
R = 6.0000

x = []
y = []
z = []

x.append(0.0)
y.append(0.0)
z.append(0.0)


def genpoly():
    genpoly.count = 0
    while genpoly.count < NPoly-1:
        np.random.seed(int(1000000*np.random.uniform()))
        theta = np.random.uniform(0,np.pi)
        phi   = np.random.uniform(0,2*np.pi)
        xtrial = b*np.sin(theta)*np.cos(phi) + x[-1]
        ytrial = b*np.sin(theta)*np.sin(phi) + y[-1]
        ztrial = b*np.cos(theta) + z[-1]
        if (np.sqrt((xtrial)**2 + (ytrial)**2 + (ztrial)**2) < R):
            x.append(xtrial)
            y.append(ytrial)
            z.append(ztrial)
            genpoly.count += 1
        else:
            pass
    
genpoly()

def gensol():
    while genpoly.count < NSol+NPoly-1:
        np.random.seed(int(1000000*np.random.uniform()))
        xtrial = np.random.uniform(-R/2,R/2)
        ytrial = np.random.uniform(-R/2,R/2)
        ztrial = np.random.uniform(-R/2,R/2)
        if (np.sqrt((xtrial)**2 + (ytrial)**2 + (ztrial)**2) < R):
            x.append(xtrial)
            y.append(ytrial)
            z.append(ztrial)
            genpoly.count += 1
        else:
            pass
    
gensol()

bonds = np.zeros(NPoly-1)
for i in range(NPoly-1):
    bonds[i] = np.sqrt((x[i+1] - x[i])**2 + (y[i+1]-y[i])**2 + (z[i+1] - z[i])**2)

r = np.zeros(NPoly+NSol)
for i in range(NPoly+NSol):
    r[i] = np.sqrt((x[i])**2 + (y[i])**2 + (z[i])**2)
    
    
#for i in range(NSol):
#    for j in range(NPoly+NSol):
#        (xs[i]-x[j])**2 + 

fig = plt.figure(figsize=(10,2))
ax1 = fig.add_subplot(141)
plt.plot(x,'o-')
plt.plot(y,'o-')
plt.plot(z,'o-')
   
ax2 = fig.add_subplot(142)
plt.plot(bonds,'o-')    
    
ax3 = fig.add_subplot(143)
plt.plot(r,'o-')

ax4 = fig.add_subplot(144, projection='3d')
ax4.scatter(x[0:NPoly], y[0:NPoly], z[0:NPoly], c='r', marker='o')
ax4.scatter(x[NPoly:NPoly+NSol], y[NPoly:NPoly+NSol], z[NPoly:NPoly+NSol], c='b', marker='o')
#ax.set_xlabel('X Label')
plt.tight_layout()
plt.show()


#############################################################################
f = open('/Users/kamal/files/pythoncodegenerated/Polymercoord.input','w')

f.write('LAMMPS data file' + '\n')
f.write( '\n')

f.write(' '+ str(NPoly+NSol) + ' atoms'+ '\n')
f.write(' '+ str(NPoly-1) + ' bonds'+ '\n')
f.write(' '+ str(NPoly-2) + ' angles'+ '\n')
f.write(' '+ str(NPoly-3) + ' dihedrals'+ '\n')
f.write(' '+ str(NPoly-4) + ' impropers'+ '\n')
f.write( '\n')

f.write(' '+ str(atomtype)  +' atom types'+ '\n')
f.write(' '+ str(bondtype)  +' bond types'+ '\n')
f.write(' '+ str(angletype)  +' angle types'+ '\n')
f.write(' '+ str(dihedraltype)  +' dihedral types'+ '\n')
f.write(' '+ str(impropertype)  +' improper types'+ '\n')
f.write( '\n')

f.write(' '+str(-R)+' ' + str(R) +' '  'xlo xhi'+ '\n')
f.write(' '+str(-R)+' ' + str(R) +' '  'ylo yhi'+ '\n')
f.write(' '+str(-R)+' ' + str(R) +' '  'zlo zhi'+ '\n')
f.write( '\n')

f.write(' Masses'+ '\n')
f.write( '\n')

f.write(' '+ str(masstype[0])+ ' ' + str(mass[0]) +'\n')
f.write(' '+ str(masstype[1])+ ' ' + str(mass[1]) +'\n')
f.write( '\n')

f.write(' Atoms'+ '\n')
f.write( '\n')

for i in range(NPoly):
    f.write(str(i+1)+' '+'1'+' '+'1'+' '+' ' + str(qm) +' '+' '+str(x[i])+' '+str(y[i])+' '+str(z[i])+' '+'\n')
    
for i in range(NPoly, NSol+NPoly):
    f.write(str(i+1)+' '+'2'+' '+'2'+' '+' ' + str(qs) +' '+' '+str(x[i])+' '+str(y[i])+' '+str(z[i])+' '+'\n')
    
f.write( '\n')
f.write(' Bonds'+ '\n')
f.write( '\n')
    
for i in range(NPoly-1):
    f.write(str(i+1)+' '+ '1' + ' '+ str(i+1)+ ' ' + str(i+2) + '\n')

if angletype != 0:
    f.write( '\n')
    f.write(' Angles'+ '\n')
    f.write( '\n')
    for i in range(NPoly-2):
        f.write(str(i+1)+' '+ '1' + ' '+ str(i+1)+ ' ' + str(i+2) + ' ' + str(i+3)+ '\n')
else:
    pass
        
f.close()





