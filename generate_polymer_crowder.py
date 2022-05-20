#Author: Kamal Tripathi, IMSc.
#Last Modified 10 March 2022
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import subprocess
import sys

qm      = float(sys.argv[1])
qs      = float(sys.argv[2])
boxsize = float(sys.argv[3])

NPoly = int(400)
NSol  = int(40)
R     = 10.0000

b     = 1.122
Ofile = 'data.input'


if (NPoly ==0 or NSol==0):
    atomtype = 1
    bondtype = 0
else:
    atomtype = 2
    bondtype = 1
    
angletype    = 1
dihedraltype = 0
impropertype = 0

masstype = []
mass     = []
for i in range(atomtype):
    masstype.append(int(i+1))
    mass.append(float(1))
    

x, y, z = [], [], []
x.append(0.0); y.append(0.0); z.append(0.0)


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

if NPoly:
    bondlen = np.zeros(NPoly-1)
    for i in range(NPoly-1):
        bondlen[i] = np.sqrt((x[i+1] - x[i])**2 + (y[i+1]-y[i])**2 + (z[i+1] - z[i])**2)
    
r = np.zeros(NPoly+NSol)
for i in range(NPoly+NSol):
    r[i] = np.sqrt((x[i])**2 + (y[i])**2 + (z[i])**2)
    
    
fig = plt.figure(figsize=(15,3))
ax1 = fig.add_subplot(141)
plt.plot(x,'o-')
plt.plot(y,'o-')
plt.plot(z,'o-')
ax1.set_xlabel('particle index')  
ax1.set_ylabel('X, Y, Z')  

ax2 = fig.add_subplot(142)
if NPoly:
    plt.plot(bondlen,'o-')     
ax2.set_xlabel('particle index')    
ax2.set_ylabel('bondlength')
        
ax3 = fig.add_subplot(143)
plt.plot(r,'o-')
ax3.set_xlabel('particle index')  
ax3.set_ylabel('$sqrt(X^2+Y^2+Z^2)$')   

ax4 = fig.add_subplot(144, projection='3d')
ax4.scatter(x[0:NPoly], y[0:NPoly], z[0:NPoly], c='r', marker='o')
ax4.scatter(x[NPoly:NPoly+NSol], y[NPoly:NPoly+NSol], z[NPoly:NPoly+NSol], c='b', marker='o')
ax4.set_xlabel('X'); ax4.set_ylabel('Y'); ax4.set_zlabel('Z')
plt.tight_layout()
#plt.show()

#############################################################################
f = open(Ofile,'w')
f.write(f'LAMMPS data file \n\n')
f.write(f" {NPoly+NSol} atoms \n")
if bondtype:
    f.write(f' {NPoly-1} bonds \n')
else:
    f.write(f' 0 bonds \n')
    
if angletype:
    f.write(f' {NPoly-2} angles \n')
else:
    f.write(f' 0 angles \n')
      
if dihedraltype:
    f.write(f' {NPoly-3} dihedrals \n')
else:
    f.write(f' 0 dihedrals \n')
    
if impropertype:
    f.write(f' {NPoly-4} impropers \n\n')
else:
    f.write(f' 0 impropers \n\n')

f.write(f' {atomtype} atom types \n')
f.write(f' {bondtype} bond types \n')
f.write(f' {angletype} angle types \n')
f.write(f' {dihedraltype} dihedral types \n')
f.write(f' {impropertype} improper types \n\n')

f.write(f' {-boxsize} {boxsize} xlo xhi \n')
f.write(f' {-boxsize} {boxsize} ylo yhi \n')
f.write(f' {-boxsize} {boxsize} zlo zhi \n\n')

f.write(f' Masses \n\n')
f.write(f' {masstype[0]} {mass[0]} \n')
if NPoly and NSol:
    f.write(f' {masstype[1]} {mass[1]} \n')


f.write(f'\n Atoms \n\n')
for i in range(NPoly):
    f.write(f'{i+1} 1 1  {qm}  {x[i]} {y[i]} {z[i]} \n')
for i in range(NPoly, NSol+NPoly):
    f.write(f'{i+1} 2 2  {qm}  {x[i]} {y[i]} {z[i]} \n')
    

if NPoly:
    f.write(f'\n Bonds \n\n')
    for i in range(NPoly-1):
        f.write(f'{i+1} 1 {i+1} {i+2} \n')

if angletype:
    f.write('\n Angles \n\n')
    for i in range(NPoly-2):
        f.write(f'{i+1} 1 {i+1} {i+2} {i+3} \n')
        
f.close()

topofile = './../../../topology-argv.tcl'
p1 = subprocess.run(['vmd', "-dispdev", "text" , "-e", topofile, '-args', 'data'])
