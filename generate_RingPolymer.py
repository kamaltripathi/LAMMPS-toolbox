#This code generates the lammps data file of multiple ring polymers which are stacked along z-axis.
#Author: Kamal Tripathi, IMSc, Chennai.
#Date: 18 Aug 2020.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
get_ipython().run_line_magic('matplotlib', 'inline')

Nr    = 20                                 #Number of ring polymers
d     = 3.0                                #distance between two ring polymers in z-plane
Nm    = 100                                #degree of polymerization
theta = 2*np.pi/Nm                         #angle made by two successive beads at the center of the polygon
b     = 1.122                              #bond length
R     = b/(2*np.sin(theta/2.0))            #radius of enclosing circle

zpos  = np.arange(-(Nr-1)*d/2.0, (Nr-1)*d/2.0 +1, d)   #z-coordinates of rings

x = np.zeros(Nm*Nr)
y = np.zeros(Nm*Nr)
z = np.zeros(Nm*Nr)

for j in range(0,Nr):
    for i in range(Nm):
        x[Nm*j+i] = R*np.cos(i*theta)
        y[Nm*j+i] = R*np.sin(i*theta)
        z[Nm*j+i] = zpos[j]

ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()

#See explanation of data_format for lammps at: https://lammps.sandia.gov/doc/2001/data_format.html
qm    = 0.00                    #charge on monomer
atomtype     = 1                #types of atom
bondtype     = 1                #types of bonds
angletype    = 1                #types of angles
dihedraltype = 0                #types of dihedral
impropertype = 0                

masstype =[1]
mass = [1.0000]

f = open('./data.input','w')

f.write('LAMMPS data file' + '\n')
f.write( '\n')

f.write(' '+ str(Nr*Nm) + ' atoms'+ '\n')
f.write(' '+ str(Nr*Nm) + ' bonds'+ '\n')
f.write(' '+ str(Nr*Nm) + ' angles'+ '\n')
f.write(' '+ str(0) + ' dihedrals'+ '\n')
f.write(' '+ str(0) + ' impropers'+ '\n')
f.write( '\n')

f.write(' '+ str(atomtype)  +' atom types'+ '\n')
f.write(' '+ str(bondtype)  +' bond types'+ '\n')
f.write(' '+ str(angletype)  +' angle types'+ '\n')
f.write(' '+ str(dihedraltype)  +' dihedral types'+ '\n')
f.write(' '+ str(impropertype)  +' improper types'+ '\n')
f.write( '\n')

f.write(' '+str(-R-1)+' ' + str(R+1) +' '  'xlo xhi'+ '\n')
f.write(' '+str(-R-1)+' ' + str(R+1) +' '  'ylo yhi'+ '\n')
f.write(' '+str(zpos[0]-1)+' ' + str(zpos[Nr-1]+1) +' '  'zlo zhi'+ '\n')
f.write( '\n')

f.write(' Masses'+ '\n')
f.write( '\n')

f.write(' '+ str(masstype[0])+ ' ' + str(mass[0]) +'\n')
f.write( '\n')

f.write(' Atoms'+ '\n')
f.write( '\n')

for i in range(Nr*Nm):
    f.write(str(i+1)+' '+str(i//Nm +1)+' '+'1'+' '+' ' + str(qm) +' '+' '+str(x[i])+' '+str(y[i])+' '+str(z[i])+' '+'\n')
    
f.write( '\n')
f.write(' Bonds'+ '\n')
f.write( '\n')
    
for i in range(Nr*Nm):
    if ((i+1)%Nm!=0):
        f.write(str(i+1)+' '+ '1' + ' '+ str(i+1)+ ' ' + str(i+2) + '\n')
    else:
        f.write(str(i+1)+' '+ '1' + ' '+ str(i+1)+ ' ' + str(i+2-Nm) + '\n')
    
if angletype != 0:
    f.write( '\n')
    f.write(' Angles'+ '\n')
    f.write( '\n')
    for i in range(Nr*Nm):
        if (0<(i+1)%Nm<(Nm-1)):
            f.write(str(i+1)+' '+ '1' + ' '+ str(i+1)+ ' ' + str(i+2) + ' ' + str(i+3)+ '\n')
        elif ((i+1)%Nm==(Nm-1)):
            f.write(str(i+1)+' '+ '1' + ' '+ str(i+1)+ ' ' + str(i+2) + ' ' + str((i+3)-Nm) + '\n')
        elif ((i+1)%Nm==0):
            f.write(str(i+1)+' '+ '1' + ' '+ str(i+1)+ ' ' + str((i+2)-Nm) + ' ' + str((i+3)-Nm) + '\n')
else:
    pass
f.close()
