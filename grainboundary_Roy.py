# -*- coding: utf-8 -*-
"""
Created on Mon May 21 13:27:22 2018

@author: Swagata Roy
"""


import numpy as np
import math
lattice_const=2.8553
cutoff=0.00001*lattice_const
vector1=[]
vector2=[]
cutrad=2.012999999999998
def delatom(cutrad,supercells1,ddx2,ddy2,ddz2):
    storesuper=[]
    for i in range(len(supercells1)):
        for j in range(i+1,len(supercells1)):
            distx=supercells1[i][0]-supercells1[j][0]
            disty=supercells1[i][1]-supercells1[j][1]
            distz=supercells1[i][2]-supercells1[j][2]
            if (distx)> 0.5*ddx2:
                distx=distx-ddx2
            if distx<=-0.5*ddx2:
                distx=distx+ddx2
            if (disty)> 0.5*ddy2:
                disty=disty-ddy2
            if disty<=-0.5*ddy2:
                disty=disty+ddy2
            if (distz)> 0.5*ddz2:
                distz=distz-ddz2
            if distz<=-0.5*ddz2:
                distz=distz+ddz2
            if (math.sqrt(((distx)**2) +((disty)**2) +((distz)**2)))<cutrad:
                storesuper.append([i,j])
            
    return storesuper
    
def distance(a,b):
   
    x=a[0]-b[0]
    y=a[1]-b[1]
    z=a[2]-b[2]
#    if abs(x)>0.5*n:
#        x=x-n
#    if abs(y)>0.5*n:
#        y=y-n
#    if abs(z)>0.5*n:
#        z=z-n
    dist=math.sqrt(((x)**2) +((y)**2) +((z)**2))
    return dist
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    if np.dot(axis,axis)!=0:
        axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

axis= [float(x) for x in input('rotation axisas [   ] ').split()]
theta=math.atan(math.sqrt(3)/3)
costheta=(np.dot(axis,[axis[0],axis[1],0]))/(distance(axis,[0,0,0])*distance([axis[0],axis[1],0],[0,0,0]))
theta1=math.acos(costheta)
costheta1=(np.dot([axis[0],axis[1],0],[axis[0],0,0]))/(distance([axis[0],axis[1],0],[0,0,0])*distance([axis[0],0,0],[0,0,0]))
theta2=math.acos(costheta1)
    
supercell=[10, 10, 10]
def create_bcc(supercell,lattice_const):
    vector=[]
    for i in range(0,supercell[0]):
        for j in range(0,supercell[1]):
            for k in range(0,supercell[2]):
                vector.append([i*lattice_const,j*lattice_const,k*lattice_const])
                vector.append([(i+0.5)*lattice_const,(j+0.5)*lattice_const,(k+0.5)*lattice_const])
    return vector
vector1=create_bcc(supercell,lattice_const)
vector2=create_bcc(supercell,lattice_const)
ccentre=[5*lattice_const,5*lattice_const,5*lattice_const]
for i in range(len(vector1)):
    vector1[i][0]=vector1[i][0]-ccentre[0]
    vector1[i][1]=vector1[i][1]-ccentre[1]
    vector1[i][2]=vector1[i][2]-ccentre[2]
rotate=[np.dot(vector1,rotation_matrix(np.cross(axis,[axis[0],axis[1],0]),-theta1))]
z=[np.dot([[0,0,1]],rotation_matrix(np.cross(axis,[axis[0],axis[1],0]),theta1))]
if min(z[0][0])!=0:
    z=z/min(z[0][0])
rotate=[np.dot(rotate[0],rotation_matrix(np.cross([axis[0],axis[1],0],[axis[0],0,0]),-theta2))]
y=[np.dot([[0,1,0]],rotation_matrix(np.cross([axis[0],axis[1],0],[axis[0],0,0]),theta2))]
if min(y[0][0])!=0:    
    y=y/min(y[0][0])
centre=[np.dot([lattice_const*5,lattice_const*5,lattice_const*5],rotation_matrix(np.cross(axis,[1,0,0]),-theta1))]
"""for i in range(len(rotate[0])):
    rotate[0][i][0]=rotate[0][i][0]-centre[0][0]
    rotate[0][i][1]=rotate[0][i][1]-centre[0][1]
    rotate[0][i][2]=rotate[0][i][2]-centre[0][2]"""
        
storerotatex=[]
storerotatey=[]
storerotatez=[]
for i in range(len(rotate[0])):
    if abs(rotate[0][i][1])<=cutoff and abs(rotate[0][i][2])<=cutoff:
        storerotatex.append(rotate[0][i])
    elif abs(rotate[0][i][0])<=cutoff and abs(rotate[0][i][2])<=cutoff:
        storerotatey.append(rotate[0][i])
    elif abs(rotate[0][i][0])<=cutoff and abs(rotate[0][i][1])<=cutoff:
        storerotatez.append(rotate[0][i])
storedx=[]
storedx1=[]
storedx2=[]
storedy=[]
storedy1=[]
storedy2=[]
storedz=[]
storedz1=[]
storedz2=[]
for i in range(0,len(storerotatex)):
    if storerotatex[i][0]<-cutoff:
        storedx1.append(storerotatex[i][0])
    elif storerotatex[i][0]>cutoff:
        storedx2.append(storerotatex[i][0])
dx=[]
dx.insert(0,max(storedx1))
dx.insert(1,min(storedx2))
for i in range(0,len(storerotatey)):
   if storerotatey[i][1]<=cutoff:
       m=i
for i in range(0,len(storerotatey)):
    if storerotatey[i][1]<-cutoff:
        storedy1.append(storerotatey[i][1])
    elif storerotatey[i][1]>cutoff:
        storedy2.append(storerotatey[i][1])
dy=[]
dy.insert(0,max(storedy1))
dy.insert(1,min(storedy2))
for i in range(0,len(storerotatez)):
   if storerotatez[i][2]<=cutoff:
       m=i
for i in range(0,len(storerotatez)):
    if storerotatez[i][2]<-cutoff:
        storedz1.append(storerotatez[i][2])
    elif storerotatez[i][2]>cutoff:
        storedz2.append(storerotatez[i][2])
dz=[]
dz.insert(0,max(storedz1))
dz.insert(1,min(storedz2))
unitcell=[]
for i in range(len(rotate[0])):
    if rotate[0][i][0]<=(dx[1]+cutoff)and rotate[0][i][0]>=(dx[0]-cutoff) and rotate[0][i][1]<=(dy[1]+cutoff)and rotate[0][i][1]>=(dy[0]-cutoff) and rotate[0][i][2]<=(dz[1]+cutoff)and rotate[0][i][2]>=(dz[0]-cutoff):
        unitcell.append(rotate[0][i])    
for i in range(len(unitcell)):
    unitcell[i][0]=unitcell[i][0]-dx[0]
    unitcell[i][1]=unitcell[i][1]-dy[0]
    unitcell[i][2]=unitcell[i][2]-dz[0]
#delstore=[]
#for i in range(len(unitcell)):
#    if i in delstore:
#        continue
#    for j in range(i+1,len(unitcell)):
#        if (abs(unitcell[i][0]-unitcell[j][0])<cutoff and abs(unitcell[i][1]-unitcell[j][1])<cutoff and abs(unitcell[i][2]-unitcell[j][2])<cutoff):
#            delstore.append(j)
#for k in sorted(delstore,reverse=True):
#    del unitcell[k]  
size=[5,5,5]
supercells=[]
for i in range(size[0]):
    for j in range(size[1]):
        for k in range(size[2]):
            for m in range(len(unitcell)):
                supercells.append([unitcell[m][0]+2*dx[1]*i,unitcell[m][1]+2*dy[1]*j,unitcell[m][2]+2*dz[1]*k])
#delstore=[]
#for i in range(len(supercells)):
#   if i in delstore:
#       continue
#   for j in range(i+1,len(supercells)):
#       if (abs(supercells[i][0]-supercells[j][0])<cutoff and abs(supercells[i][1]-supercells[j][1])<cutoff and abs(supercells[i][2]-supercells[j][2])<cutoff):
#           delstore.append(j)
#for k in sorted(delstore,reverse=True):
#   del supercells[k]  
y=np.dot(y,rotation_matrix(axis,theta))
z=np.dot(z,rotation_matrix(axis,theta))
axis=[1,0,0]
rotate1=[np.dot(supercells,rotation_matrix(axis,theta))]
rotate2=[np.dot(supercells,rotation_matrix(axis,-theta))]

dmx=[]
dmy=[]
dmz=[]
for i in range(len(supercells)):
    dmx.append(supercells[i][0])
    dmy.append(supercells[i][1])
    dmz.append(supercells[i][2])
ddmx=min(dmx)
ddmx1=max(dmx)
ddmy=min(dmx)
ddmy1=max(dmy)
ddmz=min(dmz)
ddmz1=max(dmz)
centre1=[np.dot([(ddmx+ddmx1)/2.0,(ddmy+ddmy1)/2.0,(ddmz+ddmz1)/2.0],rotation_matrix(axis,theta))]
centre2=[np.dot([(ddmx+ddmx1)/2.0,(ddmy+ddmy1)/2.0,(ddmz+ddmz1)/2.0],rotation_matrix(axis,-theta))]
for i in range(len(rotate1[0])):
    rotate1[0][i][0]=rotate1[0][i][0]-centre1[0][0]
    rotate1[0][i][1]=rotate1[0][i][1]-centre1[0][1]
    rotate1[0][i][2]=rotate1[0][i][2]-centre1[0][2]
for i in range(len(rotate2[0])):
    rotate2[0][i][0]=rotate2[0][i][0]-centre2[0][0]
    rotate2[0][i][1]=rotate2[0][i][1]-centre2[0][1]
    rotate2[0][i][2]=rotate2[0][i][2]-centre2[0][2]
store1=[]
store2=[]
for i in range(0,len(rotate1[0])):
    for j in range(0,len(rotate2[0])):
        if (abs(rotate1[0][i][0]-rotate2[0][j][0])<cutoff and abs(rotate1[0][i][1]-rotate2[0][j][1])<cutoff and abs(rotate1[0][i][2]-rotate2[0][j][2])<cutoff):
            store1.append(rotate1[0][i])
            store2.append(rotate2[0][j])
storedx=[]
storedx1=[]
storedx2=[]
storedy=[]
storedy1=[]
storedy2=[]
storedz=[]
storedz1=[]
storedz2=[]
for i in range(0,len(store1)):
    if abs(store1[i][1])<=cutoff and abs(store1[i][2])<=cutoff :
        storedx.append(store1[i][0])
    if  abs(store1[i][0])<=cutoff and abs(store1[i][2])<=cutoff :
        storedy.append(store1[i][1])
    if  abs(store1[i][1])<=cutoff and abs(store1[i][0])<=cutoff :
        storedz.append(store1[i][2])
for i in range(0,len(storedx)):
   if storedx[i]<=cutoff:
       m=i
for i in range(0,len(storedx)):
    if storedx[i]<-cutoff:
        storedx1.append(storedx[i])
    elif storedx[i]>cutoff:
        storedx2.append(storedx[i])
dx=[]
dx.insert(0,max(storedx1)) 
dx.insert(1,min(storedx2))
if  min(storedx2) < lattice_const-cutoff:
    dx=[]
    dx.insert(0,2*max(storedx1)) ##################
    dx.insert(1,2*min(storedx2))
for i in range(0,len(storedy)):
   if storedy[i]<=cutoff:
       m=i
for i in range(0,len(storedy)):
    if storedy[i]<-cutoff:
        storedy1.append(storedy[i])
    elif storedy[i]>cutoff:
        storedy2.append(storedy[i])
dy=[]
dy.insert(0,max(storedy1))
dy.insert(1,min(storedy2))
if  min(storedy2) < lattice_const-cutoff:
    dy=[]
    dy.insert(0,2*max(storedy1)) ##################
    dy.insert(1,2*min(storedy2))
for i in range(0,len(storedz)):
   if storedz[i]<=cutoff:
       m=i
for i in range(0,len(storedz)):
    if storedz[i]<-cutoff:
        storedz1.append(storedz[i])
    elif storedz[i]>cutoff:
        storedz2.append(storedz[i])
dz=[]
dz.insert(0,max(storedz1))
dz.insert(1,min(storedz2))
if  min(storedz2) < lattice_const-cutoff:
    dz=[]
    dz.insert(0,2*max(storedz1)) ##################
    dz.insert(1,2*min(storedz2))
unitcell=[]
supercellso1=[]
supercellso2=[]

for i in range(len(rotate1[0])):
    if rotate1[0][i][0]<=(0.5*dx[1]+cutoff)and rotate1[0][i][0]>=(0.5*dx[0]-cutoff) and rotate1[0][i][1]<=(0.5*dy[1]+cutoff)and rotate1[0][i][1]>=(0.5*dy[0]-cutoff) and rotate1[0][i][2]<=(0.5*dz[1]+cutoff)and rotate1[0][i][2]>=(0.5*dz[0]-cutoff):
        unitcell.append(rotate1[0][i])
        supercellso1.append(rotate1[0][i])
for i in range(len(rotate2[0])):
    if rotate2[0][i][0]<=(0.5*dx[1]+cutoff)and rotate2[0][i][0]>=(0.5*dx[0]-cutoff) and rotate2[0][i][1]<=(0.5*dy[1]+cutoff)and rotate2[0][i][1]>=(0.5*dy[0]-cutoff) and rotate2[0][i][2]<=(0.5*dz[1]+cutoff)and rotate2[0][i][2]>=(0.5*dz[0]-cutoff):
        unitcell.append(rotate2[0][i])
        supercellso2.append(rotate2[0][i])
"""delstoreunit=[]
for i in range(len(unitcell)):
    if unitcell[i][0]==dx[0] or unitcell[i][0]==dx[1] or unitcell[i][1]==dy[0] or unitcell[i][1]==dy[1] or unitcell[i][2]==dz[0] or unitcell[i][2]==dz[1]:
        delstoreunit.append(i)
for k in sorted(delstoreunit,reverse=True):
    del unitcell[k]"""
#unitcell.append([dx[0],dy[0],dz[0]])
#supercellso1.append([dx[0],dy[0],dz[0]])
#supercellso2.append([dx[0],dy[0],dz[0]])
for i in range(len(unitcell)):
    unitcell[i][0]=unitcell[i][0]-0.5*dx[0]
    unitcell[i][1]=unitcell[i][1]-0.5*dy[0]
    unitcell[i][2]=unitcell[i][2]-0.5*dz[0]
"""for i in range(len(supercellso1)):
    supercellso1[i][0]=supercellso1[i][0]-dx[0]
    supercellso1[i][1]=supercellso1[i][1]-dy[0]
    supercellso1[i][2]=supercellso1[i][2]-dz[0]
for i in range(len(supercellso2)):
    supercellso2[i][0]=supercellso2[i][0]-dx[0]
    supercellso2[i][1]=supercellso2[i][1]-dy[0]
    supercellso2[i][2]=supercellso2[i][2]-dy[0]"""
supercells=[]
supercellsso1=[]
supercellsso2=[]
size=[2,2,2]
for i in range(size[0]):
    for j in range(size[1]):
        for k in range(size[2]):
            for m in range(len(unitcell)):
                supercells.append([unitcell[m][0]+dx[1]*i,unitcell[m][1]+dy[1]*j,unitcell[m][2]+dz[1]*k])
            for l in range(len(supercellso1)):
                supercellsso1.append([supercellso1[l][0]+dx[1]*i,supercellso1[l][1]+dy[1]*j,supercellso1[l][2]+dz[1]*k])
            for h in range(len(supercellso2)):   
                supercellsso2.append([supercellso2[h][0]+dx[1]*i,supercellso2[h][1]+dy[1]*j,supercellso2[h][2]+dz[1]*k])
ddx=[]
ddy=[]
ddz=[]
for i in range(0,len(supercells)):
    ddx.append(supercells[i][0])
    ddy.append(supercells[i][1])
    ddz.append(supercells[i][2])
ddx1=min(ddx)
ddx2=max(ddx)
ddy1=min(ddy)
ddy2=max(ddy)
ddz1=min(ddz)
ddz2=max(ddz)     
gb=(ddy2+ddy1)/2
supercells1=[]
for i in range(len(supercells)):
    if supercells[i][1]<(gb-cutoff):
        for j in range(len(supercellsso1)):
            if supercells[i][0]==supercellsso1[j][0] and supercells[i][1]==supercellsso1[j][1] and supercells[i][2]==supercellsso1[j][2]:
                supercells1.append(supercells[i])
    elif supercells[i][1]>(gb+cutoff):
        for k in range(len(supercellsso2)):
            if supercells[i][0]==supercellsso2[k][0] and supercells[i][1]==supercellsso2[k][1] and supercells[i][2]==supercellsso2[k][2]:
                supercells1.append(supercells[i])
    else:
        supercells1.append(supercells[i])
delstore=[]
for i in range(len(supercells1)):
    if i in delstore:
        continue
    for j in range(i+1,len(supercells1)):
        if (abs(supercells1[i][0]-supercells1[j][0])<cutoff and abs(supercells1[i][1]-supercells1[j][1])<cutoff and abs(supercells1[i][2]-supercells1[j][2])<cutoff):
            delstore.append(j)
for k in sorted(delstore,reverse=True):
    del supercells1[k]     
boundarycheck=[]
for i in range(len(supercells1)):
    if abs(supercells1[i][0]-ddx2)<cutoff or abs(supercells1[i][1]-ddy2)<cutoff or abs(supercells1[i][2]-ddz2)<cutoff:# change ddz2 to 1 
        boundarycheck.append(i)
for k in sorted(boundarycheck,reverse=True):
    del supercells1[k]
#boundarycheck=[]
#for i in range(len(supercells1)):
#    if abs(supercells1[i][0]-ddx1)<cutoff or abs(supercells1[i][1]-ddy1)<cutoff or abs(supercells1[i][2]-ddz1)<cutoff:
#        boundarycheck.append(i)
#for k in sorted(boundarycheck,reverse=True):
#    del supercells1[k]
#storesuper=delatom(cutrad,supercells1,ddx2,ddy2,ddz2)
#storesupers=[]     
#storesuper1=[]
#for i in range(len(storesuper)):
#     storesuper1.append(storesuper[i][0])
#     storesuper1.append(storesuper[i][1])
#    
##for i in range(len(storesuper)):
##    if abs(supercells1[storesuper[i][0]][2]-supercells1[storesuper[i][1]][2])>0.5*ddz2:
##        storesupers.append([0.5*(supercells1[storesuper[i][0]][0]+supercells1[storesuper[i][1]][0]),0.5*(supercells1[storesuper[i][0]][1]+supercells1[storesuper[i][1]][1]),0.5*(supercells1[storesuper[i][0]][2]+supercells1[storesuper[i][1]][2]+ddz2)])
##    else:
##         storesupers.append([0.5*(supercells1[storesuper[i][0]][0]+supercells1[storesuper[i][1]][0]),0.5*(supercells1[storesuper[i][0]][1]+supercells1[storesuper[i][1]][1]),0.5*(supercells1[storesuper[i][0]][2]+supercells1[storesuper[i][1]][2])])
##    
#for k in sorted(storesuper1,reverse=True):
#   del supercells1[k]                      
#for i in range(len(storesupers)):
#    supercells1.append(storesupers[i])   
delstoreunit=[]
"""for i in range(len(supercells1)):
    if (supercells1[i][0]>ddx1-cutoff and supercells1[i][0]<ddx1+cutoff) or (supercells1[i][0]>ddx2-cutoff and supercells1[i][0]<ddx2+cutoff) or (supercells1[i][1]>ddy1-cutoff and supercells1[i][1]<ddy1+cutoff) or (supercells1[i][1]>ddy2-cutoff and supercells1[i][1]<ddy2+cutoff) or (supercells1[i][2]>ddz1-cutoff and supercells1[i][2]<ddz1+cutoff) or (supercells1[i][2]>ddz2-cutoff and supercells1[i][2]<ddz2+cutoff):
        delstoreunit.append(i)
for k in sorted(delstoreunit,reverse=True):
    del supercells1[k]"""
inputs=open('gb111del.dump','w')
'''output.write(str(len(supercells1)))
output.write('   ')
output.write('Atoms')
output.write('\n')
output.write(' 1  atom types')
output.write('\n')
output.write(str(ddx1))
output.write('   ')
output.write(str(ddx2))        
output.write('   ')
output.write('xlo        xhi')
output.write('\n')
output.write(str(ddy1))
output.write('   ')
output.write(str(ddy2))        
output.write('   ')
output.write('ylo        yhi')
output.write('\n')
output.write(str(ddz1))
output.write('   ')
output.write(str(ddz2))        
output.write('   ')
output.write('zlo        zhi')
output.write('\n')
output.write('Atoms#atoms\n')'''
inputs.write('# position_fraction -0.000676 \n')
inputs.write(' \n')
inputs.write(str(len(supercells1))+' atoms \n')
inputs.write('2	 atom types\n')
inputs.write('\n')
inputs.write(str(ddx1)+' '+str(ddx2)+' '+'  xlo xhi \n')#change accordingly
inputs.write(str(ddy1)+' '+str(ddy2)+' '+'  ylo yhi \n')
inputs.write(str(ddz1)+' '+str(ddz2)+' '+'  zlo zhi \n')
inputs.write('\n')
inputs.write(' Masses\n')
inputs.write('\n')
inputs.write(' 1	 56.000000\n')
inputs.write(' 2	  1.000000\n')
inputs.write('\n')
inputs.write(' Atoms\n')
inputs.write('\n')
for i in range(len(supercells1)):
    inputs.write(str(i+1)+'    '+' 1    '+str(supercells1[i][0])+'    '+str(supercells1[i][1])+'    '+str(supercells1[i][2])+'    '+'\n')
             
inputs.close()             

