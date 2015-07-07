R= 1.83586
import pylab as pl

def convert(polar,R):
    cartesian=pl.zeros((polar.shape[0],4))
    cartesian[:,0]=R * pl.sin(polar[:,1]) * pl.cos(polar[:,0]);
    cartesian[:,1]=R * pl.sin(polar[:,1]) * pl.sin(polar[:,0]) * pl.cos(polar[:,2]);
    cartesian[:,2]=R * pl.cos(polar[:,1]);
    cartesian[:,3]=R * pl.sin(polar[:,1]) * pl.sin(polar[:,0]) * pl.sin(polar[:,2]);
    return cartesian


polar=pl.loadtxt("start.txt")

cartesian=convert(polar, R)

def LJ(r):
    return 4*(r**-12-r**-6);

energies=[]
for i in range(0,120-1):
    for j in range(i+1,120):
        dx=cartesian[i,0]-cartesian[j,0]
        dy=cartesian[i,1]-cartesian[j,1]
        dz=cartesian[i,2]-cartesian[j,2]
        dw=cartesian[i,3]-cartesian[j,3]
        dr=pl.sqrt(dx*dx+dy*dy+dz*dz+dw*dw)
        if(dr<2.5):
            energies.append(LJ(dr))
        else: 
            energies.append(0)

print pl.sum(energies)/120.
# +0.016316891136000006