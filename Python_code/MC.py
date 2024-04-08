import time
import argparse
from spheretools import *
import os

parser=argparse.ArgumentParser(description="3-sphere MC simulation")

parser.add_argument("T", type=float, help="temperature")
parser.add_argument("iterations", type=int, help="number of iterations")
parser.add_argument("--N", type=int, help="number of particles. Default is 120." )
parser.add_argument("--R", type=float, help="radius of the sphere. By default, 0.5*(1.+sqrt(5)) is taken.")
parser.add_argument("--step", type=float, help="Angular step" )
parser.add_argument("--snapshots", type=int, help="Spacing between snapshots" )
parser.add_argument("--output", type=str, help="Output Folder")
parser.add_argument("--input", type=str, help="Starting configuration")

parser.set_defaults(step=0.01, R=0.5*(1.+pl.sqrt(5)),N=120, input='none', output="3Sphere-"+time.strftime("%y-%m-%d-%Hh%Mm%Ss"),snapshots=10000)
args=parser.parse_args()

N=args.N
duration=args.iterations
snaps=args.snapshots
R= args.R
domega=args.step
T=args.T
folder=args.output
inputfile=args.input

os.system("mkdir "+folder)

info="""## 3-Sphere MC simulation 
* LJ
* Number of particles\t%d
* Radius\t%g
* Temperature\t%g 
* Angular step\t%g 
* Length of a step\t%g
"""%(N,R,T,domega, domega*R)

with open(folder+"/Readme.md",'w') as fw:
    fw.write(info)
print (info)

print ("Starting the simulation...")
# setup particle positions

if(inputfile != 'none'):
    polar=pl.load(inputfile)
else:
    # random start
    polar=pl.uniform(0,2*pl.pi, size=(N,3))

cartesian=convert(polar, R)

fout=open(folder+"/Arclog.txt",'w')


print ("Iteration\tEnergy")
for x in range(1,duration):
    # attemp move for all the N particles
    OldEnergy=energy(cartesian,R)
    oldpolar=pl.copy(polar)
    polar+=pl.uniform(-domega,domega ,size=(N,3))
    cartesian=convert(polar, R)
    NewEnergy=energy(cartesian,R)
    dE=NewEnergy-OldEnergy
    # print dE
    if(pl.uniform(0,1)>pl.exp(-dE/T)):
        polar=oldpolar
        cartesian=convert(polar, R)
        Energy=OldEnergy
    else:
        Energy=NewEnergy

    fout.write( "%i %g\n"%(x, Energy/N))
    # save results:
    if x%snaps==0:
        print ("%d\t%g"%(x,Energy/N))
        pl.save(folder+"/snap%04d"%x, polar)
fout.close()
    
