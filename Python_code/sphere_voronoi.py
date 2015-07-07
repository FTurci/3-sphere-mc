import pylab as pl

from scipy.spatial import ConvexHull

# On the 2-Sphere, the Voronoi tesselation is equivalent to the convex hull projected on the sphere
# (Sugihara, Journal for Geometry and Graphics Volume 6 (2002), No. 1, 69-81.)
# I assume that the same is true in 4D.... [This has to be checked!]

R= 1.6180339887498949 #magic number by straley for 120 particles

import sys
polar=pl.load(sys.argv[1])
from spheretools import *
cartesian=convert(polar, R)

CHull=ConvexHull(cartesian)

with open("bonds.txt",'w') as fw:
    for p in range(cartesian.shape[0]):
        # print p
        which_simplex,position=pl.where(CHull.simplices==p)
        # print which_simplex
        all_neighs=pl.unique(CHull.simplices[which_simplex].flatten())
        # print "all_neighs",all_neighs
        index_of_p=pl.where(all_neighs==p)
        # print "p is at",index_of_p
        neighs=pl.delete(all_neighs,index_of_p)
        # print"neighs after ", neighs
        fw.write(str(len(neighs))+" "+" ".join(map(str, neighs))+"\n" )