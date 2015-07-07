from spheretools import *

# parser=argparse.ArgumentParser(description="3-sphere plotting")
# parser.add_argument("input", type=str, help="Input file")
# parser.add_argument("--resolution", type=float, help="sphere resolution")
# parser.add_argument("--save",  type =str, help="save figure with given name")
# parser.add_argument("--R", type=float, help="radius of the sphere. By default, 0.5*(1.+sqrt(5)) is taken.")
# parser.set_defaults(resolution=13,save='none', R=0.5*(1.+pl.sqrt(5)))
# args=parser.parse_args()
# inputfile=args.input
# Res=args.resolution
# R=args.R
# save=args.save

R=0.5*(1.+pl.sqrt(5))
rcut=1.4
import sys
polar=pl.load("3Sphere-15-07-02-11h55m47s/snap%04d.npy"%int(sys.argv[1]))
# polar=pl.load("3Sphere-15-07-03-10h52m55s/snap%04d.npy"%int(sys.argv[1]))
cartesian=convert(polar, R)

dr_data=arcdistance(pdist(cartesian, 'euclidean').flatten(),R)

adjacency=squareform(dr_data)
pl.fill_diagonal(adjacency, pl.inf)

test=adjacency<rcut

import os
os.system("rm bonds.txt")
with open("bonds.txt",'wa') as bondfile:
    for p in range(adjacency.shape[0]):
        pl.savetxt(bondfile, pl.nonzero(test[p,:]), fmt="%d")

# from scipy.spatial import Voronoi

# grep 'b)' logs | sed 's/[^0-9]//g' > f.txt