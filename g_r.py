import argparse

from spheretools import *

parser=argparse.ArgumentParser(description="3-sphere g(r)")
parser.add_argument("folder", type=str, help="Input Folder")
parser.add_argument("start", type=int, help="first frame")
parser.add_argument("end", type=int, help="last frame")
parser.add_argument("step", type=int, help="frame step")
parser.add_argument("--R", type=float, help="radius of the sphere. By default, 0.5*(1.+sqrt(5)) is taken.")
parser.add_argument("--nbins", type=int, help="number of bins. Default 50" )
parser.add_argument("--save",  dest='save',action='store_true', help="save result")
parser.add_argument("--plot",  dest='plot',action='store_true', help="plot result")
parser.set_defaults(nbins=50, R=0.5*(1.+pl.sqrt(5)), plot=False,save=False)
args=parser.parse_args()


import time
from scipy import stats, integrate

import os

R= args.R
folder=args.folder
start=args.start
end=args.end
step=args.step
nbins=args.nbins
plot=args.plot
save=args.save

bins=pl.linspace(0,pl.pi*R,nbins)
bincenters = 0.5*(bins[1:]+bins[:-1])


g_r=pl.zeros(nbins-1)
count=0
for i in xrange(start, end,step ):
    print "Frame",i,"sample",count
    polar=pl.load(folder+"snap%04d.npy"%i)
    cartesian=convert(polar, R)
    #calculate distances 
    dr_data=arcdistance(pdist(cartesian, 'euclidean').flatten(),R)
    # compute the ideal gas
    ideal=pl.uniform(0,2*pl.pi, size=(polar.shape[0],3))
    ideal_cartesian=convert(ideal, R)
    #calculate distances 
    dr_Ideal=arcdistance(pdist(ideal_cartesian, 'euclidean').flatten(),R)

    H,bins,binnumbers=stats.binned_statistic(dr_data, dr_data, statistic='count', bins=bins)
    H_Ideal,binsID,binnumbersID=stats.binned_statistic(dr_Ideal, dr_Ideal, statistic='count', bins=bins)

    g_r+=H/H_Ideal
    count+=1


g_r/=count
if(save):
    pl.savetxt(folder+"radial_distr.txt", zip(bincenters,g_r))
if(plot):
    pl.clf()
    pl.figure(figsize=(4,4))
    pl.subplot(2,1,1)
    pl.plot(bincenters,g_r)
    pl.ylabel("g(r)")
    pl.subplot(2,1,2)
    pl.plot(bincenters[:-1],integrate.cumtrapz(g_r,x=bincenters))
    pl.ylabel("IRDF(r)")
    pl.xlabel("r")
    pl.show()
