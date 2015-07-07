from mayavi import mlab
from spheretools import *
import pylab as pl
import argparse

parser=argparse.ArgumentParser(description="3-sphere plotting")
parser.add_argument("input", type=str, help="Input file")
parser.add_argument("--resolution", type=float, help="sphere resolution")
parser.add_argument("--save",  type =str, help="save figure with given name")
parser.add_argument("--show",  type =bool, help="show the figure")
parser.add_argument("--R", type=float, help="radius of the sphere. By default, 0.5*(1.+sqrt(5)) is taken.")
parser.set_defaults(resolution=13,save='none', R=0.5*(1.+pl.sqrt(5)), show=True)
args=parser.parse_args()
inputfile=args.input
Res=args.resolution
R=args.R
save=args.save
show=args.show

polar=pl.load(inputfile)
cartesian=convert(polar, R)
x,y,z,w=pl.hsplit(cartesian,cartesian.shape[1])
mlab.clf()
mlab.points3d(x,y,z,w, resolution=Res)
if save is not 'none':
    mlab.savefig(save)
if show:
    mlab.show()
