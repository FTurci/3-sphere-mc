import pylab as pl
from scipy.spatial.distance import pdist,squareform

def hyperarea(R):
    return 2*pl.pi**2*R**3
def radius(rho,N):
    V=N/rho
    return (V/(2*pl.pi**2))**(1./3.)

def convert(polar,R):
    cartesian=pl.zeros((polar.shape[0],4))
    cartesian[:,0]=R * pl.sin(polar[:,1]) * pl.cos(polar[:,0]);
    cartesian[:,1]=R * pl.sin(polar[:,1]) * pl.sin(polar[:,0]) * pl.cos(polar[:,2]);
    cartesian[:,2]=R * pl.cos(polar[:,1]);
    cartesian[:,3]=R * pl.sin(polar[:,1]) * pl.sin(polar[:,0]) * pl.sin(polar[:,2]);
    return cartesian

def arcdistance(d,R):
    return R * pl.arccos(1 - 0.5*d*d / (R*R))

def energy(cartesian,R):
    dists=pdist(cartesian, metric='euclidean')
    # dists=arcdistance(dists, R)
    # a la Straley
    # return (4.*dists**-12).sum()
    return (4.*(dists**-12-dists**-6)).sum()