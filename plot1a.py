#!/usr/bin/env python3
from numpy import *
import pylab as plt
from scipy import genfromtxt
from os.path import exists
i=1
insidearr=genfromtxt("data/inside-1.dat".format(i)).T
outsidearr=genfromtxt("data/outside-1.dat".format(i)).T
N=len(insidearr)-1
x=arange(N+1)/N
z=arange(N+1)/N
arr=hstack((insidearr,outsidearr))
xout=1/(1-x)
xall=hstack((x,xout))
zall=1/(1-z)-1
nx=len(xall[xall<1e50])
nz=len(zall[zall<1e50])
cs3=plt.contour(xall[:nx],zall[:nz],arr[:nz,:nx],(1-cos(arange(0,pi/2+0.01,pi/30)))*1.74)
plt.clabel(cs3)
plt.savefig("data/img2a.png".format(i),dpi=300)
plt.clf()
