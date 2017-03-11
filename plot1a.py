#!/usr/bin/env python3
from numpy import *
import pylab as plt
from scipy import genfromtxt
i=1
insidearr=genfromtxt("data/inside-1.dat".format(i)).T
outsidearr=genfromtxt("data/outside-1.dat".format(i)).T
N=len(insidearr)-1
arr=hstack((insidearr,outsidearr))
xall=arange(N+1)/N*10
zall=arange(N+1)/N*10
cs3=plt.contour(xall,zall,arr,(1-cos(arange(0,pi/2+0.01,pi/30)))*1.74)
plt.clabel(cs3)
plt.savefig("data/img2a.png".format(i),dpi=300)
plt.clf()
