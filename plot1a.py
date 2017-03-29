#!/usr/bin/env python3
from numpy import *
import pylab as plt
from scipy import genfromtxt
arr=genfromtxt("data/all-1.dat").T
N=len(arr)-1
xall=arange(N+1)/N*10
zall=arange(N+1)/N*10
cs3=plt.contour(xall,zall,arr,(1-cos(arange(0,pi/2+0.01,pi/30)))*1.74)
plt.clabel(cs3)
plt.savefig("data/img2a.png",dpi=300)
plt.clf()
