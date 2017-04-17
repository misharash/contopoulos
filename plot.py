#!/usr/bin/env python3
from numpy import *
import pylab as plt
from scipy import genfromtxt
from os.path import exists
i=1
while exists("data/inside-{0:d}.dat".format(i)):
    insidearr=genfromtxt("data/inside-{0:d}.dat".format(i)).T
    outsidearr=genfromtxt("data/outside-{0:d}.dat".format(i)).T
    Psis,AAs=genfromtxt("data/table-{0:d}.dat".format(i)).T
    N=len(insidearr)-1
    x=arange(N+1)/N
    z=arange(N+1)/N
    plt.subplot(221)
    #plt.imshow(insidearr)
    cs=plt.contour(x,z,insidearr,arange(0,3,0.25))
    plt.clabel(cs)
    plt.subplot(222)
    #plt.imshow(outsidearr)
    cs2=plt.contour(x,z,outsidearr,arange(0,3,0.25))
    plt.clabel(cs2)
    plt.subplot(212)
    plt.plot(Psis,AAs)
    plt.savefig("data/img-{0:02d}.png".format(i),dpi=300)
    plt.clf()
    arr=hstack((insidearr,outsidearr))
    xout=1/(1-x)
    xall=hstack((x,xout))
    zall=1/(1-z)-1
    nx=len(xall[xall<3])
    nz=len(zall[zall<3])
    cs3=plt.contour(xall[:nx],zall[:nz],arr[:nz,:nx],arange(0,3,0.25))
    plt.clabel(cs3)
    plt.savefig("data/img2-{0:02d}.png".format(i),dpi=300)
    plt.clf()
    print(i)
    i+=1
