#!/usr/bin/env python3
from numpy import *
import pylab as plt
from scipy import genfromtxt
from os.path import exists
from sys import argv
i=int(argv[1])
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
plt.xlabel("Xin")
plt.ylabel("Zin")
plt.grid()
plt.subplot(222)
#plt.imshow(outsidearr)
cs2=plt.contour(x,z,outsidearr,arange(0,3,0.25))
plt.clabel(cs2)
plt.xlabel("Xout")
plt.ylabel("Zout")
plt.grid()
plt.subplot(223)
plt.plot(Psis,AAs)
plt.xlabel("Psi")
plt.ylabel("AA'")
plt.minorticks_on()
plt.grid(which='both')
#integrate
As=zeros(len(AAs))
As[-1]=0;
for j in range(len(AAs)-1,1,-1):
    As[j-1]=As[j]+(Psis[j-1]-Psis[j])*(AAs[j-1]+AAs[j])
As=sqrt(As)
plt.subplot(224)
plt.plot(Psis,As)
plt.xlabel("Psi")
plt.ylabel("-A")
plt.minorticks_on()
plt.grid(which='both')
plt.tight_layout()
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
plt.xlabel("x")
plt.ylabel("z")
plt.grid()
plt.savefig("data/img2-{0:02d}.png".format(i),dpi=300)
plt.clf()
