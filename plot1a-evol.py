#!/usr/bin/env python3
from numpy import *
import pylab as plt
from scipy import genfromtxt
from os.path import exists
from os import popen
i=0
s=""
while exists("data/all-1-{0:d}.dat".format(i)):
    print(i)
    arr=genfromtxt("data/all-1-{0:d}.dat".format(i)).T
    N=len(arr)-1
    xall=arange(N+1)/N*10
    zall=arange(N+1)/N*10
    cs3=plt.contour(xall,zall,arr,(1-cos(arange(0,pi/2+0.01,pi/30)))*1.74)
    plt.clabel(cs3)
    plt.title("Step {0:d}".format(i))
    s+="data/img2a-{0:d}.png ".format(i)
    plt.savefig("data/img2a-{0:d}.png".format(i),dpi=300)
    plt.clf()
    i+=1
print("make gif")
popen("convert -delay 200 "+s+"data/evol.gif").read()
print("finished")
