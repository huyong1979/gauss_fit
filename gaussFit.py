# -*- coding: utf-8 -*-

import sys, time
import numpy as np
import scipy as sp
import pylab
import cothread

from cothread.catools import *
from scipy.optimize import leastsq

cam  = sys.argv[1]
A={'CAM':cam}

def gauss_fit(p, x):
    return ((p[0]-p[3])*np.exp(-(x-p[1])**2/(2*p[2]**2))+p[3])

def e_gauss_fit(p, x, y):
    return (y - gauss_fit(p,x))

def initguess(data):
    X=np.arange(data.size)
    max=data.max()
    mean=sum(X*data)/sum(data)
    sigma=np.sqrt(abs(sum((X-mean)**2*data)/sum(data)))
    offset=data[10]
    P0=[0]*4
    P0[0]=max
    P0[1]=mean
    P0[2]=sigma
    P0[3]=offset
    return np.asarray(P0)

def gaussfit(initguess, data):
    Pf,_,_,sts,ret = leastsq(e_gauss_fit, initguess, args=(np.arange(data.size), data), full_output=1)
    if ret not in [1,2,3,4]:
        print('No solution %d: %s'%(ret,sts))
    return Pf

def callbackX(value):
    size=caget("%(CAM)sROI1:SizeX_RBV"%A)
    data=+value[1:size-1]
    initial = initguess(data)
    final = gaussfit(initial,data)
    fitdata = gauss_fit(final,np.arange(data.size))
    caput("%(CAM)sX-Gauss:Max-I"%A, final[0])
    caput("%(CAM)sX-Gauss:Mean-I"%A, final[1])
    caput("%(CAM)sX-Gauss:Sigma-I"%A, abs(final[2]))
    caput("%(CAM)sX-Gauss:Offset-I"%A, final[3])
    caput("%(CAM)sX-Gauss:Data-I"%A, fitdata)

def callbackY(value):
    size=caget("%(CAM)sROI1:SizeY_RBV"%A)
    data=+value[1:size-1]
    initial = initguess(data)
    final = gaussfit(initial,data)
    fitdata = gauss_fit(final,np.arange(data.size))
    caput("%(CAM)sY-Gauss:Max-I"%A, final[0])
    caput("%(CAM)sY-Gauss:Mean-I"%A, final[1])
    caput("%(CAM)sY-Gauss:Sigma-I"%A, abs(final[2]))
    caput("%(CAM)sY-Gauss:Offset-I"%A, final[3])
    caput("%(CAM)sY-Gauss:Data-I"%A, fitdata)

def main():
    camonitor("%(CAM)sStats1:ProfileAverageX_RBV"%A, callbackX, format=FORMAT_RAW)
    camonitor("%(CAM)sStats1:ProfileAverageY_RBV"%A, callbackY, format=FORMAT_RAW)
    cothread.WaitForQuit()

if __name__ == "__main__":
    main()

