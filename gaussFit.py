# -*- coding: utf-8 -*-

import sys, time
import numpy as np
import scipy as sp
import cothread

from cothread.catools import *
from scipy.optimize import leastsq
from numpy import exp, arange, sqrt, asarray

if len(sys.argv) != 2:
    print('must use camera prefix such as LTB-BI{VF:1BD1} for argv[1]')
    exit()

cam  = sys.argv[1]

def peval(p, x):
    return ((p[0]-p[3])*exp(-(x-p[1])**2/(2*p[2]**2))+p[3])

def residuals(p, x, y):
    return (y - peval(p,x))

def initguess(data):
    #x: axis X; array[0,data.size)
    x=arange(0,data.size,1)
    max=data.max()
    min=data.min()
    #make the array as list
    datalist=list(data)
    peak=datalist.index(max)
    mean=sum(x*data)/sum(data)
    sigma=sqrt(abs(sum((x-mean)**2*(data-min))/sum(data-min)))
#    sigma=sqrt(abs(sum((x-mean)**2*(data))/sum(data)))
#    offset=data[10]
    offset=min
    p=[0]*5
    p[0]=max
    p[1]=mean
    p[2]=sigma
    p[3]=offset
    p[4]=peak
    return asarray(p)

def gaussfit(p0, data):
    bp,_,_,sts,ret = leastsq(residuals, p0, args=(arange(data.size), data), full_output=1)
    if ret not in [1,2,3,4]:
        print('No solution %d: %s: %s'%(ret,sts, str(p0)))
    return bp

def callback(value):
    recname=value.name
    #dir (direction/axis) is X or Y
    dir=recname[-5]
    size=caget('%sROI1:Size%s_RBV'%(cam,dir))
    #data=+value[1:size-1]
    #data: array/waveform data; image profile/intensity
    wf = +value[0:size]
    initp = initguess(wf)
    #print(initp)
    bestp = gaussfit(initp,wf)
    fittedwf = peval(bestp,arange(wf.size))
    #fittederr = sum((fittedwf-wf)**2)
    fittederr = sum(((fittedwf-wf)/bestp[0])**2)/size
    #print('%s fitted error: %f'%(dir, fittederr))
    caput("%s%s-Gauss:Max-I"%(cam,dir), bestp[0])
    caput("%s%s-Gauss:Mean-I"%(cam, dir), bestp[1])
    caput("%s%s-Gauss:Sigma-I"%(cam, dir), abs(bestp[2]))
    caput("%s%s-Gauss:Offset-I"%(cam, dir), bestp[3])
    caput("%s%s-Gauss:FittedErr-I"%(cam, dir), fittederr)
    caput("%sStats1:Peak%s_RBV"%(cam, dir), initp[4])
    #print('initial %s peak index: %d'%(dir,initp[4]))
    #print('initial/fitted %s offsets: %d / %d'%(dir,initp[3],bestp[3]))
    caput("%s%s-Gauss:Data-I"%(cam,dir), fittedwf)

def main():
    camonitor("%sStats1:ProfileAverageX_RBV"%(cam), callback, format=FORMAT_RAW)
    camonitor("%sStats1:ProfileAverageY_RBV"%(cam), callback, format=FORMAT_RAW)
    cothread.WaitForQuit()

if __name__ == "__main__":
    main()

