# -*- coding: utf-8 -*-

from pkg_resources import require
require('cothread')

import sys, os
#import numpy as np
#import scipy as sp
os.environ["EPICS_CA_MAX_ARRAY_BYTES"] = '40000000'
from cothread import WaitForQuit
from cothread.catools import caget, caput, camonitor, FORMAT_TIME
from scipy.optimize import leastsq
from numpy import exp, arange, sqrt, asarray

tsBuffer = [0]*2

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
    maxV=data.max()
    minV=data.min()
    #make the array as list
    datalist=list(data)
    peak=datalist.index(maxV)
    mean=sum(x*data)/sum(data)
    sigma=sqrt(abs(sum((x-mean)**2*(data-minV))/sum(data-minV)))
#    sigma=sqrt(abs(sum((x-mean)**2*(data))/sum(data)))
#    offset=data[10]
    offset=minV
    p=[0]*5
    p[0]=maxV
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
    tsBuffer.append(value.timestamp)
    #print(tsBuffer)
    if (tsBuffer[-1] - tsBuffer[-2])>0.02:
        updateRate = 1.0/(tsBuffer[-1] - tsBuffer[-2])
        #print("system update rate: %.3f"%updateRate)
        try:
            caput('%sSysUpdateRate-I'%(cam), updateRate)
        except:
            print('can not caput SysUpdateRate')
            return
    tsBuffer.pop(0)
    #print(tsBuffer)
    
    recname=value.name
    #dire (direction/axis) is X or Y
    dire=recname[-5]
    try:
        size=caget('%sROI1:Size%s_RBV'%(cam,dire))
    except:
        print("can't get size of ROI")
        return
    #data=+value[1:size-1]
    #data: array/waveform data; image profile/intensity
    wf = +value[0:size]
    initp = initguess(wf)
    #print(initp)
    bestp = gaussfit(initp,wf)
    fittedwf = peval(bestp,arange(wf.size))
    #fittederr = sum((fittedwf-wf)**2)
    fittederr = sum(((fittedwf-wf)/bestp[0])**2)/size
    #print('%s fitted error: %f'%(dire, fittederr))
    try:
        caput("%s%s-Gauss:Max-I"%(cam,dire), bestp[0])
        caput("%s%s-Gauss:Mean-I"%(cam, dire), bestp[1])
        caput("%s%s-Gauss:Sigma-I"%(cam, dire), abs(bestp[2]))
        caput("%s%s-Gauss:Offset-I"%(cam, dire), bestp[3])
        caput("%s%s-Gauss:FittedErr-I"%(cam, dire), fittederr)
        caput("%sStats1:Peak%s_RBV"%(cam, dire), initp[4])
        #print('initial %s peak index: %d'%(dire,initp[4]))
        #print('initial/fitted %s offsets: %d / %d'%(dire,initp[3],bestp[3]))
        caput("%s%s-Gauss:Data-I"%(cam,dire), fittedwf)
    except:
        print("can't caput the results out")
        return

def main():
    camonitor("%sStats1:ProfileAverageX_RBV"%(cam), callback, format=FORMAT_TIME)
    camonitor("%sStats1:ProfileAverageY_RBV"%(cam), callback, format=FORMAT_TIME)
    WaitForQuit()

if __name__ == "__main__":
    main()

