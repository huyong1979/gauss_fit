# -*- coding: utf-8 -*-

import sys, time
import numpy as np
import scipy as sp
#import pylab
import cothread

from cothread.catools import *
from scipy.optimize import leastsq
from numpy import exp, arange, sqrt, asarray

#print('argc: %d'%(len(sys.argv)))
if len(sys.argv) != 2:
    print('must use camera prefix such as LTB-BI{VF:1BD1} for argv[1]')
    exit()

cam  = sys.argv[1]
#print('argv[0]: %s; argv[1]: %s'%(sys.argv[0],cam))
#A={'CAM':cam}
#B={'DIR':X, Y}

def peval(p, x):
    return ((p[0]-p[3])*exp(-(x-p[1])**2/(2*p[2]**2))+p[3])

def residuals(p, x, y):
    return (y - peval(p,x))

def initguess(data):
    x=arange(data.size)
    #print('Axis size: %d'%(data.size))
    #print('intensity array data: ', data)
    max=data.max()
    #make the array as list
    datalist=list(data)
    #print('convert data to list: ', datalist)
    peak=datalist.index(max)
    mean=sum(x*data)/sum(data)
    sigma=sqrt(abs(sum((x-mean)**2*data)/sum(data)))
    offset=data[10]
    p=[0]*5
    p[0]=max
    p[1]=mean
    p[2]=sigma
    p[3]=offset
    p[4]=peak
    return asarray(p)

def gaussfit(p0, data):
    pf,_,_,sts,ret = leastsq(residuals, p0, args=(arange(data.size), data), full_output=1)
    if ret not in [1,2,3,4]:
        print('No solution %d: %s'%(ret,sts))
    return pf

def callback(value):
    #print('record name is: %s'%(value.name))
    recname=value.name
    #print('record name: %s, size: %d, dir: %s'%(recname,len(recname),recname[-5]))
    #dir (direction/axis) is X or Y
    dir=recname[-5]
    #size=caget("%(CAM)sROI1:SizeX_RBV"%A)
    size=caget('%sROI1:Size%s_RBV'%(cam,dir))
    print("Axis %s size: %d"%(dir,size))
    #data=+value[1:size-1]
    data=+value[0:size]
    initial = initguess(data)
    final = gaussfit(initial,data)
    fitdata = peval(final,arange(data.size))
    caput("%s%s-Gauss:Max-I"%(cam,dir), final[0])
    caput("%s%s-Gauss:Mean-I"%(cam, dir), final[1])
    caput("%s%s-Gauss:Sigma-I"%(cam, dir), abs(final[2]))
    caput("%s%s-Gauss:Offset-I"%(cam, dir), final[3])
    caput("%sStats1:Peak%s_RBV"%(cam, dir), initial[4])
    #print('initial %s peak index: %d'%(dir,initial[4]))
    caput("%s%s-Gauss:Data-I"%(cam,dir), fitdata)

#def callbacky(value):
    #size=caget("%(CAM)sROI1:SizeY_RBV"%A)
    #data=+value[1:size-1]
    #initial = initguess(data)
    #final = gaussfit(initial,data)
    #fitdata = peval(final,arange(data.size))
    #caput("%(CAM)sY-Gauss:Max-I"%A, final[0])
    #caput("%(CAM)sY-Gauss:Mean-I"%A, final[1])
    #caput("%(CAM)sY-Gauss:Sigma-I"%A, abs(final[2]))
    #caput("%(CAM)sY-Gauss:Offset-I"%A, final[3])
    #caput("%(CAM)sStats1:PeakY_RBV"%A, initial[4])
    #print('initialY', initial[4])
    #caput("%(CAM)sY-Gauss:Data-I"%A, fitdata)

def main():
    camonitor("%sStats1:ProfileAverageX_RBV"%(cam), callback, format=FORMAT_RAW)
    camonitor("%sStats1:ProfileAverageY_RBV"%(cam), callback, format=FORMAT_RAW)
    cothread.WaitForQuit()

if __name__ == "__main__":
    main()

