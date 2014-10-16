# -*- coding: utf-8 -*-

from pkg_resources import require
require('cothread')

import sys, os
import traceback
os.environ["EPICS_CA_MAX_ARRAY_BYTES"] = '40000000'
from cothread import WaitForQuit
from cothread.catools import caget, caput, camonitor, FORMAT_TIME
from scipy.optimize import leastsq
from numpy import exp, pi, arange, sqrt, asarray, ravel, indices, cos, sin

tsBuffer = [0]*2

if len(sys.argv) != 2:
    print('must use camera prefix such as LTB-BI{VF:1BD1} for argv[1]')
    exit()

cam  = sys.argv[1]
fitResultPVs = ["%sX-Gauss:Mean-I"%(cam),"%sY-Gauss:Mean-I"%(cam),
                "%sX-Gauss:Sigma-I"%(cam),"%sY-Gauss:Sigma-I"%(cam)]
print fitResultPVs

def gaussian(height, off, center_x, center_y, width_x, width_y, rota):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    rota = float(rota)
    return lambda x,y: (height-off)*exp(-(((center_x * cos(rota) - center_y * sin(rota)-x * cos(rota) + y * sin(rota))/width_x)**2+((center_x * sin(rota) + center_y * cos(rota)-x * sin(rota) - y * cos(rota))/width_y)**2)/2) + off

def moments(data):
    """Returns (height, offset, centroidX, contrenoidY, sigmaX, sigmaY, rota)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    height = data.max()
    offset = data.min()
    (centroidX,centroidY,sigmaX,sigmaY) = caget(fitResultPVs)

#compute initial theta (rotation angle)
    maxValue = data.max()
    data0 = data/maxValue
    total = data0.sum()
    X, Y = indices(data0.shape)
    x = (X*data0).sum()/total
    y = (Y*data0).sum()/total
    col = data0[:, int(y)]
    row = data0[int(x), :]
    if col.size == row.size:
           rota = sqrt(abs((arange(col.size)-x)*(arange(row.size)-y)*data0).sum()/total)
           print "rota in square image: %f"%rota
    else:
           smin = min(col.size,row.size)
           data1=data0[0:smin,0:smin]
           rota = sqrt(abs((arange(smin)-x)*(arange(smin)-y)*data1).sum()/total)
           print "rota: %f"%rota
    rota1 = int(rota/(2*pi))
    rota = rota-rota1*2*pi
    print "fraction of rota: %f"%rota

    return height, offset, centroidX, centroidY, sigmaX, sigmaY,rota

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y, rota)
Returns (height, offset, centroidX, contrenoidY, sigmaX, sigmaY, rota)
    the gaussian parameters of a 2D distribution found by a fit"""
    initParams = moments(data)
#    params = [ 4202.,    478,   605,    10,    100,     0.]
    print 'initial params:'
    print initParams
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape))-data)
    p, success = leastsq(errorfunction, initParams)
    return p

def callback(value):
    try:
        #DA=caget('BTS-BI:BD1{VF:2}image1:ArrayData')
        (W, H) = caget(['%simage1:ArraySize0_RBV'%(cam), '%simage1:ArraySize1_RBV'%(cam)])
        D=+value[0:W*H]
        A=D.reshape((H, W))
        aa = A
        amax=aa.max()
        amin=aa.min()
        aa = aa - amin
        data = aa
        print data.size,data.shape

        fitParams = fitgaussian(data)
        print 'fitted params:'
        print fitParams
        print "***\n"
        caput("%sX-2DGauss:Sigma-I"%(cam), fitParams[4])
        caput("%sY-2DGauss:Sigma-I"%(cam), fitParams[5])
        caput("%sX-2DGauss:Center-I"%(cam), fitParams[2])
        caput("%sX-2DGauss:Center-I"%(cam), fitParams[3])
        caput("%s2DGauss:Tilt-I"%(cam), fitParams[6])

    except:
	traceback.print_exc()
        return

def main():
    #data=caget('%simage1:ArrayData'%(cam))
    #print data
    camonitor("%simage1:ArrayData"%(cam), callback, format=FORMAT_TIME)
    WaitForQuit()

if __name__ == "__main__":
    main()

