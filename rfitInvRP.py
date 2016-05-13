#!/usr/bin/python

################################
#
#   Sript to calculate the r fit of the iverted Rotne Prager for a specific combination of a and p
#   THIS SCRIPT RETURNS A RESULT IN UNITS OF THE REDUCED VARIABLES s = 4*r/(a+p) AND lambda = a/p
#
################################

import numpy as np
import os
import numpy.linalg as LA
import sys
from scipy.optimize import curve_fit

def RotnePrager(rvec, p, a):
    _pradius = p/2.
    _polyrad = a/2.
    asq = (_polyrad * _polyrad + _pradius * _pradius)/2
    I = np.identity(3)
    r = LA.norm(rvec)
    rsq = r*r
    return _pradius * 3. / ( 4 * r ) * ( ( 1. + 2. * asq / (3. * rsq )) * I + ( 1. - 2.* asq / rsq ) * np.outer(rvec,rvec) / rsq )

def getInvRP(rvec,p,a,i=0,j=0,full=False):
    RP2p = np.identity(6)
    RP2p[0:3,3:] = RotnePrager(rvec, p, a=a)
    RP2p[3:,0:3] = RP2p[0:3,3:]
    RP2p[3:,3:] = np.identity(3) * p/a
    if full:
        return (LA.inv(RP2p)[:3,:3])
    return (LA.inv(RP2p)[:3,:3])[i,j]

def getPrefactors(rij,p,a):
    rr = np.outer(rij,rij) / rij.dot(rij)
    pref_rr = getInvRP(rij,p,a,0,1) / rr[0,1]   # This is correct!
    pref_I = (getInvRP(rij,p,a,0,0) - pref_rr * rr[0,0])
    return pref_I, pref_rr

def fpoly(x, a, b, c, d, e, f, g, h):
    x = np.asarray(x)
    return a + b*x**-1 + c*x**-2 + d*x**-3 + e*x**-4 + f*x**-5 + g*x**-6 + h*x**-7

def fpoly2(x, a, b, c, d, e, f, g, h):
    x = np.asarray(x)
    return a + b*x**-2 + c*x**-4 + d*x**-6 + e*x**-8 + f*x**-10 + g*x**-12 + h*x**-14

# fpoly = np.vectorize(fpoly)

def fitRPinv(a,p):
    # First, we need to find the prefactors pI and prr of the identity Matrix I and the outer product rr
    # for different sphere distances r.
    pI_store = []
    prr_store = []
    s_store = []
    r_store = []
    gets = lambda r: 4.*r/(a+p)
    lam = a/p
    #TODO change this so simply an r array or something, since the results only depend on r not on the vector
    for x in np.arange(0.001,15*(a+p)/2,0.025*(a+p)):
        vec = np.array([x,0.001,(a+p)/2])
        r = LA.norm(vec)
        s = gets(r)
        pI,prr = getPrefactors(vec,a=a,p=p)
        pI_store.append(pI)
        prr_store.append(prr)
        r_store.append(r)
        s_store.append(s)
    # Then we need to fit the result. For now only a polyfit.
    fitI, pcov = curve_fit(fpoly, s_store, pI_store)
    fitrr, pcov = curve_fit(fpoly, s_store, prr_store)
    return lam, fitI, fitrr

# This code writes r fit result to a file.
p = sys.argv[1] # Store p and a as strings
a = sys.argv[2]
directory = 'fits/'
filename = 'fitp' + p + 'a' + a + '.txt'
if not os.path.exists(directory):
    os.makedirs(directory)
Fitfile = open(directory + filename,'w')
lam, rfitpI, rfitprr = fitRPinv(a=float(a),p=float(p))
for i in range(len(rfitpI)):
    Fitfile.write( str(rfitpI[i]) + ' ' + str(rfitprr[i]) + '\n' )
Fitfile.close()
