from neuron import h
h.load_file('stdlib.hoc')
h.load_file('stdrun.hoc')
from neuroAlter import prepCell
from sv import bproceed
import multiprocessing as mp
import numpy as np
from n128 import *
import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot, ticker

def getK(v1,v2,v,nstep,nS):
    K = np.empty(nstep)
    vp = np.empty((nstep,nS))
    for i in xrange(nstep):
        dv = v[i,:]-v1[i,:]-v2[i,:] 
        A = v1[i,:]*v2[i,:]
        At = np.transpose(A)
        K[i] = np.dot(At,dv)/np.dot(At,A) 
        if np.isnan(K[i]):
            K[i] = 0
        vp[i,:] = v1[i,:]+v2[i,:]+A*K[i]
    return K, vp

if __name__ == '__main__':
    sel = [0,1]
    if sel[0] == sel[1]:
        spikeTrain = [np.array([0,dt])]
    else:
        spikeTrain = [np.array([0]),np.array([dt])]

    # Bilinear
    bv = np.empty((run_nt,nv))
    bv0 = np.empty((run_nt,nv))

    for i in xrange(nv):
        v0 = vRange[i]
        fign = 'bi_'+str(v0)+'-'+str(sel[0])+str(sel[1])
        v, fired, _ = bproceed(cell, v0, synList, gList, vSL, spikeTrain, n, sel, run_t, tstep, fign, pos, loc, alphaR)
        bv[:,i] = v - v0

    for i in xrange(nv):
        v0 = vRange[i]
        cell.set_volt(v0)
        fign = 'bi0_'+str(v0)+'-'+str(sel[0])+str(sel[1])
        v, fired, _ = bproceed(cell, v0, synList, gList, vSL, spikeTrain, n, sel, run_t, tstep, fign, pos, loc, alphaR, cpi)
        bv0[:,i] = v - leakyDendV0[:,i]
