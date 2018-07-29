from neuron import h
h.load_file('stdlib.hoc')
h.load_file('stdrun.hoc')
from neuroAlter import prepCell
from sv import bproceed, write_one, gather_and_distribute_results
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

def generateBvData():

